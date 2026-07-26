// ===========================================================================
//  kmer_count — Stage 1 of the k-mer GWAS pipeline
// ===========================================================================
//
//  Reads the paired FASTQ/FASTA files of ONE accession and writes that
//  accession's canonical k-mers, with per-accession counts, as a single
//  self-indexed pack file: <output>/<accession>.kbin
//
//  Pipeline inside this program:
//
//    SequenceReader -> batches of reads -> thread pool -> BinStore -> pack
//       (streaming)       (bounded by       (count       (bounded    (one
//                          BatchGate)        k-mers)      memory)     file)
//
//  TWO INVARIANTS TO PRESERVE WHEN EDITING:
//
//  1. MEMORY IS BOUNDED, NOT DATA-PROPORTIONAL.
//     Reads are streamed and discarded batch by batch; k-mer accumulation is
//     capped by a budget and spills to one temporary file if exceeded. Peak
//     memory therefore does not scale with genome size or sequencing depth.
//     Previously the program held ~3 copies of every read (coverage x genome
//     x 3 bytes: ~34 GB for rice, ~1.5 TB for wheat) and accumulated k-mers
//     without any cap, which made large genomes impossible. Tasks #4 and #17.
//
//  2. INODES ARE O(1) PER ACCESSION.
//     One pack file, plus one spill file only if the budget was exceeded.
//     Previously: 1 + num_bins directories + 2 x num_bins files, measured at
//     4,502 peak inodes for 1500 bins versus 2 now. That coupling is why bin
//     count could not be raised for large genomes without an inode explosion.
// ===========================================================================

#include "thread_pool.hpp"
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <string>
#include <cstring>
#include <chrono>
#include <memory>
#include <thread>
#include "kmer_key.hpp"
#include "pack_io.hpp"
#include "bin_store.hpp"
#include <zlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <filesystem>

using namespace std;

// Report a fatal I/O error and exit. Previously lived in mmap_io.cpp, which was
// otherwise entirely dead code (its get_reads/read_vector/write_vector were
// unused by both programs).
static void handle_error(const string &msg)
{
   perror(msg.c_str());
   exit(-1);
}

const int k = kmer::K;

using Key = kmer::Key;

// Reads handed to a worker in one go.
//
// This is a throughput/memory trade-off, and the direction is not the obvious
// one. LARGER batches are both faster and *lighter*, because a worker's hash
// map collapses duplicate k-mers before they are spilled: a bigger batch dedups
// more, so fewer records reach disk and far fewer short-lived maps are
// allocated. Measured on a 200 kb genome at high coverage (1.5M read pairs):
//
//     batch  1024 -> 21.5 s, 5.14 GB      batch 16384 ->  4.2 s, 1.69 GB
//     batch  2048 -> 18.2 s, 5.02 GB      batch 32768 ->  2.5 s, 1.15 GB
//     batch  4096 -> 14.8 s, 3.36 GB      batch 65536 ->  1.7 s, 1.00 GB
//
// The ceiling is genome complexity, which that benchmark understates. On a real
// genome nearly every k-mer in a batch is distinct, so a worker's map holds
// about (BATCH_READS x kmers_per_read) entries at ~60 B each — at 16384 reads of
// 150 bp that is ~98 MB per in-flight task. 16384 keeps the worst case bounded
// while capturing most of the speed-up.
//
// Task #17 replaces this heuristic with a real global memory budget, which is
// the correct fix; until then this constant is the safety margin.
static const size_t BATCH_READS = 16384;

// Default k-mer accumulation budget, overridable as the 6th argument. Peak
// memory is this figure regardless of genome size — exceeding it spills to a
// single temporary file rather than failing.
static const double DEFAULT_BUDGET_GB = 32.0;

// Default minimum per-accession count for a k-mer to be kept, overridable as
// the 7th argument.
//
// The previous code hardcoded "drop count == 1", i.e. this value fixed at 2. At
// typical sequencing depth, count-2 and count-3 k-mers are still overwhelmingly
// sequencing errors, and this method family normally filters at 2-5. Raising it
// compounds through everything downstream: Stage 1 output size, Stage 2 read
// volume, matrix_merge peak memory, and the union row count of the matrix.
//
// It is left at 2 by default so behaviour is unchanged unless asked for.
static const kmer::Count DEFAULT_MIN_COUNT = 2;

struct membuf : streambuf
{
   membuf(char *begin, char *end)
   {
      this->setg(begin, begin, end);
   }
};

class GzipStreambuf : public streambuf
{
   gzFile gz_;
   char buf_[262144];  // 256 KB inflate buffer
public:
   GzipStreambuf(const string& path) : gz_(gzopen(path.c_str(), "rb")) {}
   ~GzipStreambuf() { if (gz_) gzclose(gz_); }
   int underflow() override {
      int n = gzread(gz_, buf_, sizeof(buf_));
      if (n <= 0) return EOF;
      setg(buf_, buf_, buf_ + n);
      return (unsigned char)buf_[0];
   }
};

// ---------------------------------------------------------------------------
// SequenceReader — streams sequences out of a FASTQ or FASTA file.
//
// Handles plain and gzip-compressed input transparently (gzip is detected by
// the 0x1f 0x8b magic bytes, so the file extension is irrelevant), and both
// FASTQ and multi-line FASTA.
//
// The file is never loaded in full: next_batch() pulls at most `max_reads`
// sequences and returns, so the caller controls how much is resident.
// ---------------------------------------------------------------------------
class SequenceReader
{
public:
   enum Format { FASTQ, FASTA };

   explicit SequenceReader(const string &path) : path_(path)
   {
      // gzip magic bytes — independent of the file extension
      unsigned char m[2] = {0, 0};
      if (FILE *mf = fopen(path.c_str(), "rb")) {
         if (fread(m, 1, 2, mf) != 2) { /* short file; falls through to error below */ }
         fclose(mf);
      } else {
         handle_error("open " + path);
      }
      is_gz_ = (m[0] == 0x1f && m[1] == 0x8b);

      if (is_gz_) {
         gzbuf_ = make_unique<GzipStreambuf>(path);
         in_    = make_unique<istream>(gzbuf_.get());
      } else {
         // mmap the plain file and wrap it in an istream. The mapping stays
         // alive for the reader's lifetime and is released in the destructor.
         fd_ = open(path.c_str(), O_RDONLY);
         if (fd_ == -1) handle_error("open " + path);
         struct stat sb;
         if (fstat(fd_, &sb) == -1) handle_error("fstat " + path);
         size_ = sb.st_size;
         if (size_ == 0) handle_error("empty input file: " + path);
         map_ = static_cast<char *>(mmap(NULL, size_, PROT_READ, MAP_PRIVATE, fd_, 0u));
         if (map_ == MAP_FAILED) handle_error("mmap " + path);
         membuf_ = make_unique<membuf>(map_, map_ + size_);
         in_     = make_unique<istream>(membuf_.get());
      }

      // Peek at the first non-empty line to establish the format, then keep it:
      // for FASTA it is a header we must not lose, for FASTQ it is discarded.
      string first;
      while (getline(*in_, first) && first.empty()) {}
      if      (!first.empty() && first[0] == '>') format_ = FASTA;
      else if (!first.empty() && first[0] == '@') format_ = FASTQ;
      else handle_error("unrecognized file format (expected FASTA '>' or FASTQ '@'): " + path);

      if (format_ == FASTA) pending_header_ = true;   // 'first' was a record header
      line_no_ = 1;
   }

   ~SequenceReader()
   {
      in_.reset(); membuf_.reset(); gzbuf_.reset();
      if (map_ && map_ != MAP_FAILED) munmap(map_, size_);
      if (fd_ != -1) close(fd_);
   }

   Format format() const { return format_; }
   const char *format_name() const { return format_ == FASTQ ? "FASTQ" : "FASTA"; }

   // Fill `out` with up to max_reads sequences. Returns false once the input is
   // exhausted and nothing was produced.
   bool next_batch(vector<string> &out, size_t max_reads)
   {
      out.clear();
      out.reserve(max_reads);
      return format_ == FASTQ ? next_fastq(out, max_reads)
                              : next_fasta(out, max_reads);
   }

   size_t reads_seen() const { return reads_seen_; }

private:
   // FASTQ: fixed 4-line records; the sequence is the 2nd line of each.
   // line_no_ counts lines already consumed, with line 0 being the first header,
   // so sequence lines are exactly those where line_no_ % 4 == 1.
   bool next_fastq(vector<string> &out, size_t max_reads)
   {
      string line;
      while (out.size() < max_reads && getline(*in_, line)) {
         if (line_no_ % 4 == 1) { out.push_back(line); reads_seen_++; }
         line_no_++;
      }
      return !out.empty();
   }

   // FASTA: a header line starting with '>' followed by one or more sequence
   // lines, which must be concatenated. Previously this format was detected and
   // then parsed with the FASTQ rule, silently yielding garbage (task #7).
   bool next_fasta(vector<string> &out, size_t max_reads)
   {
      string line;
      while (out.size() < max_reads) {
         if (!pending_header_ && !getline(*in_, line)) break;   // EOF
         if (!pending_header_ && (line.empty() || line[0] != '>')) continue;
         pending_header_ = false;

         // Accumulate sequence lines until the next header or EOF.
         seq_.clear();
         while (in_->peek() != EOF && in_->peek() != '>') {
            if (!getline(*in_, line)) break;
            seq_ += line;
         }
         if (!seq_.empty()) { out.push_back(seq_); reads_seen_++; }
         if (in_->peek() == EOF) break;
      }
      return !out.empty();
   }

   string                    path_;
   bool                      is_gz_  = false;
   int                       fd_     = -1;
   char                     *map_    = nullptr;
   size_t                    size_   = 0;
   unique_ptr<GzipStreambuf> gzbuf_;
   unique_ptr<membuf>        membuf_;
   unique_ptr<istream>       in_;
   Format                    format_ = FASTQ;
   size_t                    line_no_ = 0;
   size_t                    reads_seen_ = 0;
   bool                      pending_header_ = false;
   string                    seq_;
};

// ---------------------------------------------------------------------------
// BatchGate — backpressure for the reader.
//
// thread_pool::push_task() queues without limit, so a reader that ran ahead of
// the workers would rebuild exactly the unbounded read buffer this design is
// meant to remove. The gate blocks the reader once MAX_INFLIGHT batches are
// queued or executing, capping resident read memory at
//   MAX_INFLIGHT x BATCH_READS x average_read_length.
// ---------------------------------------------------------------------------
class BatchGate
{
public:
   explicit BatchGate(size_t max_inflight) : max_(max_inflight) {}

   void acquire()
   {
      unique_lock<mutex> lk(m_);
      cv_.wait(lk, [&] { return inflight_ < max_; });
      inflight_++;
   }
   void release()
   {
      { lock_guard<mutex> lk(m_); inflight_--; }
      cv_.notify_one();
   }

private:
   mutex              m_;
   condition_variable cv_;
   size_t             inflight_ = 0;
   size_t             max_;
};

// ---------------------------------------------------------------------------
// kmers_obj — counts the k-mers of one batch of reads and hands them to the
// BinStore.
//
// The batch is MOVED in, not copied. It previously took `vector<string>&` but
// stored a by-value member, silently duplicating every read once per in-flight
// task — one of the three copies removed in task #4.
// ---------------------------------------------------------------------------
class kmers_obj
{
public:
   unordered_map<Key, kmer::Count, kmer::KeyHash> kmers_;
   vector<string> chunk;

   explicit kmers_obj(vector<string> &&batch) : chunk(std::move(batch)) {}

   void index(kmer::BinStore &store, size_t num_files)
   {
      // Rolling canonical encode: the forward and reverse-complement keys are
      // updated in O(1) per base instead of re-encoding each overlapping k-mer.
      //
      // Windows containing a non-ACGT base are SKIPPED, not counted.
      //
      // They used to be folded into the all-zero key (canonical() returned ""
      // and bit_encode("") gave 0) and discarded later. That silently destroyed
      // a real k-mer: A encodes as 00, so the genuine 51-mer AAAA...A is also
      // all-zero and was indistinguishable from the "invalid" marker. Poly-A is
      // its own canonical form (its reverse complement TTT...T sorts higher), so
      // it could not escape by canonicalising to something else — it was simply
      // lost, without warning, wherever a poly-A tract occurred. Task #6.
      //
      // Skipping is also less work: a read with one N no longer counts all k
      // spanning windows into a bucket that is thrown away.
      const uint8_t *T = kmer::base_table();
      for (auto &seq : chunk)
      {
         if (seq.size() < (size_t)k) continue;
         Key f, r;
         int valid = 0;
         for (size_t i = 0; i < seq.size(); i++)
         {
            uint8_t v = T[(unsigned char)seq[i]];
            if (v == 0xFF) { valid = 0; f = Key{}; r = Key{}; continue; }
            f.push_fwd(v); r.push_rev(v);
            if (valid < k) valid++;

            if (i + 1 >= (size_t)k && valid >= k)    // a fully-valid window ends here
               kmer::sat_inc(kmers_[(r < f) ? r : f]);
         }
      }

      // Hand the batch's counts to the store, grouped by bin.
      //
      // Sorting by bin first means ONE lock acquisition per (batch, bin) instead
      // of one per record. The previous code took a mutex for every single entry
      // and did an ofstream write while holding it.
      staging_.clear();
      staging_.reserve(kmers_.size());
      for (auto &pair_ : kmers_)
         staging_.push_back({(uint32_t)(kmer::legacy_bin_hash(pair_.first) % num_files),
                             kmer::Record{pair_.first, pair_.second}});

      std::sort(staging_.begin(), staging_.end(),
                [](const Staged &a, const Staged &b) { return a.bin < b.bin; });

      for (size_t i = 0; i < staging_.size(); )
      {
         size_t j = i;
         uint32_t bin = staging_[i].bin;
         while (j < staging_.size() && staging_[j].bin == bin) j++;
         run_.clear();
         run_.reserve(j - i);
         for (size_t t = i; t < j; t++) run_.push_back(staging_[t].rec);
         store.append(bin, run_.data(), run_.size());
         i = j;
      }
   }

private:
   struct Staged { uint32_t bin; kmer::Record rec; };
   vector<Staged>       staging_;
   vector<kmer::Record> run_;
};

vector<string> split(const string &str, const char &sep)
{
   string segment;
   vector<string> ret;
   istringstream ss(str);
   while (getline(ss, segment, sep))
      ret.push_back(segment);
   return ret;
}

int main(int argc, char *argv[])
{
   if (argc < 6 || argc > 8)
   {
      cout << "usage: " << argv[0]
           << " <accession> <num_bins> <output folder path> <R1> <R2> [memory_budget_GB] [min_count]\n"
           << "\n"
           << "  Writes <output>/<accession>.kbin — one self-indexed pack file holding\n"
           << "  this accession's canonical " << k << "-mers for every bin.\n"
           << "\n"
           << "  memory_budget_GB bounds k-mer accumulation (default "
           << DEFAULT_BUDGET_GB << "). Peak memory is\n"
           << "  the budget regardless of genome size; exceeding it spills to a single\n"
           << "  temporary file rather than failing.\n"
           << "\n"
           << "  min_count drops k-mers seen fewer than this many times in this accession\n"
           << "  (default " << DEFAULT_MIN_COUNT << "). Low counts are overwhelmingly sequencing errors;\n"
           << "  raising it shrinks Stage 1 output and everything downstream.\n";
      return -1;
   }

   try
   {
      string accession = argv[1];
      size_t NUM_FILES = atoi(argv[2]);
      string output_path = argv[3];
      if (output_path[output_path.size()-1] != '/')
         output_path += '/';

      vector<string> accession_list = { argv[4], argv[5] };
      double budget_gb = (argc >= 7) ? atof(argv[6]) : DEFAULT_BUDGET_GB;
      if (budget_gb <= 0) budget_gb = DEFAULT_BUDGET_GB;
      kmer::Count min_count = (argc >= 8) ? (kmer::Count)atoi(argv[7]) : DEFAULT_MIN_COUNT;
      if (min_count < 1) min_count = 1;

      cout << "***************************** " << endl;
      cout << "PROCESSING ACCESSION: " << accession << endl;
      cout << "  k = " << k << ", bins = " << NUM_FILES
           << ", record = " << kmer::RECORD_BYTES << " B"
           << ", memory budget = " << budget_gb << " GB"
           << ", min_count = " << min_count << endl;
      auto start = chrono::steady_clock::now();

      filesystem::create_directories(output_path);
      const string pack_path  = output_path + accession + ".kbin";
      const string spill_path = output_path + accession + ".spill";

      // Bounded-memory accumulator. Replaces the 1 + B dirs + 2B files scheme
      // and the separate dedup pass entirely (task #17).
      kmer::BinStore store(NUM_FILES, (size_t)(budget_gb * 1e9), spill_path);

      thread_pool pool;

      // Cap on batches queued-or-running. Bounds BOTH resident read memory and
      // the number of live per-worker hash maps, which is the larger term on a
      // complex genome (~98 MB per task at BATCH_READS=16384, 150 bp reads).
      // 2x the worker count keeps every thread fed with one batch in reserve.
      BatchGate gate(max<size_t>(4, 2 * thread::hardware_concurrency()));

      cout << "Building Kmer index (streaming) ... This will take a while" << endl;

      // -- Stream both mates, dispatching batches as they are parsed ---------
      // Reads are consumed and discarded batch by batch; the full read set is
      // never resident. Batching does not affect the result: counts for a k-mer
      // are summed wherever it occurs, so how reads are grouped changes only how
      // often compaction runs, not the totals.
      size_t total_reads = 0;
      for (auto &accession_path : accession_list)
      {
         cout << "+++ PROCESSING " << accession_path << endl;
         SequenceReader reader(accession_path);
         cout << "    format: " << reader.format_name() << endl;

         vector<string> batch;
         while (reader.next_batch(batch, BATCH_READS))
         {
            gate.acquire();                       // blocks if too many in flight
            pool.push_task([batch = std::move(batch), &store, NUM_FILES, &gate]() mutable
                           {
                              kmers_obj ko(std::move(batch));
                              ko.index(store, NUM_FILES);
                              gate.release();
                           });
         }
         total_reads += reader.reads_seen();
         cout << "    num reads: " << reader.reads_seen() << endl;
         cout << "+++++++++++++++++++++" << endl;
      }

      pool.wait_for_tasks();

      auto end = chrono::steady_clock::now();
      cout << "Streaming & counting: "
           << chrono::duration_cast<chrono::seconds>(end - start).count()
           << " sec  (" << total_reads << " reads)" << endl;
      if (store.spilled())
         cout << "  (memory budget exceeded; spilled " << store.spill_rounds()
              << " block(s) to disk)" << endl;
      start = chrono::steady_clock::now();

      // -- Finalise: one sorted, deduplicated, filtered slice per bin --------
      // Bins are written in ascending order into a single self-indexed pack.
      cout << "Writing pack " << pack_path << " ..." << endl;
      {
         kmer::PackWriter pack(pack_path, (uint32_t)NUM_FILES);
         vector<kmer::Record> slice;
         size_t total_kmers = 0;
         for (size_t b = 0; b < NUM_FILES; b++)
         {
            store.finalize_bin(b, slice, min_count);
            pack.write_bin((uint32_t)b, slice);
            total_kmers += slice.size();
         }
         pack.finish();
         cout << "  " << total_kmers << " distinct k-mers across "
              << NUM_FILES << " bins" << endl;
      }
      store.close_spill();

      end = chrono::steady_clock::now();
      cout << "Finalising & packing: "
           << chrono::duration_cast<chrono::seconds>(end - start).count()
           << " sec" << endl;

      cout << "FINISHED ACCESSION: " << accession << endl;
      cout << "***************************** " << endl;
   }

   catch (exception const &e)
   {
      cout << e.what() << endl;
   }

   return 0;
}

