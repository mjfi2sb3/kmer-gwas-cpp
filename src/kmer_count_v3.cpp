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
//     Previously the program held roughly three copies of every read and
//     accumulated k-mers without any cap, so peak memory tracked input size
//     directly rather than being something you could set. Tasks #4 and #17.
//
//  2. INODES ARE O(1) PER ACCESSION.
//     One pack file, plus one spill file only if the budget was exceeded.
//     Previously: 1 + num_bins directories + 2 x num_bins files, measured at
//     4,502 peak inodes for 1500 bins versus 2 now. That coupling is why bin
//     count could not be raised for large genomes without an inode explosion.
// ===========================================================================

#include "thread_pool.hpp"
#include "mem_limit.hpp"
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <atomic>
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
// Gzip decompression. zlib-ng's inflate is ~3x faster than stock zlib on FASTQ
// (measured 23.2 s -> 7.4 s per rice mate) and is API-shaped the same, so it is
// used when available (-DHAVE_ZLIBNG) and stock zlib is the portable fallback so
// the non-container profiles still build on clusters without it.
#ifdef HAVE_ZLIBNG
  #define WITH_GZFILEOP           // exposes zlib-ng's gz* file ops in the header
  #include <zlib-ng.h>
  #define KGZ_open  zng_gzopen
  #define KGZ_read  zng_gzread
  #define KGZ_close zng_gzclose
#else
  #include <zlib.h>
  #define KGZ_open  gzopen
  #define KGZ_read  gzread
  #define KGZ_close gzclose
#endif
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
   GzipStreambuf(const string& path) : gz_(KGZ_open(path.c_str(), "rb")) {}
   ~GzipStreambuf() { if (gz_) KGZ_close(gz_); }
   int underflow() override {
      int n = KGZ_read(gz_, buf_, sizeof(buf_));
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
   // then parsed with the FASTQ rule, silently yielding garbage.
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
// task — one of three redundant copies since removed.
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
   if (argc < 6 || argc > 9)
   {
      cout << "usage: " << argv[0]
           << " <accession> <num_bins> <output folder path> <R1> <R2> [memory_budget_GB] [min_count] [read_threads]\n"
           << "\n"
           << "  Writes <output>/<accession>.kbin — one self-indexed pack file holding\n"
           << "  this accession's canonical " << k << "-mers for every bin.\n"
           << "\n"
           << "  memory_budget_GB bounds k-mer accumulation. 0 (default) auto-sizes it\n"
           << "  from the memory actually enforced on this process (cgroup / SLURM /\n"
           << "  MemTotal); a positive value is honoured but capped at that limit. Peak\n"
           << "  memory is the budget regardless of genome size; exceeding it spills to a\n"
           << "  single temporary file rather than failing.\n"
           << "\n"
           << "  min_count drops k-mers seen fewer than this many times in this accession\n"
           << "  (default " << DEFAULT_MIN_COUNT << "). Low counts are overwhelmingly sequencing errors;\n"
           << "  raising it shrinks Stage 1 output and everything downstream.\n"
           << "\n"
           << "  read_threads decompresses the mates concurrently (default 2, one per file).\n"
           << "  Set 1 to serialise, e.g. on a single spinning disk where concurrent reads\n"
           << "  seek-thrash; on NVMe / parallel filesystems 2 roughly halves the read phase.\n";
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

      // Accumulation budget. A positive 6th argument is an explicit request,
      // honoured but capped at the memory actually enforced on this process; 0
      // (the default the pipeline passes) means AUTO — a fraction of that
      // enforced limit, so the budget tracks the real node whether it is a
      // shared allocation, an --mem=0 exclusive node, or a bare workstation.
      // See mem_limit.hpp. Output is identical whatever the budget; it only
      // changes how often accumulation spills.
      static const double BUDGET_FRACTION = 0.7;   // ~30% left for read batches + worker maps
      double budget_req_gb = (argc >= 7) ? atof(argv[6]) : 0.0;
      kmer::MemLimit lim   = kmer::detect_memory_limit();
      double limit_gb      = lim.bytes ? lim.bytes / 1e9 : 0.0;
      double budget_gb;
      string budget_note;
      if (budget_req_gb > 0.0) {
         if (limit_gb > 0.0 && budget_req_gb > limit_gb) {
            budget_gb   = limit_gb;
            budget_note = "explicit request reduced to the enforced limit";
         } else {
            budget_gb   = budget_req_gb;
            budget_note = "explicit request";
         }
      } else if (limit_gb > 0.0) {
         budget_gb   = BUDGET_FRACTION * limit_gb;
         budget_note = "auto (70% of enforced limit)";
      } else {
         budget_gb   = DEFAULT_BUDGET_GB;
         budget_note = "fallback default (no limit detected)";
      }
      if (budget_gb < 1.0) budget_gb = 1.0;

      kmer::Count min_count = (argc >= 8) ? (kmer::Count)atoi(argv[7]) : DEFAULT_MIN_COUNT;
      if (min_count < 1) min_count = 1;

      int read_threads = (argc >= 9) ? atoi(argv[8]) : 2;
      if (read_threads < 1) read_threads = 1;

      cout << "***************************** " << endl;
      cout << "PROCESSING ACCESSION: " << accession << endl;
      cout << "  k = " << k << ", bins = " << NUM_FILES
           << ", record = " << kmer::RECORD_BYTES << " B"
           << ", min_count = " << min_count << endl;
      cout << "  memory: enforced limit = "
           << (limit_gb > 0.0 ? std::to_string((long)limit_gb) + " GB" : string("unknown"))
           << " [" << lim.source << "]" << endl;
      cout << "  memory: budget = " << (long)budget_gb << " GB (" << budget_note << ")" << endl;
      auto start = chrono::steady_clock::now();

      filesystem::create_directories(output_path);
      const string pack_path  = output_path + accession + ".kbin";
      const string spill_path = output_path + accession + ".spill";

      // Bounded-memory accumulator. Replaces the 1 + B dirs + 2B files scheme
      // and the separate dedup pass entirely.
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
      // Decompress + parse one file, dispatching its read batches to the pool.
      // Safe to run concurrently for the two mates: push_task locks the queue,
      // BatchGate is mutex-guarded, and BinStore::append takes per-bin locks, so
      // nothing is shared unsynchronised. Read order does not matter — counts are
      // summed wherever a k-mer occurs.
      atomic<size_t> total_reads{0};
      auto drain_file = [&](const string &accession_path)
      {
         SequenceReader reader(accession_path);
         cout << "+++ PROCESSING " << accession_path
              << " (" << reader.format_name() << ")" << endl;
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
      };

      // Decompression is the read-phase bottleneck (CPU-bound inflate), so the
      // mates are read on separate threads by default — one per file — which on
      // NVMe / parallel HPC filesystems roughly halves the read phase. On a single
      // spinning disk concurrent streams can seek-thrash; --kmer_count_read_threads 1
      // serialises them. I/O is a small fraction of the inflate cost on the target.
      if (read_threads >= 2 && accession_list.size() >= 2)
      {
         vector<thread> readers;
         for (auto &p : accession_list) readers.emplace_back(drain_file, std::cref(p));
         for (auto &t : readers) t.join();
      }
      else
      {
         for (auto &p : accession_list) drain_file(p);
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

         // Phase A — sort / compact / filter every bin IN PARALLEL. The bins are
         // independent (finalize_bin_inplace touches only its own bin), and the
         // per-bin sorts dominate this stage, so we run them across the same pool
         // that did the counting instead of one bin at a time. Sorting is done in
         // place, so this adds no memory beyond what accumulation already held.
         for (size_t b = 0; b < NUM_FILES; b++)
            pool.push_task([&store, b, min_count]() { store.finalize_bin_inplace(b, min_count); });
         pool.wait_for_tasks();

         // Phase B — write bins in ascending order (the pack's offset table
         // requires it). Each slice is already finalised, so this is just I/O and
         // stays serial; free each bin's memory as soon as it is written.
         size_t total_kmers = 0;
         for (size_t b = 0; b < NUM_FILES; b++)
         {
            const vector<kmer::Record>& slice = store.bin_slice(b);
            pack.write_bin((uint32_t)b, slice);
            total_kmers += slice.size();
            store.release_bin(b);
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

