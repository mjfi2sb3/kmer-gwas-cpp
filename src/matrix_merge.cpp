//#include "thread_pool.hpp"
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>
#include <bitset>
#include <mutex>
#include <thread>
#include <atomic>
#include <cerrno>
#include <functional>
#include <string>
#include <cstring>
#include <chrono>
#include <stdexcept>
#include "kmer_key.hpp"
#include "pack_io.hpp"
#include <queue>
#include <memory>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <filesystem>
#include <getopt.h>
#include <unistd.h>


using namespace std;

const int k = kmer::K;

using Key = kmer::Key;

inline string bit_decode(const Key &code) { return kmer::decode(code); }

// const size_t NUM_ACC = 100;
vector<string> split(const string &str, const char &sep)
{
    string segment;
    vector<string> ret;
    istringstream ss(str);
    while (getline(ss, segment, sep))
        ret.push_back(segment);
    return ret;
}

/*vector<string> get_accessions(string accessions_path)
{
        vector<string> accessions;
	cout << "accessions path: " << accessions_path << endl;
	try {
        	ifstream stream(accessions_path);
        	string line;
        	while (getline(stream, line))
        	    accessions.push_back(line);
        	stream.close();
	}
    catch (exception const &e)
        	{cout << e.what() << endl;}
        	
	if (accessions.size() == 0) 
		cout << "could not read accessions, make sure file is not empty." << endl;
	return accessions;
}*/

vector<string> get_accessions(const string& accessions_path) {
    vector<string> accessions;
    cout << "Accessions path: " << accessions_path << endl;

    try {
        ifstream stream(accessions_path);
        if (!stream) {
            throw runtime_error("Failed to open accessions file: " + accessions_path);
        }

        string line;
        while (getline(stream, line)) {
            accessions.push_back(line);
        }

        stream.close();
    } catch (const exception& e) {
        cerr << e.what() << endl;
        throw;
    }

    if (accessions.empty()) {
        throw runtime_error("Could not read accessions. Make sure the file is not empty.");
    }

    return accessions;
}

// ===========================================================================
// Parallelising the merge by RANGE-SHARDING the key space.
//
// The merge below is a k-way streaming merge over N sorted per-accession slices.
// It parallelises cleanly because the key space can be cut into S contiguous
// shards that are INDEPENDENT: a shard's rows depend only on that shard's keys,
// and every accession's occurrence of a given key lands in the SAME shard.
//
// Why the same shard, deterministically: a shard is the top B bits of the key
// (S = 2^B). The first base of a k-mer sits in the key's HIGHEST bits and
// A<C<G<T maps to 0<1<2<3, so numeric key order == lexicographic k-mer order and
// the top B bits are the sort prefix. Two things follow:
//   1. keys sharing a shard form a CONTIGUOUS run in every pack's sorted bin, so
//      a shard is a byte range [lo,hi) found once per accession by a linear scan;
//   2. a k-mer encodes to the identical canonical Key at the identical k in every
//      accession (the footer pins k), so it maps to the identical shard in every
//      accession. Shard-s therefore sees ALL occurrences of every shard-s k-mer,
//      so its per-k-mer occurrence count (used by the MAF and core filters) is the
//      true panel-wide count, and NUM_ACC stays the GLOBAL accession count.
//
// Each shard runs an independent merge into its own part file; the parts are
// concatenated in shard order. Because shard 0 < shard 1 < ... in key order and
// no key straddles a shard, the concatenation is the same globally-sorted stream
// the single-threaded merge produced — byte-for-byte. S=1 collapses to exactly
// the original path. This is a performance change only; output is identical.
// ===========================================================================

// Top B bits of the 2k-bit key: the sort prefix, hence the shard id. B must be
// in [1, BITS-1]; B==0 means a single shard (id 0). Shifts are kept in range for
// both the one-word (k<=32) and two-word (k<=63) key forms.
static inline uint64_t key_top_bits(const Key& key, int B)
{
    if (B <= 0) return 0;
    const int sh = kmer::BITS - B;            // in [1, BITS-1]
#if KMER_K * 2 <= 64
    return key.w0 >> sh;                       // sh in [1,63]: well-defined
#else
    if (sh >= 64) return key.w1 >> (sh - 64);  // top bits live entirely in w1
    return (key.w1 << (64 - sh)) | (key.w0 >> sh);
#endif
}

// Run fn(i) for i in [0,n) across up to `threads` workers pulling off a shared
// counter. The first exception thrown in any worker is rethrown after all join,
// so a failing shard/accession surfaces as a normal error rather than terminate.
template <class F>
static void parallel_for(size_t n, uint threads, F&& fn)
{
    if (n == 0) return;
    uint t = (uint)min<size_t>(threads < 1 ? 1 : threads, n);
    atomic<size_t> next{0};
    exception_ptr err;
    mutex em;
    auto worker = [&]() {
        try {
            size_t i;
            while ((i = next.fetch_add(1)) < n) fn(i);
        } catch (...) {
            lock_guard<mutex> lk(em);
            if (!err) err = current_exception();
        }
    };
    if (t <= 1) {
        worker();
    } else {
        vector<thread> ths;
        ths.reserve(t);
        for (uint k = 0; k < t; k++) ths.emplace_back(worker);
        for (auto& x : ths) x.join();
    }
    if (err) rethrow_exception(err);
}

// Streams a byte range [start_off, start_off + nrec*RECORD_BYTES) of a sorted bin
// via pread() on a SHARED read-only fd. pread is positioned — it neither uses nor
// moves the fd's file offset — so many RangeCursors over the same fd (different
// shards of the same accession, on different threads) read concurrently and
// safely. One fd per accession is opened for the whole merge, N total, not N x S.
class RangeCursor {
public:
    // 2048 records = 36 KB/pread: large enough that syscall overhead is
    // negligible, small enough that the N-per-shard x T-concurrent-shards buffer
    // footprint stays modest (see the matrix_merge_memory note in nextflow.config).
    RangeCursor(int fd, uint64_t start_off, size_t nrec, size_t buf_records = 2048)
        : fd_(fd), next_off_(start_off), remaining_(nrec),
          buf_cap_(buf_records ? buf_records : 1)
    {
        buf_.resize(buf_cap_);
        refill();
    }

    bool valid() const { return pos_ < filled_; }
    const kmer::Record& current() const { return buf_[pos_]; }
    void advance() { if (++pos_ >= filled_) refill(); }

private:
    void refill()
    {
        pos_ = 0;
        filled_ = 0;
        if (!remaining_) return;
        size_t n    = remaining_ < buf_cap_ ? remaining_ : buf_cap_;
        size_t want = n * kmer::RECORD_BYTES;
        char*  dst  = reinterpret_cast<char*>(buf_.data());
        size_t got  = 0;
        // pread may return a short count without error; loop until the records
        // are fully read (or a real error / premature EOF is hit).
        while (got < want) {
            ssize_t r = pread(fd_, dst + got, want - got, (off_t)(next_off_ + got));
            if (r < 0) {
                if (errno == EINTR) continue;
                throw runtime_error("RangeCursor: pread failed");
            }
            if (r == 0) throw runtime_error("RangeCursor: unexpected EOF in bin slice");
            got += (size_t)r;
        }
        next_off_  += want;
        filled_     = n;
        remaining_ -= n;
    }

    int                 fd_;
    uint64_t            next_off_;
    size_t              remaining_;
    size_t              buf_cap_;
    vector<kmer::Record> buf_;
    size_t              pos_    = 0;
    size_t              filled_ = 0;
};

// ===========================================================================
// merge_chunk — build one bin's matrix by k-way merging the accessions' sorted
// pack slices, sharded across threads (see the range-sharding note above).
//
// MEMORY: O(number of accessions), NOT O(union k-mers x accessions).
//
// An earlier implementation loaded the whole bin into
// `unordered_map<Key, vector<Count>>`, needing union_rows x (N+1) x 2 bytes:
// two bytes of count per accession for every distinct k-mer in the bin. At tens
// of thousands of accessions a single bin's union is gigabytes, which is why
// num_bins used to grow with genome size, and num_bins is still what drives the
// inode count.
//
// Because Stage 1 emits each bin slice SORTED, the same result comes from a
// streaming k-way merge: advance one cursor per accession through a heap, pop
// the smallest key, gather the counts from every accession holding it, emit the
// row. The only resident state is a read buffer per accession (in each active
// shard worker) plus one reusable row per worker; nothing scales with the k-mer
// union.
//
// The dense row is still inherently O(N) to WRITE: two bytes per accession, so
// tens of KB per k-mer at large panels, however it was computed. That is the
// output format's cost, not the merge's.
// ===========================================================================
void merge_chunk(const uint file_index, const uint min_occur, string input_path,
                 string accessions_path, string delimiter, bool show_count,
                 string ouput_dir, bool write_core, bool bit_packed, uint num_threads = 1)
{
    auto accessions = get_accessions(accessions_path);
    size_t NUM_ACC = accessions.size();

    // The accession-occurrence counter is a 16-bit kmer::Count.
    if (NUM_ACC > kmer::COUNT_MAX)
        throw runtime_error("Too many accessions (" + to_string(NUM_ACC) + "): the "
                            "occurrence counter is 16-bit, maximum is "
                            + to_string(kmer::COUNT_MAX) + ". Widen kmer::Count to raise this.");

    // One descriptor per accession stays open for the whole merge. Raise the
    // soft limit to the hard limit so large panels work without the caller
    // having to remember `ulimit -n`.
    {
        struct rlimit rl;
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur < rl.rlim_max) {
            rl.rlim_cur = rl.rlim_max;
            setrlimit(RLIMIT_NOFILE, &rl);
        }
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur < NUM_ACC + 16)
            throw runtime_error("Need at least " + to_string(NUM_ACC + 16)
                                + " file descriptors for " + to_string(NUM_ACC)
                                + " accessions, but the limit is " + to_string(rl.rlim_cur)
                                + ". Raise it with `ulimit -n`.");
    }

    // Check every pack exists before opening any, so a missing accession fails
    // immediately rather than part-way through a long merge.
    vector<string> paths(NUM_ACC);
    for (size_t i = 0; i < NUM_ACC; i++) {
        paths[i] = input_path + accessions[i] + ".kbin";
        if (!filesystem::exists(paths[i]))
            throw runtime_error("Pack file does not exist: " + paths[i]);
    }

    // -- Open one read-only fd per accession and read its bin geometry once ---
    // The fd stays open for the whole merge and is shared by every shard worker
    // via pread (positioned, offset-free), so we need N descriptors, not N x S.
    // PackReader validates magic/version/k and gives the bin's byte offset and
    // record count without reading the payload.
    struct Acc { int fd; uint64_t off; size_t nrec; };
    vector<Acc> acc(NUM_ACC);
    for (size_t i = 0; i < NUM_ACC; i++) {
        kmer::PackReader r(paths[i]);
        acc[i].off  = r.bin_offset(file_index);
        acc[i].nrec = r.bin_records(file_index);
        acc[i].fd   = ::open(paths[i].c_str(), O_RDONLY);
        if (acc[i].fd < 0) {
            for (size_t j = 0; j < i; j++) ::close(acc[j].fd);
            throw runtime_error("Cannot open pack file: " + paths[i]);
        }
    }
    auto close_all = [&]() { for (auto& a : acc) if (a.fd >= 0) { ::close(a.fd); a.fd = -1; } };

    // -- Decide the shard count S = 2^B --------------------------------------
    // Over-partition beyond the thread count (~4x) so a dynamic queue can even
    // out the GC-skew in how many k-mers fall in each key range. S=1 reproduces
    // the single-threaded path exactly. Tiny bins are not worth splitting.
    size_t total_records = 0;
    for (auto& a : acc) total_records += a.nrec;

    uint S = 1;
    int  B = 0;
    const size_t SMALL_BIN = 1u << 16;   // below this, sharding overhead does not pay
    if (num_threads > 1 && total_records >= SMALL_BIN) {
        uint target = num_threads * 4;
        if (target > 128) target = 128;             // cap: bounds part files and B
        while ((1u << (B + 1)) <= target) B++;      // B = floor(log2(target)), >= 1
        if (B < 1) B = 1;
        if (B > kmer::BITS - 1) B = kmer::BITS - 1;
        S = 1u << B;
    }

    // Report the plan before the (potentially long) merge, so the effective
    // parallelism — and the reason when it collapses to single-threaded — is
    // visible in the task log up front, not only inferable from the end summary.
    // Only min(num_threads, S) workers are ever dispatched.
    uint eff_threads = (uint)min<size_t>(num_threads, S);
    if (S > 1) {
        cout << "  merge plan: " << S << " shards over " << eff_threads << " threads" << endl;
    } else if (num_threads <= 1) {
        cout << "  merge plan: 1 shard, single-threaded (--threads 1)" << endl;
    } else {
        cout << "  merge plan: 1 shard, single-threaded (input below sharding threshold, "
             << total_records << " records)" << endl;
    }

    // -- Find each accession's S+1 shard boundaries (record indices) ----------
    // A single linear scan of the sorted bin: top_bits(key) is monotonic
    // non-decreasing, so shard s occupies [bnd[s], bnd[s+1]). bnd[0]=0 and
    // bnd[S]=nrec always. Done in parallel across accessions; each reads its bin
    // once (transient, freed immediately), so peak extra memory is one bin per
    // worker, not the whole panel.
    vector<vector<size_t>> bnd(NUM_ACC);
    try {
        if (S == 1) {
            for (size_t i = 0; i < NUM_ACC; i++) bnd[i] = {0, acc[i].nrec};
        } else {
            parallel_for(NUM_ACC, num_threads, [&](size_t i) {
                size_t nrec = acc[i].nrec;
                vector<size_t> b(S + 1, nrec);
                b[0] = 0;
                if (nrec) {
                    vector<kmer::Record> recs(nrec);
                    // Sequential pread of the whole bin slice for this accession.
                    char*  dst = reinterpret_cast<char*>(recs.data());
                    size_t want = nrec * kmer::RECORD_BYTES, got = 0;
                    while (got < want) {
                        ssize_t r = pread(acc[i].fd, dst + got, want - got, (off_t)(acc[i].off + got));
                        if (r < 0) { if (errno == EINTR) continue; throw runtime_error("boundary pread failed on " + paths[i]); }
                        if (r == 0) throw runtime_error("boundary pread hit EOF on " + paths[i]);
                        got += (size_t)r;
                    }
                    size_t r = 0;
                    for (uint s = 1; s < S; s++) {
                        while (r < nrec && key_top_bits(recs[r].key, B) < s) r++;
                        b[s] = r;
                    }
                }
                bnd[i] = move(b);
            });
        }
    } catch (...) { close_all(); throw; }

    auto part_matrix = [&](uint s) { return ouput_dir + to_string(file_index) + "_matrix.tsv.p" + to_string(s); };
    auto part_core   = [&](uint s) { return ouput_dir + to_string(file_index) + "_core.txt.p"   + to_string(s); };
    const string final_matrix = ouput_dir + to_string(file_index) + "_matrix.tsv";
    const string final_core   = ouput_dir + to_string(file_index) + "_core.txt";

    atomic<size_t> total_rows{0};

    // -- One shard: an independent k-way merge over its key range -------------
    // Each shard owns its output part and its own row/packed/touched scratch, so
    // there is no cross-shard shared state. row and the heap idx are indexed by
    // the GLOBAL accession index i, so counts land in the right column and occ
    // (compared against the GLOBAL NUM_ACC) is the true panel-wide count.
    auto merge_shard = [&](uint s) {
        ofstream m_stream(part_matrix(s));
        if (!m_stream) throw runtime_error("Cannot create " + part_matrix(s));
        ofstream ck_stream;
        if (write_core) {
            ck_stream.open(part_core(s));
            if (!ck_stream) throw runtime_error("Cannot create " + part_core(s));
        }

        vector<kmer::Count>   row(NUM_ACC, 0);
        vector<unsigned char> packed;
        vector<uint32_t>      touched;
        touched.reserve(256);

        auto write_row = [&](const Key& key, size_t occ) {
            // Optionally record every k-mer present in ALL accessions to the
            // separate _core.txt file. This is independent of the matrix: it
            // neither adds nor removes a row here. Whether such k-mers appear in
            // the matrix is the MAF filter's decision below (via --threshold).
            if (write_core && occ == NUM_ACC)
                ck_stream << bit_decode(key) << "\n";
            // Two-sided MAF filter: discard both rare and near-ubiquitous k-mers
            // (both carry little association signal). min_occur == 0 disables it.
            if (occ < min_occur || occ > NUM_ACC - min_occur) return;

            // Row layout: <kmer> TAB <v0> D <v1> D ... D <vN-1>
            // The k-mer is always tab-separated from the values so the first column
            // can be split off even when --delimiter none packs the values
            // together; `delimiter` separates values from EACH OTHER only.
            m_stream << bit_decode(key) << "\t";
            if (bit_packed) {
                // One bit per accession, LSB-first within each byte, written as
                // lowercase hex: N accessions pack into ceil(N/8) bytes, i.e. about
                // N/4 hex characters, against roughly 2N characters for the delimited
                // text form — about 8x smaller. Presence/absence only, so there is
                // nothing to lose by packing.
                packed.assign((NUM_ACC + 7) / 8, 0);
                for (size_t i = 0; i < NUM_ACC; i++)
                    if (row[i]) packed[i >> 3] |= (unsigned char)(1u << (i & 7));
                static const char HEX[] = "0123456789abcdef";
                for (unsigned char b : packed) { m_stream << HEX[b >> 4] << HEX[b & 15]; }
            } else {
                for (size_t i = 0; i < NUM_ACC; i++) {
                    kmer::Count freq = row[i];
                    if (!show_count && freq > 0) freq = 1;
                    if (i) m_stream << delimiter;
                    m_stream << freq;
                }
            }
            m_stream << "\n";
        };

        // -- Open one range cursor per accession over this shard's slice -------
        vector<unique_ptr<RangeCursor>> cur(NUM_ACC);
        for (size_t i = 0; i < NUM_ACC; i++) {
            size_t lo = bnd[i][s], hi = bnd[i][s + 1];
            if (hi > lo)
                cur[i] = make_unique<RangeCursor>(acc[i].fd, acc[i].off + lo * kmer::RECORD_BYTES, hi - lo);
        }

        // -- Min-heap over the cursors' current keys ---------------------------
        // Greater-than comparison so priority_queue yields the SMALLEST key.
        struct HeapItem { Key key; uint32_t idx; };
        auto worse = [](const HeapItem& a, const HeapItem& b) {
            if (a.key == b.key) return a.idx > b.idx;
            return b.key < a.key;
        };
        priority_queue<HeapItem, vector<HeapItem>, decltype(worse)> heap(worse);

        for (size_t i = 0; i < NUM_ACC; i++)
            if (cur[i] && cur[i]->valid()) heap.push({cur[i]->current().key, (uint32_t)i});

        size_t rows = 0;
        while (!heap.empty())
        {
            Key key = heap.top().key;
            touched.clear();

            // Drain every accession positioned on this key.
            while (!heap.empty() && heap.top().key == key) {
                uint32_t i = heap.top().idx;
                heap.pop();
                row[i] = cur[i]->current().count;
                touched.push_back(i);
                cur[i]->advance();
                if (cur[i]->valid()) heap.push({cur[i]->current().key, i});
            }

            write_row(key, touched.size());
            rows++;
            for (uint32_t i : touched) row[i] = 0;    // clear only what was set
        }

        m_stream.close();
        if (write_core) ck_stream.close();
        total_rows.fetch_add(rows);
    };

    // -- Run the shards on the thread pool, then stitch the parts in order ----
    // Shards are dispatched off a shared counter (S >= threads), so a thread that
    // draws a light shard immediately picks up the next one.
    try {
        parallel_for(S, num_threads, merge_shard);
    } catch (...) { close_all(); throw; }

    close_all();

    // Concatenate the per-shard parts in shard order. Shard 0 < shard 1 < ... in
    // key order and no key straddles a shard, so this is the same globally-sorted
    // stream the single-threaded merge produced. Streamed (not slurped) because a
    // bin's matrix can be many GB; parts are removed once folded in.
    auto stitch = [&](const string& final_path, function<string(uint)> part_of) {
        ofstream out(final_path, ios::binary | ios::trunc);
        if (!out) throw runtime_error("Cannot create " + final_path);
        vector<char> buf(1u << 20);
        for (uint s = 0; s < S; s++) {
            ifstream in(part_of(s), ios::binary);
            if (!in) throw runtime_error("Cannot read part " + part_of(s));
            while (in) {
                in.read(buf.data(), (streamsize)buf.size());
                streamsize n = in.gcount();
                if (n > 0) out.write(buf.data(), n);
            }
            in.close();
            filesystem::remove(part_of(s));
        }
        out.close();
    };
    stitch(final_matrix, part_matrix);
    if (write_core) stitch(final_core, part_core);

    cout << "  merged " << total_rows.load() << " distinct k-mers across "
         << NUM_ACC << " accessions"
         << " (" << S << (S == 1 ? " shard" : " shards")
         << ", " << num_threads << (num_threads == 1 ? " thread)" : " threads)") << endl;
}

int main(int argc, char *argv[])
{
	string input_path, accessions_path;
	uint file_index = 0;
	bool have_index = false;      // --index is required; see the guard below
	uint min_occur = 0;
	std::string delimiter = "";  // Default: none (concatenated, presence/absence)
	bool show_count = false;  // Default is to show presence/absence
	bool write_core = false;  // Default is to skip core k-mers file
	bool bit_packed = false;  // --encoding bits: 1 bit per accession (else text)
	uint num_bins = 0;
	uint available_threads;
	if (const char* s = getenv("SLURM_CPUS_PER_TASK")) {
	    uint cpus = (uint)atoi(s);
	    available_threads = cpus > 1 ? cpus : max(1u, thread::hardware_concurrency());
	} else {
	    available_threads = max(1u, thread::hardware_concurrency());
	}
	uint num_threads = available_threads;

	// Define the long options
	static struct option long_options[] = {
		{"input",    required_argument, 0, 'i'},
		{"accessions", required_argument, 0, 'a'},
		{"index",    required_argument, 0, 'f'},
		{"threshold", required_argument, 0, 't'},
		{"delimiter", required_argument, 0, 'd'},
		{"encoding", required_argument, 0, 'e'},
		{"count",    required_argument, 0, 'c'},
		{"core",     required_argument, 0, 'r'},
		{"bins",     required_argument, 0, 'b'},
		{"threads",  required_argument, 0, 'T'},
		{0, 0, 0, 0}
	};
    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "i:a:f:t:d:e:c:r:b:T:", long_options, &option_index)) != -1)
    {
        switch (c) 
        {
            case 'i':
                input_path = optarg;
                break;
            case 'a':
                accessions_path = optarg;
                break;
            case 'f':
            	if (*optarg == '\0' || !all_of(optarg, optarg + strlen(optarg), ::isdigit)) {
					cerr << "Error: --index requires a valid non-negative integer argument." << endl;
					return -1;
				}
				file_index = (uint)stoul(optarg);
				have_index = true;
				break;
            case 't':
                min_occur = stoi(optarg);
                break;
            case '?':
                // getopt_long will print an error message
                break;
            case 'd':
				if (string(optarg) == "tab")
				{
				    delimiter = "\t";
				}
				else if (string(optarg) == "space")
				{
				    delimiter = " ";
				}
				else if (string(optarg) == "none")
				{
				    delimiter = "";
				}
				else if (string(optarg) == "bits")
				{
				    // 'bits' is an ENCODING, not a delimiter (it was under
				    // --delimiter in v3.0.0; that was a naming error). Reject it
				    // here so the mistake is caught rather than silently ignored.
				    cerr << "Error: 'bits' is not a delimiter — it is an output "
				            "encoding. Use --encoding bits (with --delimiter for "
				            "the text encoding only)." << endl;
				    return -1;
				}
				else
				{
				    cerr << "Invalid delimiter option. Use 'tab', 'space' or 'none'." << endl;
				    return -1;
				}
				break;
            case 'e':
				if (string(optarg) == "text")
				{
				    bit_packed = false;
				}
				else if (string(optarg) == "bits")
				{
				    // Bit-packed presence/absence: 1 bit per accession written as
				    // hex, instead of one delimited character. See write_row.
				    bit_packed = true;
				}
				else
				{
				    cerr << "Invalid encoding option. Use 'text' or 'bits'." << endl;
				    return -1;
				}
				break;
            case 'c':
				if (string(optarg) == "y") 
				{
				    show_count = true;
				} 
				else if (string(optarg) == "n") 
				{
				    show_count = false;
				} 
				else 
				{
				    cerr << "Invalid display option. Use 'y' or 'n'." << endl;
				    return -1;
				}
				break;
            case 'r':
				if (string(optarg) == "y")
				    write_core = true;
				else if (string(optarg) == "n")
				    write_core = false;
				else {
				    cerr << "Invalid core option. Use 'y' or 'n'." << endl;
				    return -1;
				}
				break;
            case 'b':
                num_bins = stoi(optarg);
                break;
            case 'T':
                num_threads = stoi(optarg);
                break;
            default:
                abort();
        }
    }
    // Cap threads at available resources
    if (num_threads > available_threads) {
        cerr << "Warning: --threads " << num_threads << " exceeds available CPUs ("
             << available_threads << "), capping at " << available_threads << "." << endl;
        num_threads = available_threads;
    }

    // The two content/encoding constraints, both for the same underlying reason:
    // a value vector is only recoverable if the values can be told apart.
    //
    //   --encoding bits   packs presence/absence into single bits, so counts
    //                     cannot be represented at all.
    //   --delimiter none  concatenates values with no separator, so multi-digit
    //                     counts run together and cannot be split.
    //
    // Both are therefore rejected in combination with --count y. --delimiter tab
    // and --delimiter space carry counts fine.
    if (bit_packed && show_count) {
        cerr << "Error: --encoding bits is presence/absence only and cannot be "
                "combined with --count y. Use --encoding text (with --delimiter "
                "tab or space) for counts." << endl;
        return -1;
    }
    if (delimiter.empty() && !bit_packed && show_count) {
        cerr << "Error: --delimiter none cannot be combined with --count y — "
                "multi-digit counts would be concatenated with no separator and "
                "the matrix could not be parsed. Use --delimiter tab or space for "
                "counts, or --count n for presence/absence." << endl;
        return -1;
    }

    // Check if all required options are provided.
    // --index must be included: it was previously read from an uninitialised
    // variable when omitted, silently merging an arbitrary bin.
    if (input_path.empty() || accessions_path.empty() || !have_index)
    /*{
        cout << "usage: " << argv[0] << " --input <input path> --accessions <accessions path> --index <file index> --threshold <min occurence threshold>\n";
        return -1;
    }*/
    {
		cout << "\nusage: " << argv[0] << "\n"
		     << "\t\t--input <input path> \n"
		     << "\t\t--accessions <accessions path> \n"
		     << "\t\t--index <file index which corresponds to bin> \n"
		     << "\t\t--threshold <value> (two-sided MAF filter, default: " << min_occur << ")\n"
		     << "\t\t            Keeps k-mers carried by [value, n_accessions - value] accessions.\n"
		     << "\t\t            Filters on accession occurrence, NOT on per-accession count. 0 = off.\n"
		     << "\t\t            Note: Use --threshold 20 or --threshold=20 (both work)\n"
		     << "\t\t--delimiter <tab|space|none> (default: " << (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none")) << ")\n"
		     << "\t\t            Separator between values in the text encoding.\n"
		     << "\t\t            'none' concatenates and is presence/absence only.\n"
		     << "\t\t--encoding <text|bits> (default: " << (bit_packed ? "bits" : "text") << ")\n"
		     << "\t\t            'bits' packs presence/absence 1 bit per accession as\n"
		     << "\t\t            hex, ~8x smaller than tab at large cohorts. Ignores\n"
		     << "\t\t            --delimiter and requires --count n.\n"
		     << "\t\t--count <print matrix as absence/presence or actual k-mer counts; type: y|n> (default: " << (show_count ? "y" : "n") << ")\n"
		     << "\t\t--core    <write k-mers present in ALL accessions to _core.txt; type: y|n> (default: n)\n"
		     << "\t\t            Independent of the matrix: it only adds the file. Whether such\n"
		     << "\t\t            k-mers appear in the matrix is set by --threshold, not this flag.\n"
		     << "\t\t--bins    <number of bins (used in output folder name)> (default: 0)\n"
		     << "\t\t--threads <parallel merge threads> (default: SLURM_CPUS_PER_TASK if set, else hardware concurrency)\n"
		     << "\t\t            The bin's key space is split into contiguous shards merged in\n"
		     << "\t\t            parallel; output is byte-identical to a single-threaded run.\n"
		     << "\n"
		     << "\tHow --encoding, --delimiter and --count combine:\n"
		     << "\n"
		     << "\t  encoding  count  delimiter   result\n"
		     << "\t  --------  -----  ---------   --------------------------------------\n"
		     << "\t  text      n      tab         KMER <tab> 1 <tab> 0 <tab> 1\n"
		     << "\t  text      n      space       KMER <tab> 1 0 1\n"
		     << "\t  text      n      none        KMER <tab> 101   (concatenated)\n"
		     << "\t  text      y      tab         KMER <tab> 5 <tab> 0 <tab> 3\n"
		     << "\t  text      y      space       KMER <tab> 5 0 3\n"
		     << "\t  text      y      none        REJECTED  (multi-digit counts merge)\n"
		     << "\t  bits      n      (ignored)   KMER <tab> 05    (1 bit/accession, hex)\n"
		     << "\t  bits      y      -           REJECTED  (bits is presence/absence only)\n"
		     << "\n"
		     << "\t  The k-mer is always tab-separated from the values; --delimiter\n"
		     << "\t  separates the values from each other. 'none' and 'bits' both\n"
		     << "\t  require --count n, for the same reason: their values cannot be\n"
		     << "\t  told apart once a value needs more than one character.\n\n";
		return -1;
	}

    /*if (argc != 5)
    {
        cout << "usage: " << argv[0] << " <input path> <accessions path> <file index> <min occurence threshold (across panel)>\n";
        return -1;
    }
	*/
    try
    {
      	/*string input_path = argv[1];
		string accessions_path = argv[2];
        uint file_index = stoi(argv[3]);
		uint min_occur = stoi(argv[4]);*/
		
      	if (input_path[input_path.size()-1] != '/')
		input_path += '/';
    	auto accessions = get_accessions(accessions_path);
   		size_t NUM_ACC = accessions.size();

   		std::string ouput_dir = "";
   		std::string delim = bit_packed ? "bits" : (delimiter == "\t" ? "tab" : (delimiter == " " ? "space" : "none"));
		if (show_count)
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_k"+to_string(kmer::K)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_count_delim-"+delim+"/";
		}
		else
		{
			ouput_dir = "matrix_acc"+to_string(NUM_ACC)+"_k"+to_string(kmer::K)+"_bins"+to_string(num_bins)+"_minOcc"+to_string(min_occur)+"_pres-abs_delim-"+delim+"/";
		}
   		
		if (!filesystem::exists(ouput_dir)){
			filesystem::create_directory(ouput_dir);
		}

        cout << "***************************** " << endl;
        cout << "PROCESSING MATRIX CHUNK: " << file_index << endl;
        cout << "Number of accessions: " << NUM_ACC << std::endl;
        cout << "Output Folder: " << ouput_dir << std::endl;
        
        auto start = chrono::steady_clock::now();

        merge_chunk(file_index, min_occur, input_path, accessions_path, delimiter, show_count, ouput_dir, write_core, bit_packed, num_threads);

        auto end = chrono::steady_clock::now();
        cout << "processing index: " << file_index << " took "
             << chrono::duration_cast<chrono::seconds>(end - start).count()
             << " sec" << endl;

        cout << "FINISHED MATRIX CHUNK: " << file_index << endl;
        cout << "***************************** " << endl;
    }

    catch (exception const &e)
    {
        cerr << e.what() << endl;
        return EXIT_FAILURE; // This ensures a non-zero exit code
    }
    return EXIT_SUCCESS; // This ensures a zero exit code for success
}

