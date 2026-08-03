#ifndef KMER_KEY_HPP
#define KMER_KEY_HPP

// ===========================================================================
//  Packed k-mer key and per-accession count type.
// ===========================================================================
//
//  k IS A COMPILE-TIME CONSTANT, set with -DKMER_K=<odd number>.
//
//  The Nextflow modules compile these programs at job start (they already do,
//  to get -march=native for the actual node CPU), so k can be a pipeline
//  parameter without costing anything at run time. Measured alternative: making
//  k a runtime value costs 11% (2-word) to 25% (1-word) in the codec. That
//  would in fact be tolerable — the codec runs ~30x faster than the downstream
//  accumulation — but compile-time is free here, so we take it.
//
//  TWO REPRESENTATIONS, selected automatically:
//
//    k <= 32   ->  one uint64  (2k <= 64 bits).  Record = 10 bytes.
//    k <= 63   ->  two uint64  (2k <= 126 bits). Record = 18 bytes.
//
//  The single-word form is ~1.23x faster and 44% smaller on disk, so small k is
//  worth having rather than padding everything to 16 bytes.
//
//  LAYOUT (2-word form) is deliberately identical to libstdc++'s bitset<2k>,
//  so files written before the bitset -> packed-key change remain readable:
//    w0 holds bits 0..63, w1 holds bits 64..2k-1
//    the FIRST base of the k-mer occupies the two HIGHEST bits
//    each base is a 2-bit value A=0 C=1 G=2 T=3
//
//  Because the first base sits in the highest bits and A<C<G<T maps to 0<1<2<3,
//  numeric ordering of the key is identical to lexicographic ordering of the
//  decoded string. canonical() therefore reduces to min(forward, reverse).
//  Sorted key order == sorted k-mer order, which the Stage 2 merge relies on.
// ===========================================================================

#include <cstdint>
#include <cstddef>
#include <string>

#ifndef KMER_K
#define KMER_K 51
#endif

namespace kmer {

static constexpr int K    = KMER_K;
static constexpr int BITS = 2 * K;

// k must be ODD: with an even k a sequence can be its own reverse complement,
// which makes the canonical form ambiguous and would collapse distinct k-mers.
static_assert(K % 2 == 1, "KMER_K must be odd (even k makes canonical k-mers ambiguous)");
static_assert(K >= 15,    "KMER_K must be at least 15 to be informative");
static_assert(K <= 63,    "KMER_K must be <= 63 (the key holds 2k <= 126 bits in two words)");

static constexpr bool SINGLE_WORD = (BITS <= 64);

// ---------------------------------------------------------------------------
// Per-accession k-mer counts.
//
// 16 bits. Accumulation SATURATES at 65535 instead of wrapping — the previous
// behaviour silently wrapped, so a k-mer occurring 65,536 times was recorded as
// occurring 0 times and then dropped entirely by the singleton filter.
//
// Saturating rather than widening is deliberate: counts that high come
// from repeats, organellar sequence or contamination, where an exact figure
// carries no association signal. Widening would cost +11% on all Stage 1 output
// and would double the matrix row's per-accession width (two bytes to four).
// If exact high counts are ever needed, widen Count here and the record layout
// follows.
// ---------------------------------------------------------------------------
using Count = uint16_t;
static constexpr Count COUNT_MAX = 65535;

inline void sat_inc(Count& c)               { if (c < COUNT_MAX) ++c; }
inline void sat_add(Count& c, uint32_t add) {
    uint32_t s = (uint32_t)c + add;
    c = (s > COUNT_MAX) ? COUNT_MAX : (Count)s;
}

// ---------------------------------------------------------------------------
// Key
// ---------------------------------------------------------------------------
#if KMER_K * 2 <= 64
// ---- single-word form (k <= 32) ----
static constexpr uint64_t LO_MASK = (BITS >= 64) ? ~0ULL : ((1ULL << BITS) - 1);
static constexpr int      REV_TOP = 2 * (K - 1);

struct Key {
    uint64_t w0 = 0;

    bool operator==(const Key& o) const { return w0 == o.w0; }
    bool operator!=(const Key& o) const { return w0 != o.w0; }
    bool operator<(const Key& o) const  { return w0 < o.w0; }
    bool is_zero() const                { return w0 == 0; }

    inline void push_fwd(uint8_t v) { w0 = ((w0 << 2) | v) & LO_MASK; }
    inline void push_rev(uint8_t v) { w0 = (w0 >> 2) | ((uint64_t)(3 - v) << REV_TOP); }

    uint64_t bin_hash() const { return w0; }
    uint64_t low()      const { return w0; }
    uint64_t bit_at(int p) const { return (w0 >> p) & 1ULL; }
    uint64_t bits_at(int p) const { return (w0 >> p) & 3ULL; }
};

#else
// ---- two-word form (33 <= k <= 63) ----
static constexpr uint64_t HI_MASK = (BITS - 64 >= 64) ? ~0ULL : ((1ULL << (BITS - 64)) - 1);
static constexpr int      HI_TOP  = BITS - 64 - 2;

struct Key {
    uint64_t w0 = 0, w1 = 0;

    bool operator==(const Key& o) const { return w0 == o.w0 && w1 == o.w1; }
    bool operator!=(const Key& o) const { return !(*this == o); }
    bool operator<(const Key& o) const  { return w1 != o.w1 ? w1 < o.w1 : w0 < o.w0; }
    bool is_zero() const                { return w0 == 0 && w1 == 0; }

    inline void push_fwd(uint8_t v) {
        w1 = ((w1 << 2) | (w0 >> 62)) & HI_MASK;
        w0 = (w0 << 2) | v;
    }
    inline void push_rev(uint8_t v) {
        w0 = (w0 >> 2) | (w1 << 62);
        w1 = (w1 >> 2) | ((uint64_t)(3 - v) << HI_TOP);
    }

    // Legacy bin hash: the original took the top 64 bits of the 2k-bit string
    // via to_string().substr(0,64) then to_ulong(). Preserved exactly so that
    // k-mers keep landing in the same bins as before.
    uint64_t bin_hash() const { return (w1 << (128 - BITS)) | (w0 >> (BITS - 64)); }
    uint64_t low()      const { return w0; }
    uint64_t bit_at(int p) const { return (p >= 64) ? ((w1 >> (p - 64)) & 1ULL) : ((w0 >> p) & 1ULL); }
    uint64_t bits_at(int p) const { return (p >= 64) ? ((w1 >> (p - 64)) & 3ULL) : ((w0 >> p) & 3ULL); }
};
#endif

// On-disk / in-memory record: key followed by a 16-bit count, tightly packed.
// 10 bytes for k <= 32, 18 bytes otherwise.
#pragma pack(push, 1)
struct Record {
    Key   key;
    Count count;
};
#pragma pack(pop)

static constexpr size_t RECORD_BYTES = sizeof(Key) + sizeof(Count);

struct KeyHash {
    size_t operator()(const Key& k) const {
        uint64_t h = k.low() * 0x9E3779B97F4A7C15ULL;
#if KMER_K * 2 > 64
        h ^= (k.w1 + 0x165667B19E3779F9ULL);
#endif
        h ^= h >> 29; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 32;
        return (size_t)h;
    }
};

// 'A'->0 'C'->1 'G'->2 'T'->3, anything else -> 0xFF (invalid)
inline const uint8_t* base_table() {
    static uint8_t t[256];
    static bool init = [] {
        for (int i = 0; i < 256; i++) t[i] = 0xFF;
        t[(unsigned char)'A'] = 0; t[(unsigned char)'C'] = 1;
        t[(unsigned char)'G'] = 2; t[(unsigned char)'T'] = 3;
        return true;
    }();
    (void)init;
    return t;
}

// Encode a single k-mer, canonicalised. Returns false if it contains any
// non-ACGT base; callers map that to Key{} to reproduce the legacy behaviour of
// canonical()+bit_encode().
inline bool encode_canonical(const char* s, size_t len, Key& out) {
    if (len != (size_t)K) return false;
    const uint8_t* T = base_table();
    Key f, r;
    for (size_t i = 0; i < len; i++) {
        uint8_t v = T[(unsigned char)s[i]];
        if (v == 0xFF) return false;
        f.push_fwd(v);
        r.push_rev(v);
    }
    out = (r < f) ? r : f;
    return true;
}

inline bool encode_canonical(const std::string& s, Key& out) {
    return encode_canonical(s.data(), s.size(), out);
}

inline std::string decode(const Key& k) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    std::string out(K, 'A');
    for (int i = 0; i < K; i++)
        out[i] = B[k.bits_at(BITS - 2 - 2 * i) & 3ULL];
    return out;
}

inline uint64_t legacy_bin_hash(const Key& k) { return k.bin_hash(); }

// Low bits, matching matrix_merge's shard_of(). num_shards must be a power of 2.
inline size_t shard_of(const Key& k, size_t num_shards) {
    return (size_t)(k.low() & (num_shards - 1));
}

} // namespace kmer

#endif
