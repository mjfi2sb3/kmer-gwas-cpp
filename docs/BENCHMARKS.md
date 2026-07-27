# Benchmarks and validation

Measurements for the bounded-memory redesign of the counting and merging
stages, against the previous implementation (`3297c96`). Every figure here was measured, not estimated;
where something is a projection it says so.

---

## Real-data validation

Four rice accessions from the 3000 Rice Genomes Project (`ERS467753`–`ERS467756`),
12 GB gzipped paired FASTQ, 83 bp reads, k=51, 8 bins. Both implementations were
built from source and run over identical input on one exclusive Ibex node
(40 cores, 376 GB).

Reproduce with:

```bash
tools/compare_to_baseline.sh --data <fastq_dir> --accessions <list> --bins 8 --k 51
```

### Correctness

```
rows: baseline 396,551,209   current 396,551,210

rows only in baseline            : 0
rows only in current             : 1
  of which poly-A (expected)     : 1
  unexplained                    : 0
```

**396.5 million matrix rows, one intentional difference, zero unexplained.**

That single extra row is the poly-A k-mer (`AAA…A`). The baseline silently
discarded it: `A` encodes as `00`, so a genuine poly-A k-mer is all-zero bits,
which the old code also used as its "window contained a non-ACGT base" marker.
Poly-A is its own canonical form (its reverse complement `TTT…T` sorts higher),
so it could not escape by canonicalising to something else — it was simply lost
wherever a poly-A tract occurred.

The comparison also normalises away one other intentional change: the baseline
emitted a spurious second tab after the k-mer column (`KMER\t\t0\t1`).

Row *order* differs by design — the baseline iterated a hash map, the current
version emits sorted k-mers — so rows are compared as sets.

### Peak memory

| stage | baseline | current |
|---|---|---|
| `kmer_count` | **57.1 GB** | **27.6 GB** |
| `matrix_merge` | **6.4 GB** | **< 0.1 GB** |

Sampled at 2 s intervals across the whole run. The `matrix_merge` figure for the
current version stayed below the sampler's 0.1 GB resolution.

The two numbers differ in kind, not just magnitude:

- The **baseline** scales with the input. It rose from 42.9 GB on the 1.3 GB
  accession to 57.1 GB on the 1.8 GB accession, because it held roughly three
  copies of every read in memory before counting anything.
- The **current** version is bounded by its configured budget (32 GB here) and
  spills to a single temporary file if the budget is exceeded, so peak memory
  does not track genome size or sequencing depth at all.

The practical effect is that the memory a job needs can be set in advance
rather than discovered from the data.

### Intermediate I/O

For one accession (`ERS467753`, 219,272,515 k-mer records):

| | baseline | current |
|---|---|---|
| final output | 3.95 GB (8 × `_nr.bin`) | 3.95 GB (1 × `.kbin`) |
| transient intermediates | **~12 GB written, then read back** | none |

Final output is byte-for-byte the same size. The difference is that the baseline
wrote every hash-map entry to per-bin shard files and re-read the lot in a
second pass to deduplicate — because its counting chunks were small, a k-mer
seen *n* times was written to disk close to *n* times.

---

## Inode use

Measured by sampling the working directory during a run, at 1500 bins:

| | baseline | current |
|---|---|---|
| Stage 1 peak inodes per job | **4,502** | **2** |
| Stage 2 inodes per bin job | **25,201** | **0** |

Stage 1 previously created `1 + num_bins` directories and `2 × num_bins` files
per accession; it now writes one pack file (plus one spill file only if the
memory budget is exceeded).

Stage 2 previously extracted `<acc>/<bin>_nr.bin` from every accession's tar
into a local directory — 2 inodes per accession per bin job, live concurrently
across every running job. It now seeks directly to the bin it needs inside each
pack, so there is no extraction step.

At 12,600 accessions and `queueSize = 200`, that Stage 2 fan-out was
approximately 5 million concurrent inodes.

> **Note on tar:** the extraction step was *not* a bandwidth problem. Measured
> GNU tar 1.34 seeks over member payloads rather than scanning — extracting one
> member from a 1001 MB / 500-member archive read 6 MB with 499 `lseek()` calls,
> flat regardless of member position. The tar packaging cost only ~1.3× in read
> volume. The case for replacing it is inodes and metadata operations.

---

## Component measurements

These were taken on synthetic data to isolate individual changes.

### k-mer codec

Rolling canonical encoder versus the previous `substr` + `canonical` +
`bit_encode` path, which built a 2k-character `std::string` per k-mer:

| | throughput |
|---|---|
| baseline (string building) | 1.14 M k-mer/s |
| current (rolling, bit ops) | **531 M k-mer/s** |

The codec was previously the bottleneck — slower than the ~16 M k-mer/s
downstream accumulation. It no longer is.

### Accumulation strategy

Per-bin partition + sort + run-length compaction versus `unordered_map`
insertion, single-threaded, random keys:

| | throughput |
|---|---|
| `unordered_map` insert | 9.3–10.1 M k-mer/s |
| flat records + global sort | 10.7–10.9 M k-mer/s |
| **per-bin partition + sort** | **16.3–16.9 M k-mer/s** |

Sorting is ~1.6× *faster* than hashing here: random inserts into a large map are
cache- and TLB-miss bound, while per-bin sorts stay cache-resident. Compaction
itself is ~1% of the time. Sorting also produces exactly the ordering the pack
format and the Stage 2 merge require.

### Stage 2 memory scaling

The property that matters is *which variable* memory tracks. Same panel size
(200 accessions), different union sizes:

| union k-mers per bin | baseline (hash map) | current (k-way merge) |
|---|---|---|
| ~7,600 | 6 MB | 18 MB |
| ~500,000 | **237 MB** | **18 MB** |

The baseline grows with the union; the current version does not. That is
`O(union × accessions)` → `O(accessions)`, and it is what decouples `num_bins`
from genome size — previously bin count had to rise on large genomes to bound
Stage 2 memory, and bin count was what drove inode use.

### Batch size

Larger read batches are both faster and lighter, because a worker's map collapses
duplicate k-mers before they are spilled. Measured on a 200 kb genome at high
coverage (1.5 M read pairs):

| batch | time | peak RSS | | batch | time | peak RSS |
|---|---|---|---|---|---|---|
| 1024 | 21.5 s | 5.14 GB | | 16384 | 4.2 s | 1.69 GB |
| 2048 | 18.2 s | 5.02 GB | | 32768 | 2.5 s | 1.15 GB |
| 4096 | 14.8 s | 3.36 GB | | 65536 | 1.7 s | 1.00 GB |

The default is 16384 rather than the fastest value: that benchmark uses a small
genome, where a batch's map saturates quickly. On a real genome nearly every
k-mer in a batch is distinct, so map size scales with batch size, and 65536
would mean roughly 390 MB per in-flight task.

### k configurability

`k` is a compile-time constant (`-DKMER_K`), because the binaries are already
built at job start to get `-march=native`. Making it a runtime value would cost:

| | compile-time | runtime |
|---|---|---|
| k=51, two words | 550 M k-mer/s | 493 M k-mer/s (+11%) |
| k=31, one word | 680 M k-mer/s | 544 M k-mer/s (+25%) |

Either would be acceptable — the codec runs ~30× faster than the downstream
accumulation — but compile-time is free here.

For `k ≤ 32` the key fits one 64-bit word, so records are 10 bytes instead of
18: **44% less Stage 1 output**, and ~1.23× faster.

### Output format

Row width grows linearly with the panel, so at large cohorts the matrix
dominates total output:

| accessions | `tab` | `none` | `bits` |
|---|---|---|---|
| 1,000 | 2,051 B | 1,052 B | 302 B |
| 12,600 | 25,251 B | 12,652 B | **3,202 B** |

`--delimiter bits` packs presence/absence one bit per accession as hex, ~8×
smaller than `tab` at 12,600 accessions, with no change to the data it encodes.

---

## Test coverage

- **Codec equivalence**: 400,000 random k-mers including 23,530 containing
  non-ACGT bases — encode bits, on-disk bytes, decode, bin hash and ordering all
  identical to the baseline implementation.
- **Pack format**: unit tests at k=31 and k=51 covering write/validate/read-back
  with out-of-order random access and empty bins, plus rejection of missing
  files, bad magic, out-of-range bins, out-of-order writes and incomplete packs.
- **Cross-k safety**: a pack written at one k is rejected with a specific error
  when opened by a binary built for another, in both directions.
- **Spill path**: forced with budgets down to 20 KB (16 spill blocks); packs are
  byte-identical across every budget from 0.00002 GB to 32 GB.
- **k-mer correctness**: k=31 and k=51 outputs validated against an independent
  Python implementation — correct length, all canonical, counts exact.
- **Format equivalence**: FASTQ, FASTA and gzip inputs produce identical k-mers
  from the same underlying sequences.
- **Pipeline**: 8 configurations end-to-end through Nextflow
  (k ∈ {31, 51} × `--count` y/n × `--core` y/n), all matching the baseline.
