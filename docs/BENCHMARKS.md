# Benchmarks and validation

## What is being compared

Two versions of **this same pipeline**, running the same analysis on the same
input:

| | |
|---|---|
| **"before"** | the pipeline as of commit `3297c96`, the last revision prior to the rewrite of the counting and merging stages |
| **"after"** | the current tree |

Nothing else differs. Same reads, same *k*, same number of bins, same flags,
same machine. The tables below use these two labels throughout.

The rewrite changed *how* the two stages compute their result, not *what* the
result is. So the first question these measurements answer is "does it still
produce the same matrix?", and only then "is it cheaper?".

## What the pipeline computes, and the vocabulary used below

Each **accession** is one sequenced individual — a pair of FASTQ files. Its
reads are chopped into overlapping **k-mers**, substrings of fixed length *k*,
and each distinct k-mer is counted.

Because the full k-mer table across a whole panel is far too large to build in
one piece, k-mers are hashed into `num_bins` **bins**. Each bin is processed as
an independent job, and the pipeline runs in two stages:

- **Stage 1 (`kmer_count`)** — one job per accession. Counts that accession's
  k-mers and writes them out grouped by bin.
- **Stage 2 (`matrix_merge`)** — one job per bin. Reads that bin's k-mers from
  *every* accession and produces the matrix: one row per distinct k-mer, one
  column per accession.

Three further terms appear below:

- **union** — the number of *distinct* k-mers across all accessions in a bin,
  i.e. the number of rows that bin's matrix will have. It grows with the genetic
  diversity of the panel, not with any single genome.
- **pack** (`.kbin`) — the "after" version's Stage 1 output: one file per
  accession containing every bin back-to-back, plus an offset table so any bin
  can be read by seeking. The "before" version instead wrote a directory of
  per-bin files and tarred them.
- **inodes** — filesystem entries: every file *and* every directory counts as
  one. Shared HPC filesystems usually impose a quota on these separately from
  disk space, and it is frequently the quota that binds first.

---

# Part 1 — Real-data validation

This is the measurement that matters most: the two versions run on real
sequencing data, and their outputs compared.

**Input:** four rice accessions from the 3000 Rice Genomes Project
(`ERS467753`–`ERS467756`) — 12 GB of gzipped paired FASTQ, 83 bp reads.
**Settings:** k=51, 8 bins, default filters.
**Machine:** one exclusive Ibex node, 40 cores, 376 GB RAM.
**Method:** both versions compiled from source, run over identical input,
outputs diffed.

Reproduce on your own data with:

```bash
tools/compare_to_baseline.sh --data <fastq_dir> --accessions <list> --bins 8 --k 51
```

## Does it produce the same matrix?

```
rows: before 396,551,209   after 396,551,210

rows only in before              : 0
rows only in after               : 1
  of which poly-A (expected)     : 1
  unexplained                    : 0
```

**396.5 million matrix rows compared. Zero unexplained differences.**

Two differences are expected, and the comparison accounts for both:

**1. One extra row — the poly-A k-mer (`AAA…A`).** The "before" version silently
discarded it. `A` encodes as bit pattern `00`, so a genuine poly-A k-mer is all
zero bits — and the old code used that same all-zero value as its marker for
"this window contained an ambiguous base". The two were indistinguishable, so
the real k-mer was thrown away with the markers. Poly-A is also its own
canonical form (its reverse complement `TTT…T` sorts higher), so it could not
avoid the collision. Wherever a poly-A tract occurs in the data, that k-mer was
lost. This is now fixed, hence one additional row.

**2. A spurious tab.** The "before" version emitted a second tab after the k-mer
column (`KMER\t\t0\t1` rather than `KMER\t0\t1`). The comparison normalises this
away before diffing, so it does not mask any other difference.

Row *order* also differs, by design: the old version emitted rows in hash-table
order, the new one emits them sorted. Rows are therefore compared as sets.

## Memory

Peak resident memory, sampled every 2 seconds across the whole run:

| stage | before | after |
|---|---|---|
| `kmer_count` (Stage 1) | **57.1 GB** | **27.6 GB** |
| `matrix_merge` (Stage 2) | **6.4 GB** | **< 0.1 GB** |

The headline ratio understates the change, because the two figures are
different *kinds* of number.

**"Before" was determined by the input.** Stage 1 loaded an accession's entire
read set into memory — in roughly three copies — before counting anything. Its
peak rose from 42.9 GB on the 1.3 GB accession to 57.1 GB on the 1.8 GB
accession, tracking input size directly. You could not choose it; you could only
discover it by running out.

**"After" is determined by a setting.** Reads are streamed and discarded as they
are consumed, and k-mer accumulation is capped by a configurable budget (32 GB
in this run). If the budget is reached, data spills to a single temporary file
rather than the job growing. Peak memory is therefore whatever you set it to.

The Stage 2 figure for "after" stayed below the sampler's 0.1 GB resolution.
That stage's memory now depends on the number of accessions rather than the
number of distinct k-mers — with only 4 accessions here, it is negligible. See
[Stage 2 memory scaling](#stage-2-memory-scaling) for the measurement that
separates those two variables.

## Intermediate I/O

For a single accession (`ERS467753`, 219,272,515 k-mer records):

| | before | after |
|---|---|---|
| final output | 3.95 GB (8 × `_nr.bin`) | 3.95 GB (1 × `.kbin`) |
| transient intermediates | **~12 GB written, then read back** | **none** |

The final output is the same size — the same k-mers, the same 18 bytes each.
The difference is everything in between.

The "before" version counted reads in small chunks and wrote each chunk's
results straight to per-bin files, so a k-mer occurring *n* times across
different chunks was written to disk close to *n* separate times. A second pass
then read all of it back to merge the duplicates. The "after" version compacts
duplicates in memory before anything is written, so it writes the 3.95 GB once
and reads nothing back.

---

# Part 2 — Filesystem entry (inode) use

Measured by sampling the working directory during a run, at 1500 bins:

| | before | after |
|---|---|---|
| Stage 1 peak inodes per job | **4,502** | **2** |
| Stage 2 inodes per bin job | **25,201** | **0** |

**Stage 1** previously created, per accession, one directory plus `num_bins`
subdirectories plus two files in each — `1 + 1500 + 3000 = 4,501`, plus the
tar. It now writes one pack file, plus one spill file only if the memory budget
was exceeded.

**Stage 2** previously extracted one file per accession out of every tar into a
local directory before merging — two inodes per accession per bin job, live
simultaneously across every concurrently running job. At 12,600 accessions with
`queueSize = 200` that is roughly 5 million concurrent inodes. It now seeks
directly to the bin it needs inside each pack, so there is no extraction step
and no temporary directory at all.

This is also why the two problems were entangled before: the fix for Stage 2's
memory growth was to raise `num_bins`, but `num_bins` was exactly what drove
Stage 1's inode count.

> **A correction worth recording.** An earlier analysis in this work assumed the
> tar extraction was also a *bandwidth* problem, on the reasoning that tar must
> scan an archive linearly to find a member. That is wrong. Measured: GNU tar
> 1.34 seeks over member payloads — extracting one member from a 1001 MB,
> 500-member archive read 6 MB and issued 499 `lseek()` calls, flat regardless
> of the member's position. The tar layout cost only ~1.3× in read volume. The
> case for replacing it is inodes and metadata operations, not bandwidth.

---

# Part 3 — Component measurements

These were taken on synthetic data, to isolate one change at a time. Real-data
figures above are the ones that describe end-to-end behaviour; these explain
*why* those figures came out as they did.

## k-mer encoding

Converting a k-mer string into its packed binary form, done once per k-mer per
read position — the innermost loop of Stage 1.

| | throughput |
|---|---|
| before — build a 2*k*-character string, then parse it | 1.14 M k-mer/s |
| after — rolling encoder, bit operations only | **531 M k-mer/s** |

The old encoder built a `std::string` of `"00"`/`"01"`/`"10"`/`"11"` fragments
for every k-mer and converted that back to bits. At 1.14 M k-mer/s it was slower
than everything downstream of it, making it the limiting step. The new one
updates the encoding incrementally as it slides along the read.

## Accumulation strategy: sorting versus hashing

How Stage 1 collects counts. The obvious choice is a hash map; the measurement
says otherwise. Single-threaded, random keys:

| | throughput |
|---|---|
| `unordered_map` insertion | 9.3–10.1 M k-mer/s |
| flat records, one global sort | 10.7–10.9 M k-mer/s |
| **flat records, partitioned by bin then sorted** | **16.3–16.9 M k-mer/s** |

Sorting is ~1.6× *faster* than hashing here. Random insertions into a large hash
map miss cache and TLB on almost every access; per-bin sorts work on blocks
small enough to stay cache-resident. The run-length compaction that follows the
sort costs about 1% of the time.

Sorting also produces exactly the ordering that the pack format and the Stage 2
merge require, so it is not a trade-off — it is faster *and* it is what the rest
of the design needs.

## Stage 2 memory scaling

The single most important structural result. Stage 2's job is to combine every
accession's k-mers for one bin. The question is what its memory depends on.

Holding the panel size fixed at 200 accessions and varying only the union (the
number of distinct k-mers in the bin):

| union k-mers in the bin | before | after |
|---|---|---|
| ~7,600 | 6 MB | 18 MB |
| ~500,000 | **237 MB** | **18 MB** |

The "before" version built a hash table holding the whole bin at once, so its
memory grew with the union — i.e. with how genetically diverse the panel is.
The "after" version merges the accessions' pre-sorted k-mer lists through a
heap, holding one small buffer per accession, so its memory is **flat** in the
union.

At the small end the new version is *worse* in absolute terms (18 MB versus
6 MB), because those per-accession buffers are a fixed cost. That is the right
trade: the fixed cost is bounded and known, the old growth was neither.

This is what allows `num_bins` to be chosen for parallelism and output file size
rather than as a memory control.

## Read batch size

How many reads a worker thread processes at once. Measured on a 200 kb genome at
high coverage (1.5 M read pairs):

| batch | time | peak RSS | | batch | time | peak RSS |
|---|---|---|---|---|---|---|
| 1024 | 21.5 s | 5.14 GB | | 16384 | 4.2 s | 1.69 GB |
| 2048 | 18.2 s | 5.02 GB | | 32768 | 2.5 s | 1.15 GB |
| 4096 | 14.8 s | 3.36 GB | | 65536 | 1.7 s | 1.00 GB |

Larger batches are both faster *and* lighter, which is the opposite of the usual
expectation. The reason is deduplication: a bigger batch lets a worker's map
collapse more duplicate k-mers before any of them are written out, so less
reaches disk and far fewer short-lived allocations are made.

The default is **16384, not the fastest value**. That benchmark uses a small
genome, where a batch's k-mer map saturates quickly. On a real genome nearly
every k-mer in a batch is distinct, so map size scales with batch size, and
65536 would mean roughly 390 MB per in-flight task.

## Making *k* configurable

*k* was previously fixed at 51. It is now set with `-DKMER_K` at compile time,
which is free because the binaries are already compiled at job start (to get
`-march=native` for the actual node). The alternative — reading *k* at run time
— would cost:

| | compile-time | run-time |
|---|---|---|
| k=51, key spans two 64-bit words | 550 M k-mer/s | 493 M k-mer/s (+11%) |
| k=31, key fits one 64-bit word | 680 M k-mer/s | 544 M k-mer/s (+25%) |

Either would have been acceptable — the encoder runs roughly 30× faster than the
accumulation it feeds, so it is nowhere near the limiting step any more — but
compile-time costs nothing here.

The one-word case matters independently: for `k ≤ 32` a k-mer fits in a single
64-bit word, so each record is 10 bytes instead of 18. **44% less Stage 1
output**, and ~1.23× faster.

## Matrix row width

The matrix has one column per accession, so row width grows linearly with panel
size. At large cohorts the matrix dominates everything else the pipeline
produces.

| accessions | `--delimiter tab` | `--delimiter none` | `--delimiter bits` |
|---|---|---|---|
| 1,000 | 2,051 B | 1,052 B | 302 B |
| 12,600 | 25,251 B | 12,652 B | **3,202 B** |

`bits` stores presence/absence as one bit per accession, written as hex, rather
than one character plus one separator. It encodes exactly the same information —
verified by decoding it back and comparing to the `tab` output — and is ~8×
smaller at 12,600 accessions. It is presence/absence only, so it cannot be
combined with `--count y`.

---

# Part 4 — Test coverage

What was verified beyond the end-to-end comparison above.

| what | how |
|---|---|
| k-mer encoding is unchanged | 400,000 random k-mers including 23,530 containing ambiguous bases; encoded bits, on-disk bytes, decoded string, bin assignment and sort order all identical to the old implementation |
| pack file format | unit tests at k=31 and k=51: write, validate, read back with out-of-order random access, empty bins; plus rejection of missing files, bad magic numbers, out-of-range bins, out-of-order writes and incomplete files |
| a pack cannot be misread at the wrong *k* | packs written at k=31 and k=51 each rejected with a specific error by a binary built for the other, in both directions |
| the spill path works | forced with budgets down to 20 KB (16 spill events); output packs byte-identical across every budget from 0.00002 GB to 32 GB |
| k-mers are correct in absolute terms | k=31 and k=51 output checked against an independent Python implementation — correct length, all canonical, counts exact |
| input formats agree | FASTQ, FASTA and gzipped input produce identical k-mers from the same underlying sequences |
| the workflow drives the binaries correctly | 8 configurations end-to-end through Nextflow (k ∈ {31, 51} × `--count` y/n × `--core` y/n), each matching the old implementation |
