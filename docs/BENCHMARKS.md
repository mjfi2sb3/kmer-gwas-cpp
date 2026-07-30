# Benchmarks

This document reports what the pipeline costs to run and shows that its output is
deterministic and reproducible. It is meant to stand on its own. The first
section explains what the pipeline computes; the rest reports measured time,
memory, disk, and filesystem use for each of the two stages.

## What the pipeline computes

Each accession is one sequenced individual, supplied as a pair of gzipped FASTQ
files. Its reads are cut into overlapping k-mers (substrings of a fixed length
k), and each distinct k-mer is counted.

Each distinct k-mer is assigned to one of `num_bins` bins by a hash of the k-mer.
The assignment is deterministic: a given k-mer always falls in the same bin, in
every accession. This is what lets the matrix be built in parallel. Because each
k-mer belongs to exactly one bin, the bins are independent, so the second stage
builds the matrix one bin at a time as separate jobs running concurrently, rather
than as a single serial pass over the entire k-mer set. That serial pass grows
with both the number of distinct k-mers and the number of accessions, so removing
it is what keeps matrix construction fast as the panel grows. The work then runs
in two stages:

- Stage 1 (`kmer_count`) runs one job per accession. It counts that accession's
  k-mers and writes them into a single pack file, `<accession>.kbin`, which holds
  every bin back to back with an offset table so any one bin can be read by
  seeking straight to it.
- Stage 2 (`matrix_merge`) runs one job per bin. It reads that bin's slice from
  every accession's pack and produces the matrix: one row per distinct k-mer, one
  column per accession.

Two design choices drive the measurements below. Stage 1 caps its memory with a
configurable budget instead of letting memory grow with the input. Stage 2 merges
the per-accession slices as sorted streams, so its memory grows with the number
of accessions and not with the number of distinct k-mers.

## The input data

The measurements use real sequencing data from the 3000 Rice Genomes Project.
Rice (Oryza sativa) has a genome of roughly 380 to 400 megabases, a mid-sized
plant genome. Each accession used here (`ERS467753` and its siblings) is one
sequenced rice individual, supplied as a pair of gzipped FASTQ files of about
1.3 GB each, so about 2.6 GB per accession on disk. That expands to roughly
7.4 GB of uncompressed sequence, about 33 million paired short reads, and at
k = 51 it produces about 219 million distinct k-mers per accession.

These figures set the scale for everything below: the resource numbers are the
cost of turning one such accession into its k-mer pack, and of merging a panel of
them into the matrix. The number of distinct k-mers, and so the pack size and the
merge work, grows with genome size and genetic diversity, so a larger genome such
as wheat would produce substantially more. The budget and streaming mechanisms
described below are what keep that growth bounded in memory.

## How these numbers were measured

Results come from two kinds of run, both on the rice data above at k = 51,
1500 bins, and default filters unless noted otherwise:

- Pipeline runs: a full `nextflow run` over a panel of 10 accessions, with
  per-task time, memory, and I/O taken from Nextflow's own execution report.
- Direct runs: the `kmer_count` and `matrix_merge` binaries run on their own to
  isolate one variable at a time, such as the memory budget or the accession
  count. These ran on a single node with 40 cores and 394 GB of RAM.

## Determinism and reproducibility

The pipeline's output does not depend on how much memory Stage 1 is allowed to
use. The same accession was counted with budgets from 2 GB, which forces repeated
spilling to disk, up to 64 GB, which holds everything in memory, and every run
produced a byte-identical pack file (the same md5 checksum). The budget changes
only how often k-mers are compacted or spilled, never the counts. Stage 2 reads
those packs in a fixed bin and accession order, so the matrix it writes is
deterministic in the same way.

## Stage 1: kmer_count

Measured per accession over the 10-accession pipeline run (median values), for a
rice accession of about 219 million distinct k-mers:

| measure | value |
|---|---|
| wall time | 23 s |
| pack written | one file, about 3.95 GB |
| total data written | 4.2 GB (the pack; nothing is written and then re-read) |
| total data read | 3.1 GB (the compressed input) |
| peak memory | set by the budget (see below) |
| filesystem entries created | 1 pack file, plus 1 spill file only if the budget is exceeded |

Memory is bounded by the budget, not by the input. Stage 1 streams reads and
discards them as it goes, and it caps k-mer accumulation at a budget that
defaults to 70% of the memory actually enforced on the job (read at run time from
the cgroup limit or the SLURM allocation). When a bin exceeds its share of the
budget it is sorted and compacted in place, and only if it is still over does a
compacted block spill to a single temporary file. Peak resident memory therefore
tracks the budget rather than the genome size. On one rice accession a 64 GB
budget peaked at 26 GB with no spilling, while a 2 GB budget on the same
accession peaked at 16 GB and spilled 9,392 blocks; both produced the identical
pack.

Within a Stage 1 job the two costs are decompressing and counting the reads, and
then sorting each bin before writing the pack. Both run in parallel: the two mate
files are decompressed on separate threads, and the per-bin sorts are spread
across the thread pool rather than done one bin at a time. This is why the job
keeps many cores busy instead of stalling on a single-threaded step.

## Stage 2: matrix_merge

Stage 2's memory grows with the number of accessions, not with the genetic
diversity of the panel. It was measured directly by opening N accession packs and
merging a single bin:

| accessions | peak memory |
|---|---|
| 100 | 11 MB |
| 1,000 | 79 MB |
| 5,000 | 382 MB |
| 12,600 (extrapolated) | about 1.0 GB |
| 50,000 (extrapolated) | about 3.8 GB |

That is about 76 KB per accession, and the relationship is linear. The merge for
one bin holds a small read buffer per accession and passes k-mers through a heap
one at a time, so nothing it holds scales with the number of distinct k-mers. At
panel sizes in the tens of thousands the binding limit becomes the number of open
file descriptors, one pack held open per accession, rather than memory, so a
memory request of a few GB is ample. Each bin's merge takes a few seconds and
writes one compressed matrix file, with no separate extraction step.

## Performance across releases

The pipeline has been measured on the same 10-accession, 1500-bin panel at each
released version, so the effect of each release is directly comparable. Wallclock,
CPU, memory, and I/O come from the Nextflow execution reports; the inode counts are
measured separately from the pack format, since the reports do not record them.

The baseline is v2.5.1, the original design the project started from. (The v2.5
series began at v2.5.0; these numbers are from its final release, v2.5.1.) It counts
an accession by loading its entire read set into memory, writing every k-mer
occurrence to per-bin files on disk, and then re-reading those files to remove
duplicates. Its memory grows with the input, and each accession writes and re-reads
many gigabytes of intermediate files. Each later release changes one part of that
design:

- **v3.1.1** is the redesign. Reads are streamed and discarded instead of held in
  memory, k-mer accumulation is capped by a budget, each accession is written as a
  single pack file with nothing re-read, and Stage 2 changes from a hash table over
  each bin to a streaming merge whose memory depends on the number of accessions
  rather than the number of distinct k-mers. This is where CPU, I/O, and inodes drop
  by an order of magnitude, and where RAM shifts from holding the reads to holding
  the k-mer table.
- **v3.2.0** parallelises the per-bin sort that finalises each pack, previously done
  one bin at a time.
- **v3.3.0** sizes the memory budget automatically from the memory actually enforced
  on the job. This changes correctness on exclusive nodes and workstations, not speed.
- **v3.4.0** decompresses the two mate files on separate threads with a faster gzip
  library, so the counting workers no longer wait on single-threaded decompression.

The tables below put numbers on those changes. Stage 1 (`kmer_count`) carries almost
all of the change between releases. Every
resource category, per accession, is below. Wallclock, cores, peak RAM, and I/O
are medians over the 10 jobs; CPU-time is the total across all 10; inodes is the
number of filesystem entries one job creates.

| version | wallclock | cores | CPU-time (10 jobs) | peak RAM | read | written | inodes/job |
|---|---|---|---|---|---|---|---|
| v2.5.1 (baseline) | 68 s | ~100 | 75,775 core-s | 20.6 GiB | 30.0 GB | 31.1 GB | ~4,500 |
| v3.1.1 | 126 s | ~5 | 7,719 core-s | 39.4 GiB | 3.1 GB | 4.2 GB | 2 |
| v3.2.0 | 43 s | ~16 | 7,839 core-s | 39.2 GiB | 3.1 GB | 4.2 GB | 2 |
| v3.3.0 | 44 s | ~16 | 7,949 core-s | 39.4 GiB | 3.1 GB | 4.2 GB | 2 |
| v3.4.0 | 23 s | ~38 | 9,647 core-s | 39.4 GiB | 3.1 GB | 4.2 GB | 2 |

Two columns are easy to misread:

- **Wallclock and cores.** v2.5.1's 68 seconds looks fast but is misleading: it
  ran on about 100 cores while requesting 32, so it used cores it was not
  allocated, and on a node that enforces the request it would be far slower. From
  v3.1.1 the core count is honest, and wallclock then tracks how much of the work
  is parallel: 126 seconds at about 5 cores while the heavy steps were serial, 43
  at 16 cores once sorting was parallelised, and 23 at 38 cores once decompression
  was too. The honest measure of total work is CPU-time, where every v3 release is
  8 to 10 times below the baseline.
- **RAM.** This is the one category where the current design uses more, about 39
  GiB against 20.6, and it is deliberate. The baseline held the reads in memory and
  streamed the growing k-mer table to disk; the current design discards the reads
  and holds the compacted table in memory instead, the trade that removes the 30 GB
  read-back in the I/O columns. The peak is bounded by the configured budget, not by
  the input, and a smaller budget lowers it: the same accession peaks at 16 GiB
  under a 2 GB budget, below the baseline.

The remaining columns move the same way from baseline to v3.4.0: reads fall from
30.0 to 3.1 GB and writes from 31.1 to 4.2 GB, since the current design compacts in
memory and writes one pack rather than writing and re-reading per-bin files; and
inodes fall from about 4,500 files and directories per job to 2 (one pack, plus one
spill only if the budget is exceeded). On shared HPC filesystems the inode quota
often binds before disk space, so that last column is frequently the one that matters
most.

Stage 2 (`matrix_merge`) changes less between releases. Per bin, over 1500 jobs:

| version | wallclock | cores | CPU-time (1500 jobs) | peak RAM | read | written | temp inodes/job |
|---|---|---|---|---|---|---|---|
| v2.5.1 (baseline) | 16 s | ~0.5 | 12,228 core-s | 0.1 GiB | 0.3 GB | 0.1 GB | ~25,000 |
| v3.1.1 | 2 s | ~1.8 | 6,458 core-s | 0.1 GiB | 0.1 GB | <0.1 GB | 0 |
| v3.2.0 | 3 s | ~1.7 | 6,472 core-s | 0.1 GiB | 0.1 GB | <0.1 GB | 0 |
| v3.3.0 | 3 s | ~1.6 | 6,481 core-s | 0.1 GiB | 0.1 GB | <0.1 GB | 0 |
| v3.4.0 | 3 s | ~1.7 | 6,471 core-s | 0.1 GiB | 0.1 GB | <0.1 GB | 0 |

At 10 accessions Stage 2 memory is 0.1 GiB in every version; how it scales with
accession count is covered in the Stage 2 section above. The releases differ in
wallclock (16 to 3 seconds), in CPU-time (about half), and most of all in
temporary inodes: v2.5.1 extracted one file per accession from every pack before
merging, about 25,000 temporary entries per bin job, while the current design
seeks straight into the packs and creates none.

## Summary

- Across releases, per-accession Stage 1 cost fell to about an eighth of the CPU
  work and a third of the wall time, at identical output.
- Stage 1 peak memory is a configuration value and does not grow with input size.
- Stage 1 writes one pack per accession and creates no intermediate files that are
  later re-read.
- Stage 2 memory is about 76 KB per accession and is independent of genetic
  diversity, staying under a few GB even at tens of thousands of accessions.
- The output is deterministic and byte-identical regardless of the Stage 1 memory
  budget.
