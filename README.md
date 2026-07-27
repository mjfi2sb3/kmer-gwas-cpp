# kmer-GWAS Pipeline

A high-performance Nextflow pipeline for genome-wide association studies (GWAS) using k-mer presence/absence or count matrices. Designed for large-scale deployment on HPC clusters with SLURM and Singularity/Apptainer.

---

## What's new

### What the pipeline does, in one paragraph

Every accession's sequencing reads are chopped into overlapping *k*-mers —
substrings of fixed length *k* — and counted. The k-mers are split into `num_bins`
groups by a hash, so the final table can be built one group at a time in parallel.
The output is a matrix with one row per distinct k-mer and one column per
accession, recording which accessions carry it. That matrix is the input to the
association test.

### The problem this release solves

Two things used to grow with the data instead of staying within limits you set.

**Memory.** Stage 1 read an accession's entire FASTQ into RAM — in roughly three
copies — before counting anything, so how much memory a job needed was a property
of the input rather than a number you could choose. Stage 2 then built a hash
table holding the whole bin's k-mer table at once, so *its* memory grew with how
many distinct k-mers the panel contained, i.e. with genetic diversity. The usual
workaround was to raise `num_bins` so each bin held less — which leads to the
second problem.

**File count (inodes).** Stage 1 created a directory tree per accession: one directory and
two files for every bin, so at 1500 bins that was ~4,500 filesystem entries
(inodes) per accession, all written and then read back to be deduplicated. Stage 2 then
extracted one file per accession out of each tar archive before merging, creating
two entries per accession per bin job, live simultaneously across every running
job. On a shared HPC filesystem with an inode quota this is often the binding
constraint, not disk space — and because the fix for the memory problem was more
bins, the two problems pulled against each other.

### What changed

**Stage 1 now streams.** Reads are parsed in batches, handed to worker threads,
and discarded. K-mers accumulate in memory only up to a configurable budget; when
that budget is reached they are sorted and compacted in place (which reclaims most
of the space, since a k-mer seen many times collapses to one entry), and only if
still over budget does anything go to disk. Peak memory is therefore the budget.

**Each accession's output is now a single file.** Instead of a tree of per-bin
files packed into a tar, Stage 1 writes one `.kbin` "pack" containing every bin
back to back, with an offset table and a footer. Reading bin *i* is a seek, not a
scan or an extraction.

**Stage 2 now merges instead of accumulating.** Because each pack's bins are
stored in sorted order, the matrix can be produced by advancing one cursor per
accession through a heap — take the smallest k-mer, record which accessions have
it, emit the row, move on. Only one small buffer per accession is resident, so
memory depends on how many accessions there are, not on how many distinct k-mers
they collectively contain. This is what breaks the tie between the two problems
above: `num_bins` is now purely a parallelism and output-file-size knob.

### Measured effect

Validated on real rice data against the previous implementation:
**396,551,209 vs 396,551,210 matrix rows, zero unexplained differences.**

| | before | after |
|---|---|---|
| Stage 1 peak memory | 57.1 GB, grows with input size | **27.6 GB, set by a budget** |
| Stage 2 peak memory | 6.4 GB, grows with k-mer diversity | **< 0.1 GB, grows only with accession count** |
| Stage 1 filesystem entries / inodes (1500 bins) | 4,502 per job | **2 per job** |
| Stage 2 filesystem entries / inodes | 25,201 per bin job | **0** |
| Intermediate I/O per accession | ~12 GB written, then read back | **none** |
| k-mer encoder throughput | 1.14 M k-mer/s | **531 M k-mer/s** |

### New options

- **`--kmer_size`** — *k* is no longer fixed at 51. Any odd value from 15 to 63,
  compiled in at job start so it costs nothing at run time. Values of 32 or below
  pack a k-mer into a single 64-bit word, making each record 10 bytes instead of
  18 and cutting Stage 1 output by ~44%.
- **`--encoding bits`** — writes presence/absence as one bit per accession rather
  than one character. At 12,600 accessions a matrix row shrinks from 25,251 to
  3,202 bytes, ~8× smaller, encoding exactly the same information.
- **`--min_kmer_count`** — drop k-mers seen fewer than *n* times within an
  accession. Low counts are overwhelmingly sequencing errors; this was previously
  hardcoded.

### Compatibility

Output is identical to the previous implementation apart from two deliberate
fixes: a spurious second tab after the k-mer column is gone, and poly-A k-mers
(`AAA…A`) are no longer silently discarded — the all-zero encoding of poly-A had
been colliding with the marker used for "this window contained an ambiguous base".
Anything that parses the matrix should be checked against the new format.

The k-mer length now appears in the output directory names, so results produced at
different *k* cannot be mixed up. Existing intermediate files from previous runs
are not readable by the new code and should be regenerated.

Full measurements, methodology and test coverage are in
**[docs/BENCHMARKS.md](docs/BENCHMARKS.md)**. The comparison is reproducible on
your own data with `tools/compare_to_baseline.sh`.

---

## Overview

The pipeline processes paired-end FASTQ files (plain or gzip-compressed) and produces a binary k-mer matrix across all accessions. It runs in two stages:

1. **KMER_COUNT**: one job per accession. Streams the paired FASTQ/FASTA files, extracts canonical k-mers, and writes a single self-indexed pack file (`<accession>.kbin`) holding that accession's k-mers for every bin.
2. **MATRIX_MERGE**: one job per bin. Seeks directly to that bin's slice inside each accession's pack and k-way merges the slices into a tabular matrix, applying the MAF filter and output format options, then compresses with pigz.

### Key features

- Streams gzip-compressed FASTQ directly via zlib; no temporary decompressed files on the shared filesystem
- Compiles C++ binaries at runtime with `-march=native`, targeting the actual node CPU for optimal performance
- Fully containerised: the build toolchain (GCC 12, zlib-dev, pigz) is embedded in the image; source code is compiled fresh on each job
- Works across HPC systems; no hardcoded paths or site-specific module dependencies
- Supports presence/absence or raw count output, configurable thresholds, and core k-mer extraction
- **Inode-efficient**: one pack file per accession and no intermediate files at all. Stage 1 peak inode use is 2 per job regardless of bin count (previously 1 + `num_bins` directories + 2 × `num_bins` files); Stage 2 needs no extraction step
- **Memory-bounded**: Stage 1 streams reads and caps k-mer accumulation with a budget, spilling to a single temporary file only if exceeded. Stage 2 merges sorted slices with memory proportional to the number of accessions, not to the number of distinct k-mers — so `num_bins` no longer has to grow with genome size
- Compressed output: matrix files are compressed with pigz (parallel gzip) if available, otherwise standard gzip

---

## Requirements

| Component | Minimum version |
|-----------|----------------|
| Nextflow  | 23.x (DSL2)    |
| Singularity / Apptainer | any |
| SLURM     | any            |

No compiler or zlib installation is required on compute nodes when using the `slurm_container` profile; everything is provided by the container image.

---

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/mjfi2sb3/kmer-gwas-cpp.git
cd kmer-gwas-cpp

# 2. Create accessions.txt — one accession ID per line, matching FASTQ file prefixes
cp accessions.txt.example accessions.txt
# edit accessions.txt

# 3. Place paired FASTQ files in a data directory
#    Supported naming: <accession>_1.fq, _1.fastq, _1.fq.gz, _1.fastq.gz  (and _2.*)
ls -d /path/to/fastqs/*.fastq.gz

# 4. Load needed modules
# Or Install nextflow and Singularity
module load nextflow singularity;

# 5. Print Help
nextflow run main.nf --help

# 6. Run on SLURM with Singularity (recommended)
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastqs/
```

Results are written to `./results/` by default.

---

## Input

### Accessions file

A plain-text file with one accession ID per line. Blank lines are ignored.

```
SAMPLE1_SRR0000001
SAMPLE2_SRR0000002
SAMPLE3_SRR0000003
```

### FASTQ files

Paired-end reads for each accession. Files must follow the naming convention:

```
<accession>_1.fq[.gz]   # R1
<accession>_2.fq[.gz]   # R2
```

Both plain (`.fq`, `.fastq`) and gzip-compressed (`.fq.gz`, `.fastq.gz`) files are supported. The pipeline auto-detects compression by magic bytes — file extension does not matter.

FASTA input (single- or multi-line records) is also accepted and is detected from the leading `>`. Input is **streamed**: the reads are never held in memory in full, so peak memory is independent of sequencing depth and genome size.

---

## Output

```
results/
├── run_manifest.txt              # parameters this results directory was produced with
├── bin/
│   └── bits_to_text              # compiled converter for the bit-packed matrices
├── kmer_count_k<k>/
│   └── <accession>.kbin          # one self-indexed pack file per accession
├── matrix/
│   └── matrix_acc<N>_k<k>_bins<B>_minOcc<t>_<mode>_delim-<d>/
│       ├── <bin>_matrix.tsv.gz   # k-mer matrix for this bin (gzip-compressed)
│       └── <bin>_core.txt.gz     # core k-mers (if --core y, gzip-compressed)
└── reports/
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

The k-mer length appears in both the Stage 1 directory name and the matrix
directory name, so results produced at different `k` cannot be mixed up or
silently reused.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--accessions_file` | `./accessions.txt` | Path to accessions list |
| `--data_dir` | `./data` | Directory containing paired FASTQ files |
| `--output_dir` | `./results` | Output directory |
| `--kmer_size` | `51` | k-mer length. Must be **odd** and between 15 and 63. `k ≤ 32` uses a compact 10-byte record instead of 18, cutting Stage 1 output by ~44%. See note below |
| `--num_bins` | `1500` | Number of k-mer bins. Now a parallelism / output-file-size knob: it no longer has to grow with genome size, because Stage 2 memory is independent of the number of distinct k-mers |
| `--min_kmer_count` | `2` | Drop k-mers seen fewer than this many times within an accession. Low counts are overwhelmingly sequencing errors; raising it shrinks Stage 1 output and everything downstream |
| `--threshold` | `0` | Two-sided MAF filter on the number of **accessions** carrying a k-mer. Keeps only k-mers present in `[threshold, num_accessions - threshold]` accessions; `0` disables it (see note below) |
| `--count` | `n` | `y` = raw counts, `n` = presence/absence |
| `--encoding` | `text` | Matrix encoding: `text` (delimited) or `bits` (1 bit per accession as hex, ~8× smaller; presence/absence only). See note below |
| `--delimiter` | `none` | Value separator for the **text** encoding: `tab`, `space` or `none`. `none` concatenates single characters and is presence/absence only |
| `--core` | `n` | `y` = write core k-mers file per bin. Core k-mers are **excluded from the matrix** (see note below) |
| `--matrix_merge_cpus` | `32` | Threads for the MATRIX_MERGE stage |
| `--kmer_count_memory` | `370.GB` | RAM per KMER_COUNT job. Use dot notation: `120.GB`, `370.GB`. Use `--clusterOptions='--mem=0'` to request all available node RAM instead |
| `--matrix_merge_memory` | `370.GB` | RAM per MATRIX_MERGE job. Use dot notation: `64.GB`, `120.GB`, `370.GB` |
| `--kmer_count_time` | `5h` | Wallclock time limit per KMER_COUNT job. Examples: `'2h'`, `'5h'`, `'1d'`, `'2h 30m'` |
| `--matrix_merge_time` | `10h` | Wallclock time limit per MATRIX_MERGE job. Examples: `'5h'`, `'10h'`, `'1d'`, `'2h 30m'` |
| `--cleanup` | `true` | Delete Nextflow work directory on successful completion. Pass `--cleanup false` to preserve work dirs for debugging or `-resume` |
| `--clusterOptions` | _(none)_ | Extra SLURM flags passed to every job (see note below) |
| `--singularity_cache_dir` | `./.singularity` | Local path for Singularity image cache |

> **`--clusterOptions` syntax:** because the value starts with `--`, you must use the `=` form to prevent Nextflow misinterpreting it as a flag:
> ```bash
> --clusterOptions='--account=myproject --partition=highmem'
> ```

> **`--kmer_size` is compiled in, not read at run time.** The C++ binaries are already built at the start of every job (to get `-march=native` for the actual node CPU), so `k` is passed as `-DKMER_K` and costs nothing at run time — the k-mer encoder runs at the same speed for any `k`.
>
> `k` must be **odd**: with an even `k` a sequence can be its own reverse complement, which makes the canonical form ambiguous. It must also be between 15 and 63, since the key holds `2k` bits in at most two 64-bit words. Invalid values are rejected at launch.
>
> For `k ≤ 32` the key fits in a single 64-bit word, so each record is 10 bytes rather than 18 — about 44% less Stage 1 output, and slightly faster. Choosing `k = 31` instead of `k = 51` is therefore a substantial storage saving if the shorter k-mer suits your analysis.
>
> **Changing `k` invalidates existing bin files**: they store fixed-width k-mers and cannot be read at a different `k`. Use a fresh `--output_dir`.

> **Matrix output format.** Three options combine to control the output. The
> k-mer is always tab-separated from the values; `--delimiter` separates the
> values from each other.
>
> - **`--count`** — what each cell holds: `n` presence/absence (default), or `y` a raw count.
> - **`--encoding`** — `text` (default) writes the values as characters; `bits` packs presence/absence one bit per accession, written as hex.
> - **`--delimiter`** (text encoding only) — `none` (default) concatenates, `tab`, or `space`.
>
> | encoding | count | delimiter | a row (present in acc 0 & 2) |
> |---|---|---|---|
> | text | n | none *(default)* | `KMER⇥101` |
> | text | n | tab | `KMER⇥1⇥0⇥1` |
> | text | n | space | `KMER⇥1 0 1` |
> | text | y | tab / space | `KMER⇥5⇥0⇥3` |
> | text | y | none | **rejected** — multi-digit counts would merge |
> | bits | n | *(ignored)* | `KMER⇥05` (hex) |
> | bits | y | — | **rejected** — bits is presence/absence only |
>
> `none` and `bits` both require `--count n` for the same reason: their values
> cannot be told apart once a value needs more than one character.
>
> **Size at scale**, per matrix row:
>
> | accessions | `tab` | `none` | `bits` |
> |---|---|---|---|
> | 1,000 | 2,051 B | 1,052 B | 302 B |
> | 12,600 | 25,251 B | 12,652 B | **3,202 B** |
>
> The matrix usually dwarfs every intermediate file in this pipeline, so this is
> the single largest lever on total output size.
>
> Use **`bits_to_text`** to convert between the two forms — `--decode` (bits →
> text, the default) and `--encode` (text → bits, to compress an existing text
> matrix). It needs the accession count, to strip or check padding bits. Each
> run publishes a compiled copy to `results/bin/bits_to_text`, so the output
> directory ships with the tool that reads it:
>
> ```bash
> # decode a bits matrix to tab-separated text
> results/bin/bits_to_text -a accessions.txt --delimiter tab \
>     results/matrix/matrix_*/0_matrix.tsv.gz  0_text.tsv.gz
>
> # or compress a text matrix to bits
> results/bin/bits_to_text --encode -a accessions.txt --delimiter tab \
>     0_text.tsv.gz  0_bits.tsv.gz
> ```
>
> A portable Python version, `tools/bits_to_text.py`, has the same interface and
> byte-identical output (the C++ one is ~14× faster for large matrices). Or
> decode a row by hand: split on the tab, `bytes.fromhex(rest)`, then accession
> *i* is present if `(bs[i >> 3] >> (i & 7)) & 1`.

> **`--threshold` is a minor-allele-frequency filter, not a count filter.** It is applied to the number of accessions in which a k-mer occurs — *not* to the k-mer's count within an accession. It is two-sided: `--threshold 20` keeps k-mers found in at least 20 accessions **and** in at most `num_accessions - 20`, discarding both rare and near-ubiquitous k-mers, since neither carries usable association signal. The default of `0` disables the filter entirely, on the assumption that downstream analysis applies its own MAF cutoff.

> **Core k-mers are excluded from the matrix by design.** With `--core y`, k-mers present in *every* accession are written to `<bin>_core.txt.gz` and omitted from `<bin>_matrix.tsv.gz`. They are invariant across the panel, so they carry no association signal — this file is a record of them, not an additional part of the matrix. Note the consequence: with any `--threshold` ≥ 1, the core file and the matrix are wholly disjoint.

---

## Execution profiles

### `slurm_container` (recommended)

Runs on SLURM using a pre-built Singularity/Apptainer image from GHCR. No compiler or library modules need to be loaded on compute nodes.

```bash
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

On HPC systems that require a project account or specific partition:

```bash
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq \
    --clusterOptions='--account=myproject --partition=highmem'
```

| Stage | CPUs | Memory | Time |
|-------|------|--------|------|
| KMER_COUNT | 32 (fixed) | 370.GB (`--kmer_count_memory`) | 5h (`--kmer_count_time`) |
| MATRIX_MERGE | 32 (`--matrix_merge_cpus`) | 370.GB (`--matrix_merge_memory`) | 10h (`--matrix_merge_time`) |

The container image is pulled automatically on first run and cached in `.singularity/` under the launch directory. Override the cache location with `--singularity_cache_dir /path/to/cache` (useful for sharing the cache across multiple runs).

### `slurm`

Runs on SLURM using the host environment. Requires GCC 12+ and zlib-dev to be available (e.g., via the module system).

```bash
nextflow run main.nf \
    -profile slurm \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

### `standard`

Runs locally using all available CPUs (up to 64). Useful for small-scale testing.

```bash
nextflow run main.nf \
    -profile standard \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastqs/ \
    --num_bins 5
```

---

## Container

The pipeline image is available on the GitHub Container Registry:

```
ghcr.io/mjfi2sb3/kmer-gwas-cpp:v3.1.0
```

**The image is a toolchain, not a copy of the code.** It provides GCC 12,
zlib-dev and pigz; the C++ sources compiled inside it are the ones in **your
checkout** (`src/`), which is bind-mounted into the container. Binaries are
built at job start with `-march=native`, so they are optimised for the actual
compute node CPU. The image contains no pre-compiled binaries.

This matters for a reason worth stating plainly: the container profile used to
compile from a copy of the source baked into the image at `/opt/kmer-gwas/src/`.
That meant the code that actually ran was whatever had been frozen into the
pinned image tag, so every source change needed a new git tag, a CI image build
and a pin bump before it took effect — and if any of those lagged, the
`slurm_container` profile silently ran *different code* from the `standard` and
`slurm` profiles. Compiling from the checkout removes that possibility. The
image now only needs rebuilding when the compiler or libraries change, not when
the pipeline changes.

pigz is used for parallel gzip compression of matrix output files; if pigz is
unavailable (e.g. on the `slurm` profile with a host that lacks it), standard
gzip is used as a fallback.

### Pulling the image manually

```bash
singularity pull kmer-gwas-cpp_v3.1.0.sif \
    docker://ghcr.io/mjfi2sb3/kmer-gwas-cpp:v3.1.0
```

### Rebuilding the image

```bash
docker build -t ghcr.io/mjfi2sb3/kmer-gwas-cpp:<tag> .
docker push ghcr.io/mjfi2sb3/kmer-gwas-cpp:<tag>
```

New images are built and pushed automatically via GitLab CI when a git tag is pushed (see `.gitlab-ci.yml`).

---

## Algorithm

### k-mer counting (Stage 1)

For each accession:

1. Streams R1 and R2 (FASTQ or FASTA, plain or gzip-compressed via zlib). Reads are dispatched to workers in batches and discarded as they are consumed — the full read set is never resident, so memory does not scale with coverage or genome size.
2. Extracts all canonical k-mers using a rolling encoder (O(1) per base). K-mer *windows* containing a non-ACGT base are skipped; the rest of the read is still used.
3. Workers count each batch and append the results into per-bin in-memory buffers, partitioned by a hash of the k-mer encoding.
4. When a bin exceeds its share of the memory budget it is sorted and run-length compacted (which at typical coverage reclaims most of the space); only if it is still over budget is a block spilled to one temporary file.
5. Each bin is finally merged, filtered by `--min_kmer_count`, and written as a sorted slice into the accession's pack file.

### Matrix construction (Stage 2)

For each bin:

1. Opens one streaming cursor per accession, seeking straight to that bin's slice inside each pack.
2. K-way merges the sorted slices through a heap, emitting each distinct k-mer once. Memory is proportional to the number of accessions, not to the number of distinct k-mers.
3. Outputs a tab-separated matrix: rows = k-mers, columns = accessions.
4. Applies the two-sided `--threshold` MAF filter on accession occurrence; respects `--count`, `--delimiter`, and `--core` flags.
5. Multi-threaded, controlled via `--threads` (set by `--matrix_merge_cpus`).

---

## FASTQ file discovery

The pipeline searches for R1/R2 files in `--data_dir` using these extensions in priority order:

```
_1.fq   _1.fastq   _1.fq.gz   _1.fastq.gz    (R1)
_2.fq   _2.fastq   _2.fq.gz   _2.fastq.gz    (R2)
```

The first matching file per read pair is used. If no file is found for an accession, the job fails immediately with an error identifying the missing accession.

---

## Tips

**Estimating `--num_bins` - Under Development**

A higher bin count reduces memory per `MATRIX_MERGE` job but increases job count. A rough starting point:

```
num_bins ≈ (num_accessions × genome_size_bp × coverage × 8 bytes) / available_RAM_per_node_bytes
```

**Resuming a run**

Nextflow caches completed work in the `work/` directory. The default (`--cleanup true`) deletes that directory once a run completes successfully, to save storage and inodes — which also means **a successfully cleaned-up run cannot be resumed**, because the cache it would resume from is gone.

To keep the option of resuming, run with `--cleanup false`. A run that *failed* always keeps its work directory regardless of the setting, so debugging a failure never requires it. Resume with:

```bash
nextflow run main.nf -profile slurm_container -resume --cleanup false \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```


**Large cohorts (inode management)**

Inode use no longer scales with `--num_bins`. KMER_COUNT writes exactly one file per accession (`<accession>.kbin`) and creates no intermediates: measured peak inode use during a job is 2, against 4,502 for the previous per-bin directory scheme at 1500 bins. MATRIX_MERGE seeks directly to the bin it needs inside each pack, so it has no extraction step at all — previously that step created 2 × `num_accessions` inodes per bin job, concurrently across every running job.

Stage 1 also publishes with `mode: 'move'` rather than `'copy'`, so a pack exists in exactly one place instead of being duplicated between the work and results directories for the duration of the run.

Work directories are removed on successful completion (`--cleanup true`, the default), further reducing inode usage.

To preserve work directories for debugging or `-resume`, pass `--cleanup false`.

For very large cohorts, also consider increasing `executor.queueSize` in `nextflow.config` if your SLURM cluster permits a higher concurrent job limit. The default is 200.

**Getting help**

```bash
nextflow run main.nf --help
```

---

## Repository structure

```
.
├── main.nf                  # Nextflow workflow entry point
├── nextflow.config          # Pipeline parameters and execution profiles
├── modules/
│   ├── kmer_count.nf        # KMER_COUNT process definition
│   └── matrix_merge.nf      # MATRIX_MERGE process definition
├── src/
│   ├── kmer_count_v3.cpp    # Stage 1: streaming k-mer counting -> pack file
│   ├── matrix_merge.cpp     # Stage 2: k-way merge of pack slices -> matrix
│   ├── kmer_key.hpp         # packed k-mer key + saturating count type
│   ├── pack_io.hpp          # self-indexed pack file format (reader/writer/validator)
│   ├── bin_store.hpp        # bounded-memory accumulation with sort/compact/spill
│   ├── thread_pool.hpp      # thread pool implementation
│   └── Makefile             # local build (`make`, `make KMER_K=31`, `make test`)
├── Dockerfile               # Container image definition
├── .gitlab-ci.yml           # CI/CD pipeline (auto-build on tag)
└── accessions.txt.example   # Example accessions file
```
