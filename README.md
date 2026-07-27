# kmer-GWAS Pipeline

A high-performance Nextflow pipeline for genome-wide association studies (GWAS) using k-mer presence/absence or count matrices. Designed for large-scale deployment on HPC clusters with SLURM and Singularity/Apptainer.

---

## What's new

The counting and merging stages have been rewritten so that memory and inode use
are bounded by configuration rather than by the data. Validated on real rice
data against the previous implementation: **396.5 million matrix rows, zero
unexplained differences.**

| | before | after |
|---|---|---|
| `KMER_COUNT` peak memory | 57.1 GB, scales with input | **27.6 GB, capped by a budget** |
| `MATRIX_MERGE` peak memory | 6.4 GB, scales with k-mer union | **< 0.1 GB, scales with accession count** |
| `KMER_COUNT` peak inodes (1500 bins) | 4,502 | **2** |
| `MATRIX_MERGE` inodes per bin job | 25,201 | **0** |
| intermediate I/O per accession | ~12 GB written and re-read | **none** |
| k-mer encoder | 1.14 M k-mer/s | **531 M k-mer/s** |

Practical consequences:

- **Peak memory no longer tracks input size.** The old code held roughly three
  copies of every read before counting anything; reads are now streamed and
  k-mer accumulation is capped by a configurable budget.
- **`--num_bins` no longer has to grow with genome size**, because Stage 2 memory
  no longer depends on the number of distinct k-mers — and bin count was what
  drove inode use.
- **`--kmer_size` is configurable** (odd, 15–63) at no runtime cost. `k ≤ 32`
  uses a 10-byte record instead of 18, cutting Stage 1 output by ~44%.
- **`--delimiter bits`** packs presence/absence one bit per accession, ~8×
  smaller than `tab` at 12,600 accessions.

Two intentional output changes: a spurious second tab after the k-mer column is
gone, and poly-A k-mers are no longer silently discarded. Everything else is
identical to the previous implementation.

Full measurements, methodology and test coverage: **[docs/BENCHMARKS.md](docs/BENCHMARKS.md)**.
Reproduce the comparison yourself with `tools/compare_to_baseline.sh`.

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
| `--delimiter` | `tab` | Matrix value format: `tab`, `none` or `bits`. At large cohorts `bits` is ~8× smaller — see note below |
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

> **Matrix row format and cohort size.** Each row is the k-mer, a tab, then one value per accession, so row width grows linearly with the panel:
>
> | accessions | `tab` | `none` | `bits` |
> |---|---|---|---|
> | 1,000 | 2,051 B | 1,052 B | 302 B |
> | 12,600 | 25,251 B | 12,652 B | **3,202 B** |
>
> `bits` packs presence/absence one bit per accession and writes it as hex, making it about 8× smaller than `tab` at 12,600 accessions. Because the matrix usually dwarfs every intermediate file in this pipeline, that is the single largest lever on total output size. It is presence/absence only and cannot be combined with `--count y`.
>
> Decoding a `bits` row: split on the tab, `bytes.fromhex(rest)`, then accession *i* is present if `(bs[i >> 3] >> (i & 7)) & 1`.
>
> `none` concatenates single characters with no separator; it is only meaningful with `--count n`, since multi-digit counts would be unparseable, and that combination is rejected.

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
ghcr.io/mjfi2sb3/kmer-gwas-cpp:v2.5.0
```

The image provides GCC 12, zlib-dev, pigz, and the pipeline source code at `/opt/kmer-gwas/src/`. Binaries are compiled at job start with `-march=native` so they are optimised for the actual compute node CPU. The image does not contain pre-compiled binaries. pigz is used for parallel gzip compression of matrix output files; if pigz is unavailable (e.g. on older images or the `slurm` profile with a host that lacks pigz), standard gzip is used as a fallback.

### Pulling the image manually

```bash
singularity pull kmer-gwas-cpp_v2.5.0.sif \
    docker://ghcr.io/mjfi2sb3/kmer-gwas-cpp:v2.5.0
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
