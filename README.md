# kmer-GWAS Pipeline

A high-performance Nextflow pipeline for genome-wide association studies (GWAS) using k-mer presence/absence or count matrices. Designed for large-scale deployment on HPC clusters with SLURM and Singularity/Apptainer.

---

## Overview

The pipeline processes paired-end FASTQ files (plain or gzip-compressed) and produces a binary k-mer matrix across all accessions. It runs in two stages:

1. **KMER_COUNT** — one SLURM job per accession. Reads paired-end FASTQ files, extracts canonical 51-mers, writes per-accession k-mer bin files (`*_nr.bin`), and packages them into a single tar archive.
2. **MATRIX_MERGE** — one SLURM job per bin. Extracts the relevant bin from each tar archive, merges bin files across all accessions into a tabular matrix file, applies count thresholds and output format options, and compresses the result with pigz.

### Key features

- Streams gzip-compressed FASTQ directly via zlib — no temporary decompressed files on the shared filesystem
- Compiles C++ binaries at runtime with `-march=native`, targeting the actual node CPU for optimal performance
- Fully containerised: the build toolchain (GCC 12, zlib-dev, pigz) is embedded in the image; source code is compiled fresh on each job
- Works across HPC systems — no hardcoded paths or site-specific module dependencies
- Supports presence/absence or raw count output, configurable thresholds, and core k-mer extraction
- Inode-efficient: bin files are packaged into per-accession tar archives; MATRIX_MERGE extracts only the needed bin — scales to 20,000+ accessions without exhausting filesystem inodes
- Compressed output: matrix files are compressed with pigz (parallel gzip) if available, otherwise standard gzip

---

## Requirements

| Component | Minimum version |
|-----------|----------------|
| Nextflow  | 23.x (DSL2)    |
| Singularity / Apptainer | any |
| SLURM     | any            |

No compiler or zlib installation is required on compute nodes when using the `slurm_container` profile — everything is provided by the container image.

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
mkdir -p data
cp /path/to/fastqs/*.fastq.gz data/

# 4. Run on SLURM with Singularity (recommended)
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir ./data
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

---

## Output

```
results/
├── kmer_count/
│   └── <accession>.tar          # per-accession tar archive of k-mer bin files (*_nr.bin)
├── matrix/
│   └── matrix_*/
│       ├── <bin>_matrix.tsv.gz  # k-mer matrix for this bin (gzip-compressed)
│       └── <bin>_core.txt.gz    # core k-mers (if --core y, gzip-compressed)
└── reports/
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--accessions_file` | `./accessions.txt` | Path to accessions list |
| `--data_dir` | `./data` | Directory containing paired FASTQ files |
| `--output_dir` | `./results` | Output directory |
| `--num_bins` | `1500` | Number of k-mer bins (affects memory per job) |
| `--threshold` | `0` | Minimum per-accession count to include a k-mer |
| `--count` | `n` | `y` = raw counts, `n` = presence/absence |
| `--delimiter` | `tab` | Matrix column delimiter: `tab` or `none` |
| `--core` | `n` | `y` = write core k-mers file per bin |
| `--matrix_merge_cpus` | `32` | Threads for the MATRIX_MERGE stage |
| `--cleanup` | `true` | Delete Nextflow work directory on successful completion. Pass `--cleanup false` to preserve work dirs for debugging or `-resume` |
| `--clusterOptions` | _(none)_ | Extra SLURM flags passed to every job (see note below) |
| `--singularity_cache_dir` | `./.singularity` | Local path for Singularity image cache |

> **`--clusterOptions` syntax:** because the value starts with `--`, you must use the `=` form to prevent Nextflow misinterpreting it as a flag:
> ```bash
> --clusterOptions='--account=myproject --partition=highmem'
> ```

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
| KMER_COUNT | 32 | 128 GB | 8 h |
| MATRIX_MERGE | 32 (configurable) | 64 GB | 5 h |

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
    --data_dir ./data \
    --num_bins 5
```

---

## Container

The pipeline image is available on the GitHub Container Registry:

```
ghcr.io/mjfi2sb3/kmer-gwas-cpp:v2.4.2
```

The image provides GCC 12, zlib-dev, pigz, and the pipeline source code at `/opt/kmer-gwas/src/`. Binaries are compiled at job start with `-march=native` so they are optimised for the actual compute node CPU. The image does not contain pre-compiled binaries. pigz is used for parallel gzip compression of matrix output files; if pigz is unavailable (e.g. on older images or the `slurm` profile with a host that lacks pigz), standard gzip is used as a fallback.

### Pulling the image manually

```bash
singularity pull kmer-gwas-cpp_v2.4.2.sif \
    docker://ghcr.io/mjfi2sb3/kmer-gwas-cpp:v2.4.2
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

1. Reads R1 and R2 FASTQ files (plain or gzip-compressed via zlib streaming).
2. Extracts all canonical 51-mers; reads with non-ACGT bases are discarded.
3. Uses a thread pool to process read chunks in parallel, accumulating k-mer counts in per-thread hash maps.
4. Writes k-mer keys and counts to `num_bins` binary shard files, partitioned by a hash of the 51-mer bit encoding.
5. Deduplicates each shard (merges counts, filters singletons) and writes compact `*_nr.bin` files.

### Matrix construction (Stage 2)

For each bin:

1. Reads all accession `*_nr.bin` files for that bin.
2. Builds a unified k-mer index across all accessions.
3. Outputs a tab-separated matrix: rows = k-mers, columns = accessions.
4. Applies `--threshold` filter; respects `--count`, `--delimiter`, and `--core` flags.
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

**Estimating `--num_bins`**

A higher bin count reduces memory per `MATRIX_MERGE` job but increases job count. A rough starting point:

```
num_bins ≈ (num_accessions × genome_size_bp × coverage × 8 bytes) / available_RAM_per_node_bytes
```

**Resuming a run**

Nextflow caches completed work in the `work/` directory. Resume an interrupted run with:

```bash
nextflow run main.nf -profile slurm_container -resume \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

**Large cohorts (inode management)**

The pipeline is designed to handle 20,000+ accessions without exhausting filesystem inode limits. KMER_COUNT packages all bin files for each accession into a single tar archive (`<accession>.tar`); MATRIX_MERGE extracts only the one bin it needs from each archive. Work directories are removed on successful completion (`--cleanup true`, the default), further reducing inode usage.

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
│   ├── kmer_count_v3.cpp    # k-mer counting implementation
│   ├── matrix_merge.cpp     # matrix construction implementation
│   ├── mmap_io.cpp / .hpp   # memory-mapped I/O utilities
│   └── thread_pool.hpp      # thread pool implementation
├── Dockerfile               # Container image definition
├── .gitlab-ci.yml           # CI/CD pipeline (auto-build on tag)
└── accessions.txt.example   # Example accessions file
```
