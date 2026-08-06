# kmer-GWAS

A high-performance Nextflow pipeline for generating genome-wide k-mer matrices for genome-wide association studies. It is designed for large cohorts and supports SLURM, Singularity/Apptainer, paired-end FASTQ/FASTA input, presence/absence or raw-count matrices, and bit-packed output.

> **Status:** Benchmark values refer to the rice dataset described in [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md). Performance is hardware- and dataset-dependent.

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Input](#input)
- [Output](#output)
- [Parameters](#parameters)
- [Resuming a run](#resuming-a-run)
- [Execution profiles](#execution-profiles)
- [Matrix formats](#matrix-formats)
- [Utilities](#utilities)
- [Algorithm](#algorithm)
- [Performance guidance](#performance-guidance)
- [Repository structure](#repository-structure)
- [Citation and support](#citation-and-support)

## Overview

For each accession, the pipeline streams sequencing reads, extracts canonical *k*-mers, and writes a self-indexed binary pack containing sorted k-mer slices. A second stage seeks directly to the corresponding slice in every pack and performs a k-way merge to produce the matrix.

The output contains one row per distinct k-mer and one column per accession. It can encode either presence/absence or per-accession k-mer counts for downstream association analysis.

### Pipeline stages

1. **KMER_COUNT** processes one accession per job, streams its reads, counts canonical k-mers, and writes `<accession>.kbin`.
2. **MATRIX_MERGE** processes bins in parallel, merges the sorted slices from all accession packs, applies the optional MAF filter, and writes compressed matrix files.

### Key features

- Streaming input without loading complete read sets into memory.
- Configurable Stage 1 memory budgeting with spill-to-disk support.
- One self-indexed pack per accession instead of per-bin directory trees.
- Stage 2 memory proportional primarily to the number of accessions.
- Runtime compilation with `-march=native` for the executing node.
- Singularity/Apptainer support with a containerised GCC, zlib, and pigz toolchain.
- Presence/absence, raw-count, and bit-packed matrix output.
- Optional extraction of k-mers present in every accession.
- Nextflow caching and `-resume` support.

### Benchmark and compatibility

On the rice benchmark, the new implementation produced **396,551,210 matrix rows**, compared with **396,551,209** for the previous implementation. The single additional row is the valid all-`A` k-mer, which the previous implementation incorrectly discarded because its nucleotide encoding conflated `A` with the ambiguous-base sentinel used for `N`. The difference is fully explained and reflects a correction in the new implementation.

K-mers composed entirely of `A` bases are valid and are retained. Only windows containing an actual non-ACGT character are skipped. Other benchmark measurements include 27.6 GB Stage 1 peak memory, less than 0.1 GB Stage 2 peak memory, two Stage 1 filesystem entries per job, and 531 million k-mers per second encoder throughput.

## Requirements

| Component | Requirement |
|---|---|
| Nextflow | 23.x or later, DSL2 |
| Singularity or Apptainer | Required for `slurm_container` |
| SLURM | Required for `slurm_container` and `slurm` |
| C++17 compiler and zlib development files | Required for the host-based `slurm` profile and local builds |

The `slurm_container` profile supplies the compiler and libraries through the container.

## Quick start

```bash
git clone https://github.com/mjfi2sb3/kmer-gwas-cpp.git
cd kmer-gwas-cpp
cp accessions.txt.example accessions.txt
# Edit accessions.txt: one accession ID per line.
```

Run the recommended profile:

```bash
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

Results are written to `./results` by default. Display all options with:

```bash
nextflow run main.nf --help
```

For clusters requiring an account or partition:

```bash
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq \
    --clusterOptions='--account=myproject --partition=highmem'
```

## Input

### Accession file

The accession file is plain text with one identifier per line. Blank lines are ignored, and identifiers must match the prefixes of the read files.

```text
SAMPLE1_SRR0000001
SAMPLE2_SRR0000002
SAMPLE3_SRR0000003
```

### Read files

Paired-end files must follow this pattern:

```text
<accession>_1.fq[.gz]       # R1
<accession>_2.fq[.gz]       # R2
```

The pipeline supports `.fq`, `.fastq`, `.fq.gz`, and `.fastq.gz`. Compression is detected from file content. FASTA input, including single- and multi-line records, is also supported.

The first matching pair is selected in this order:

```text
_1.fq   _1.fastq   _1.fq.gz   _1.fastq.gz
_2.fq   _2.fastq   _2.fq.gz   _2.fastq.gz
```

If either mate is missing, the relevant accession job fails with an error identifying the accession.

## Output

```text
results/
├── run_manifest.txt
├── bin/
│   ├── bits_to_text
│   └── kbin_dump
├── kmer_count_k<k>/
│   └── <accession>.kbin
├── matrix/
│   └── matrix_acc<N>_k<k>_bins<B>_minOcc<t>_<mode>_delim-<d>/
│       ├── <bin>_matrix.tsv.gz
│       └── <bin>_core.txt.gz
└── reports/
    ├── timeline.html
    ├── report.html
    └── dag.svg
```

The core-k-mer file is produced only when `--write_core_kmers true` is enabled. The k-mer length appears in output directory names so results produced with different values of `k` are not accidentally mixed.

Each `.kbin` file contains sorted bin slices and an index for direct seeking. Packs are uncompressed during processing and can optionally be compressed after both stages with `--compress_kbin_packs true`.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--accessions_file` | `./accessions.txt` | Accession list. |
| `--data_dir` | `./data` | Directory containing read files. |
| `--output_dir` | `./results` | Output directory. |
| `--kmer_size` | `51` | Odd k-mer length from 15 to 63. Values up to 32 use the compact key representation. |
| `--num_bins` | `1500` | Number of k-mer bins and Stage 2 output shards. |
| `--min_kmer_count` | `2` | Discard k-mers observed fewer than this many times within an accession. |
| `--threshold` | `0` | Two-sided minor-allele-frequency filter expressed as an accession-count threshold. |
| `--keep_kmer_counts` | `n` | `n` for presence/absence; `y` for raw per-accession counts. |
| `--encoding` | `text` | `text` or `bits`; bit encoding is for presence/absence only. |
| `--delimiter` | `none` | Text value delimiter: `none`, `tab`, or `space`. |
| `--write_core_kmers` | `false` | Write k-mers present in every accession to a separate file. |
| `--matrix_merge_cpus` | `16` | Threads used by MATRIX_MERGE. |
| `--kmer_count_memory` | `128.GB` | Scheduler memory request per KMER_COUNT job. |
| `--kmer_count_budget_gb` | `0` | Stage 1 accumulation budget. `0` uses 70% of the detected enforced memory limit. |
| `--kmer_count_read_threads` | `2` | Concurrent read/decompression threads. |
| `--matrix_merge_memory` | `16.GB` | Memory request per MATRIX_MERGE job. |
| `--kmer_count_time` | `5h` | Wall-clock limit per KMER_COUNT job. |
| `--matrix_merge_time` | `10h` | Wall-clock limit per MATRIX_MERGE job. |
| `--cleanup` | `false` | Remove the Nextflow work directory after successful completion. |
| `--kmer_count_publish_mode` | `link` | Publish Stage 1 packs using `link`, `copy`, or `move`. |
| `--matrix_publish_mode` | `link` | Publish Stage 2 matrices using `link`, `copy`, or `move`. |
| `--compress_kbin_packs` | `false` | Compress published `.kbin` files after both stages finish. |
| `--clusterOptions` | none | Additional SLURM options passed to jobs. |
| `--queue_size` | `200` | Maximum queued and running jobs submitted by SLURM. |
| `--singularity_cache_dir` | `./.singularity` | Singularity/Apptainer image cache directory. |

Because `--clusterOptions` begins with a hyphen, use the equals form:

```bash
--clusterOptions='--account=myproject --partition=highmem'
```

### K-mer length

`--kmer_size` must be odd and between 15 and 63. It is compiled into the C++ programs at the beginning of each job. Values up to 32 fit in one 64-bit word, reducing each Stage 1 record from 18 to 10 bytes. Changing `k` invalidates existing bin files; use a new `--output_dir`.

### MAF threshold

`--threshold` is a **two-sided minor-allele-frequency (MAF) filter expressed as an accession-count threshold**. It filters on the number of accessions carrying a k-mer, not on the number of times that k-mer occurs within an individual accession.

For a k-mer present in `c` of `N` accessions, it is retained when:

```text
threshold ≤ c ≤ N − threshold
```

The corresponding binary k-mer MAF is:

```text
min(c / N, 1 − c / N)
```

Thus, the filter requires MAF to be at least `threshold / N`. For example, with 100 accessions, `--threshold 20` retains k-mers present in at least 20 and at most 80 accessions. The default, `--threshold 0`, disables the filter. The within-accession count filter is controlled separately by `--min_kmer_count`.

### Matrix encoding

The k-mer is separated from the values by a tab. `--delimiter` controls the separator between values in text encoding.

| Encoding | Counts | Delimiter | Example |
|---|---|---|---|
| `text` | `n` | `none` | `KMER<TAB>101` |
| `text` | `n` | `tab` | `KMER<TAB>1<TAB>0<TAB>1` |
| `text` | `y` | `tab` or `space` | `KMER<TAB>5<TAB>0<TAB>3` |
| `bits` | `n` | ignored | `KMER<TAB>05` in hexadecimal |

`bits` and `none` require presence/absence output. Bit-packed rows store one bit per accession and are approximately eight times smaller than tab-delimited presence/absence rows at large cohort sizes.

## Resuming a run

Nextflow can reuse completed tasks only while the relevant task cache and work-directory files remain available. The two publishing parameters control this independently for the two pipeline stages.

The recommended settings are `--kmer_count_publish_mode link` and `--matrix_publish_mode link`, provided that `--output_dir` and the Nextflow work directory are on the same filesystem. Hard links avoid unnecessary data copies while retaining the work-directory files required for resumption.

Use `--kmer_count_publish_mode copy` and/or `--matrix_publish_mode copy` when the output and work directories are on different filesystems. This preserves resumability for the corresponding stage but requires an additional copy of the published data.

Use `--kmer_count_publish_mode move` and/or `--matrix_publish_mode move` to minimise storage duplication when resumption is not required. However, `move` removes the corresponding published output from the task work directory: `--kmer_count_publish_mode move` prevents reuse of completed KMER_COUNT tasks, while `--matrix_publish_mode move` prevents reuse of completed MATRIX_MERGE tasks.

To retain both caches required by `-resume`, use `--cleanup false` together with `link` or `copy` for both publishing modes. A successful run with `--cleanup true` cannot be resumed because the Nextflow work directory has been removed.

Example:

```bash
# Initial run
nextflow run main.nf \
    -profile slurm_container \
    --cleanup false \
    --kmer_count_publish_mode link \
    --matrix_publish_mode link \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq

# Resume after a failure or interruption
nextflow run main.nf \
    -profile slurm_container \
    -resume \
    --cleanup false \
    --kmer_count_publish_mode link \
    --matrix_publish_mode link \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

## Execution profiles

### `slurm_container`

The recommended SLURM profile. It uses a pre-built Singularity/Apptainer image and does not require compiler or zlib modules on compute nodes.

| Stage | CPUs | Memory | Time |
|---|---:|---:|---:|
| KMER_COUNT | 32 | `128.GB` | `5h` |
| MATRIX_MERGE | 16 | `16.GB` | `10h` |

### `slurm`

Uses the host environment on SLURM and requires GCC 12 or later, zlib development files, and pigz where available.

```bash
nextflow run main.nf \
    -profile slurm \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq
```

### `standard`

Runs locally and is intended for small-scale testing.

```bash
nextflow run main.nf \
    -profile standard \
    --accessions_file accessions.txt \
    --data_dir /path/to/fastq \
    --num_bins 5
```

The container image is published at:

```text
ghcr.io/mjfi2sb3/kmer-gwas-cpp:v3.7.5 # check for correct version
```

It supplies GCC 12, zlib, and pigz as a build toolchain. The source in the checkout is compiled inside the container at job start with `-march=native`.

## Matrix formats and utilities

Convert between text and bit-packed representations with `bits_to_text`:

```bash
results/bin/bits_to_text --encode -i input_matrix.tsv -o output_bits.tsv
results/bin/bits_to_text --decode -n 10000 -i input_bits.tsv -o output_text.tsv
results/bin/bits_to_text --decode -n 10000 -i input_bits.tsv | pigz > output_text.tsv.gz
```

`kbin_dump` prints pack metadata or exports sorted k-mers:

```bash
results/bin/kbin_dump results/kmer_count_k51/SAMPLE1.kbin --info
results/bin/kbin_dump results/kmer_count_k51/SAMPLE1.kbin --bin 42
results/bin/kbin_dump results/kmer_count_k51/SAMPLE1.kbin --all_bins ./export
```

It requires an uncompressed pack. Decompress `.kbin.gz` files first:

```bash
pigz -dc SAMPLE1.kbin.gz > SAMPLE1.kbin
```

## Algorithm

### Stage 1: k-mer counting

For each accession, the pipeline streams paired FASTQ or FASTA reads, extracts canonical k-mers with a rolling encoder, and skips only windows containing non-ACGT characters. Valid all-`A` windows are retained. K-mers are partitioned into bins, accumulated within the configured memory budget, sorted and compacted, spilled when necessary, filtered by `--min_kmer_count`, and written as sorted slices to one `.kbin` pack.

### Stage 2: matrix construction

For each bin, the pipeline opens one streaming cursor per accession, seeks directly to the relevant slice, performs a heap-based k-way merge, emits one row per distinct k-mer, applies the optional two-sided MAF threshold, and formats values as text or hexadecimal bit-packed presence/absence. Core k-mers can be written as a separate output.

## Performance guidance

`--num_bins` primarily controls parallelism and output sharding. More bins reduce the amount of data handled by each MATRIX_MERGE job but increase the number of jobs and output files.

For large cohorts:

- Set `--kmer_count_budget_gb` explicitly when a predictable Stage 1 memory ceiling is required.
- Use `--kmer_count_publish_mode link` and `--matrix_publish_mode link` when output and work share a filesystem.
- Use `--encoding bits` when downstream tools support bit-packed data or conversion through `bits_to_text`.
- Retain `--cleanup false` until validation and any required resumption are complete.
- Increase `--queue_size` only if cluster policy permits more concurrent submissions.

## Repository structure

```text
.
├── main.nf
├── nextflow.config
├── modules/
│   ├── kmer_count.nf
│   ├── matrix_merge.nf
│   └── build_tools.nf
├── src/
│   ├── kmer_count_v3.cpp
│   ├── matrix_merge.cpp
│   ├── bits_to_text.cpp
│   ├── kmer_key.hpp
│   ├── pack_io.hpp
│   ├── bin_store.hpp
│   ├── thread_pool.hpp
│   └── Makefile
├── tools/
│   ├── kbin_dump
│   ├── bits_to_text
│   └── README.md
├── Dockerfile
├── .gitlab-ci.yml
└── accessions.txt.example
```

## Citation and support

If this pipeline contributes to published work, cite the associated software release and report the pipeline version, k-mer length, `--threshold` MAF setting, `--min_kmer_count`, matrix encoding, execution profile, and relevant hardware or cluster configuration.

For usage information:

```bash
nextflow run main.nf --help
```

Please report reproducible issues with the command, profile, parameter values, input layout, Nextflow version, and relevant log or trace files.