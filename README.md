# kmer-GWAS

A high-performance Nextflow pipeline for generating genome-wide k-mer matrices for genome-wide association studies. It is designed for large cohorts and supports SLURM, Singularity/Apptainer, paired-end FASTQ/FASTA input or pre-built KMC databases, presence/absence or raw-count matrices, and bit-packed output.

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

1. **KMER_COUNT** processes one accession per job, streams its reads, counts canonical k-mers, and writes `<accession>.kbin`. With `--input_mode kmc`, **KMC_TO_PACK** replaces it and converts an existing KMC database into the same `<accession>.kbin`.
2. **MATRIX_MERGE** processes bins in parallel, merges the sorted slices from all accession packs, applies the optional MAF filter, and writes compressed matrix files.

### Key features

- Canonical k-mers only: a k-mer and its reverse complement are counted as one k-mer, stored as the smaller of the two (A<C<G<T order).
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
    --kmer_size 51 \
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
    --kmer_size 51 \
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

### KMC databases

If the reads are not available, Stage 1 can start from pre-built KMC databases. Use `--input_mode kmc`. Each accession needs one database in `--kmc_dir` (default: `--data_dir`):

```text
<accession>.kmc_pre
<accession>.kmc_suf
```

```bash
nextflow run main.nf \
    -profile slurm_container \
    --accessions_file accessions.txt \
    --input_mode kmc \
    --kmc_dir /path/to/kmc_dbs
```

The pipeline reads the k-mer length from the first accession's database, so `--kmer_size` is not required. If you give `--kmer_size` and it differs from the databases, the run stops. A database built with a different k-mer length fails its job.

The conversion produces the same packs as reads mode, so Stage 2 and the output are unchanged. KMC does not need to be installed; only the database files are read. The pipeline includes the KMC 3.2.4 API. It reads KMC1-format and KMC2-format databases. Databases built with `-cs1` or `kmc_tools compact` store no real counts. They give correct presence/absence matrices; with `--keep_kmer_counts y`, every count is 1. A database stores counts up to its `-cs` value. `pack_stats.tsv` gives this as `count_limit`: 65535 for a counter of 2 or more bytes, 255 for a 1-byte counter (the KMC default `-cs255`; KMC does not record a smaller `-cs`, so 255 is an upper limit), and 1 for a database with no counter. `kmers_at_count_limit` counts the k-mers that reached it, and the coverage estimate leaves them out.

- Both modes produce canonical k-mers only. A `-b` (non-canonical) database is converted to canonical k-mers, and the counts of a k-mer and its reverse complement are added together.
- Each database is kept as it is: every k-mer in it is used, and the `-ci` used to build it is the count filter. `--min_kmer_count` is applied only where it is larger than a database's `-ci`. If a database's counts never reach it (for example, a database built with `-cs1`), the job stops.
- A `-b` database built with `-ci` greater than 1 can be missing k-mers whose count was split between the two strands. The job writes a warning.
- Before any conversion, the pipeline reads every database header and writes `kmc_databases.tsv`. The run stops if the databases differ in k-mer length, `-ci` or `-cx` (or in counter size with `--keep_kmer_counts y`). `--kmc_allow_mixed true` turns this into a warning. A `--min_kmer_count` equal to or larger than every database's `-ci` makes different `-ci` values consistent.

## Output

```text
results/
├── run_manifest.txt
├── pack_stats.tsv
├── kmc_databases.tsv          (kmc mode)
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

Each pack records its metadata: accession, source (`reads` or `kmc`), the count filter applied, the pipeline and KMC API versions, and statistics (distinct k-mers, total k-mers, largest count, estimated k-mer coverage, and, in reads mode, reads, bases and skipped windows). `kbin_dump --info` prints it. `pack_stats.tsv` lists the metadata of every pack, and the log shows a short summary when Stage 1 finishes. Stage 2 stops if a pack's recorded accession differs from its file name.

`run_manifest.txt` records how the results directory was made: the parameters, the command line, the profile, the container image, and a `code_checksum`. The checksum is one MD5 over `main.nf`, `nextflow.config`, `modules/` and the C++ sources, so two runs used the same code only if their checksums are the same. The manifest is written at launch with `status = running`, and a completion section (status, duration, CPU hours, task counts and, on failure, the error) is added when the run ends.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--accessions_file` | `./accessions.txt` | Accession list. |
| `--input_mode` | `reads` | Stage 1 input: `reads` (FASTQ/FASTA) or `kmc` (KMC databases). |
| `--data_dir` | `./data` | Directory containing read files. |
| `--kmc_dir` | `--data_dir` | Directory containing KMC databases (`--input_mode kmc`). |
| `--kmc_allow_mixed` | `false` | kmc mode: warn instead of stopping when the databases were built with different `-ci`, `-cx` or counter sizes. |
| `--skip_accession_check` | `false` | Allow packs whose recorded accession differs from their file name. |
| `--output_dir` | `./results` | Output directory. |
| `--kmer_size` | none | Odd k-mer length from 15 to 63. Required with `--input_mode reads`. With `--input_mode kmc`, read from the databases. |
| `--num_bins` | `1500` | Number of k-mer bins and Stage 2 output shards. |
| `--min_kmer_count` | `2` (reads mode) | Discard k-mers observed fewer than this many times within an accession. In kmc mode, the databases are kept as they are (their own `-ci` applies); a larger value is enforced as well. |
| `--threshold` | `0` | Two-sided minor-allele-frequency filter expressed as an accession-count threshold. |
| `--keep_kmer_counts` | `n` | `n` for presence/absence; `y` for raw per-accession counts. |
| `--encoding` | `text` | `text` or `bits`; bit encoding is for presence/absence only. |
| `--delimiter` | `none` | Text value delimiter: `none`, `tab`, or `space`. |
| `--write_core_kmers` | `false` | Write k-mers present in every accession to a separate file. |
| `--matrix_merge_cpus` | `16` | Threads used by MATRIX_MERGE. |
| `--kmer_count_memory` | `128.GB` | Scheduler memory request per KMER_COUNT and KMC_TO_PACK job. |
| `--kmer_count_budget_gb` | `0` | Stage 1 accumulation budget. `0` uses 70% of the detected enforced memory limit. |
| `--kmer_count_read_threads` | `2` | Concurrent read/decompression threads. |
| `--matrix_merge_memory` | `16.GB` | Memory request per MATRIX_MERGE job. |
| `--kmer_count_time` | `5h` | Wall-clock limit per KMER_COUNT and KMC_TO_PACK job. |
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

`--kmer_size` must be odd and between 15 and 63. It is required in reads mode. In KMC mode, it is read from the databases. It is compiled into the C++ programs at the beginning of each job. Values up to 32 fit in one 64-bit word, reducing each Stage 1 record from 18 to 10 bytes. Changing `k` invalidates existing bin files; use a new `--output_dir`.

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

The k-mer in each row is canonical. The k-mer is separated from the values by a tab. `--delimiter` controls the separator between values in text encoding.

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
    --kmer_size 51 \
    --data_dir /path/to/fastq

# Resume after a failure or interruption
nextflow run main.nf \
    -profile slurm_container \
    -resume \
    --cleanup false \
    --kmer_count_publish_mode link \
    --matrix_publish_mode link \
    --accessions_file accessions.txt \
    --kmer_size 51 \
    --data_dir /path/to/fastq
```

## Execution profiles

### `slurm_container`

The recommended SLURM profile. It uses a pre-built Singularity/Apptainer image and does not require compiler or zlib modules on compute nodes.

| Stage | CPUs | Memory | Time |
|---|---:|---:|---:|
| KMER_COUNT | 32 | `128.GB` | `5h` |
| KMC_TO_PACK | 2 | `128.GB` | `5h` |
| MATRIX_MERGE | 16 | `16.GB` | `10h` |

### `slurm`

Uses the host environment on SLURM and requires GCC 12 or later, zlib development files, and pigz where available.

```bash
nextflow run main.nf \
    -profile slurm \
    --accessions_file accessions.txt \
    --kmer_size 51 \
    --data_dir /path/to/fastq
```

### `standard`

Runs locally and is intended for small-scale testing.

```bash
nextflow run main.nf \
    -profile standard \
    --accessions_file accessions.txt \
    --kmer_size 51 \
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
│   ├── kmc_to_pack.nf
│   ├── matrix_merge.nf
│   └── build_tools.nf
├── src/
│   ├── kmer_count_v3.cpp
│   ├── kmc_to_pack.cpp
│   ├── kmc_api/            # KMC 3 API (GPL-3)
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

## Licence

The pipeline is released under the MIT licence (see `LICENSE`). The exception is
`src/kmc_api/` (the KMC 3.2.4 API, unmodified) and `src/kmc_to_pack.cpp`, which
links it: these files, and the `kmc_to_pack` binary built from them, are under
the GNU GPL 3. See `src/kmc_api/README`.

## Citation and support

If this pipeline contributes to published work, cite the associated software release and report the pipeline version, input mode (reads or KMC databases), k-mer length, `--threshold` MAF setting, `--min_kmer_count` (or the KMC databases' `-ci`), matrix encoding, execution profile, and relevant hardware or cluster configuration. `run_manifest.txt` records these values.

For usage information:

```bash
nextflow run main.nf --help
```

Please report reproducible issues with the command, profile, parameter values, input layout, Nextflow version, and relevant log or trace files.