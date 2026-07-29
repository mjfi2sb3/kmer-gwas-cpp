# Changelog

All notable changes to this pipeline are documented here.
Versions follow [semantic versioning](https://semver.org): the major number is
bumped when output or interfaces change in a backwards-incompatible way.

## v3.3.0

Auto-sizes the Stage 1 memory budget from the RAM actually enforced on the job.
No format or interface changes; output is byte-identical to v3.2.0.

### Changed

- **The Stage 1 accumulation budget is auto-sized from the enforced memory
  limit**, not from the config request. Previously it was `0.7 × task.memory`
  computed in the Nextflow module, which was wrong under `--mem=0` / `--exclusive`
  (where `task.memory` still reads the config number while the job owns the whole
  node) and under the `standard` profile (where `task.memory` is null, so it fell
  back to a fixed 8 GB and `--kmer_count_memory` had no effect). `kmer_count` now
  detects the limit at run time — cgroup v2 (walking `/proc/self/cgroup` up to the
  root, since the leaf often reports `max` while the real limit sits on an
  ancestor), else cgroup v1, else `SLURM_MEM_PER_NODE`/`_PER_CPU`, else
  `MemTotal` — and sets the budget to 70% of it. This tracks the real node
  everywhere: a shared allocation, an exclusive node, or a bare workstation.
- **`--kmer_count_memory` is now the scheduler request** (the SLURM `--mem` and
  the cgroup ceiling detection reads back), documented as such — no longer the
  source of the budget.
- **The `standard` (local) profile declares per-process memory**, so the local
  executor gates concurrency on RAM against its pool; otherwise several
  `KMER_COUNT` tasks could run at once and each auto-size its budget to the whole
  machine, over-committing it.

### Added

- **`--kmer_count_budget_gb`** (default `0` = auto). `0` auto-sizes as above; a
  positive value forces an explicit budget, still capped at the detected limit so
  it can never exceed what the cgroup would OOM-kill.

## v3.2.0

Parallelises the Stage 1 finalize phase. No format or interface changes; output
is byte-identical to v3.1.1.

### Changed — performance

- **Stage 1 (`KMER_COUNT`) sorts its bins in parallel.** The finalize step —
  sorting, compacting and count-filtering each of the `--num_bins` bins before
  writing the pack — previously ran one bin at a time on a single core and
  dominated Stage 1 wall-time (measured 105 of 154 s on a rice accession, with
  39 of 40 cores idle throughout). The bins are independent, so they are now
  sorted across the same thread pool that does the counting; the pack write
  stays serial (bins must be written in ascending order for the offset table)
  but is cheap next to the sort. Measured on a rice accession (219 M distinct
  k-mers, 1500 bins, 40 cores): finalize **105 s → 6 s** (parallel sort 2.9 s,
  serial write 3.3 s), total Stage 1 **154 s → 67 s**. The pack is
  byte-identical in both the in-memory and spill-to-disk paths, peak memory is
  unchanged, and there is no new dependency.

## v3.1.1

A data-loss fix and safer default memory requests. No format or interface
changes; results are byte-identical to v3.1.0.

### Fixed

- **`--cleanup true` deleted the published `.kbin` packs.** The completion
  handler ran `deleteDir()` over the work tree. `MATRIX_MERGE` stages the
  published `kmer_count_k<k>/` directory into its work dir as a *directory
  symlink*, and Groovy's `deleteDir()` follows directory symlinks — so cleanup
  recursed through the link into `results/` and removed the real packs. Work-dir
  cleanup now relies solely on Nextflow's native `cleanup` directive, which is
  symlink-safe. (`--cleanup false` was never affected, which is why the loss
  only appeared on the default setting.)

### Changed

- **Default job memory is 128 GB per stage** (was 370 GB), so jobs schedule on
  clusters with smaller nodes. Stage 1 caps k-mer accumulation at 0.7× (~90 GB;
  measured 27.6 GB on rice), and Stage 2's k-way merge needs only ~1.5–2 GB at
  12,650 accessions, so 128 GB is safe headroom. Override with
  `--kmer_count_memory` / `--matrix_merge_memory`.
- The `standard` (local) profile now requests 32 cpus / 128 GB (was 64 / 256).

## v3.1.0

Separates the matrix *encoding* from the value *delimiter*, adds a converter for
the bit-packed format, and makes the default output the most compact text form.

### Breaking changes

- **`--delimiter bits` is removed; use `--encoding bits`.** In v3.0.0 the
  bit-packed format was a value of `--delimiter`, which was a category error —
  it is an encoding, not a separator. `--delimiter bits` now errors with a
  pointer to `--encoding bits`. The output of `--encoding bits` is byte-identical
  to v3.0.0's `--delimiter bits`, so only the flag name changed.
- **The default matrix format changed** from tab-delimited to the compact form
  `--encoding text --count n --delimiter none`, i.e. rows like `KMER<TAB>101`.
  Pass `--delimiter tab` for the previous default.

### Added

- **`--encoding text|bits`** — selects how the value vector is serialised,
  independent of `--delimiter` and `--count`.
- **`--delimiter space`** — space-separated values, alongside `tab` and `none`.
- **`bits_to_text`** — convert a matrix between the bit-packed and text forms,
  both ways: `--decode` (bits → text, default) and `--encode` (text → bits, to
  compress an existing text matrix). Gzip-aware, optional accession-name header.
  Two interchangeable implementations with byte-identical output:
  `tools/bits_to_text.py` (portable, no build) and a C++ version
  (`src/bits_to_text.cpp`) that is ~14× faster for large matrices.
- Each run now **compiles `bits_to_text` into `results/bin/`** (portable build,
  no `-march=native`), so a results directory ships with the tool that reads its
  bit-packed matrices.
- The `matrix_merge` `--help`, `main.nf --help`, and the README now include the
  full `--encoding` × `--count` × `--delimiter` combination table.

### Changed

- Both illegal output combinations (`--encoding bits --count y` and
  `--delimiter none --count y`) are validated at launch in `main.nf`, not only
  in the binary, so a bad request fails immediately rather than in every job.
- A `manifest{}` block records the pipeline version; it appears in
  `run_manifest.txt`.

## v3.0.0

Bounded-memory redesign of the counting and merging stages, plus the
correctness and documentation fixes found alongside. Validated on real rice
data against v2.5.1: **396,551,209 vs 396,551,210 matrix rows, zero unexplained
differences.** Full measurements in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

This is a **major** release because the on-disk formats and the matrix output
changed; see *Breaking changes* below.

### Changed — performance and scaling

- **Stage 1 (`KMER_COUNT`) streams reads and bounds memory.** Reads are parsed
  in batches and discarded; k-mer accumulation is capped by a configurable
  budget and spills to a single temporary file only if exceeded. Peak memory is
  set by the budget rather than by the input. Measured 57.1 GB → 27.6 GB on real
  rice data, and no longer scaling with input size.
- **Stage 1 writes one self-indexed pack file per accession** (`<accession>.kbin`)
  instead of a per-bin directory tree packed into a tar. Peak filesystem entries
  (inodes) per job: 4,502 → 2 at 1500 bins. No intermediate files are written
  and re-read (~12 GB per accession, eliminated).
- **Stage 2 (`MATRIX_MERGE`) k-way merges sorted pack slices** instead of
  building a hash table over the whole bin. Memory is now proportional to the
  number of accessions, not the number of distinct k-mers: 6.4 GB → <0.1 GB in
  testing, and independent of genetic diversity. `--num_bins` is consequently a
  parallelism / output-size knob, no longer a memory control.
- **The k-mer encoder is ~460× faster** (1.14 → 531 M k-mer/s) via a rolling
  bit-operation codec, and is no longer the limiting step.

### Added

- **`--kmer_size`** — *k* is configurable (odd, 15–63), compiled in via
  `-DKMER_K` at no runtime cost. `k ≤ 32` uses a compact 10-byte record instead
  of 18, cutting Stage 1 output by ~44%. Previously fixed at 51.
- **`--delimiter bits`** — presence/absence packed one bit per accession as hex,
  ~8× smaller than `tab` at large cohorts. Presence/absence only.
- **`--min_kmer_count`** — drop k-mers seen fewer than *n* times within an
  accession (default 2). Previously hardcoded.
- **FASTA input** is parsed correctly (previously detected then parsed as FASTQ,
  silently producing garbage).
- **A run manifest** (`results/run_manifest.txt`) recording the parameters,
  pipeline version and *k* a results directory was produced with.
- **`tools/compare_to_baseline.sh`** — reproduce the validation against any
  earlier revision on your own data.
- **`src/Makefile`** — local build (`make`, `make KMER_K=31`, `make test`).

### Fixed

- **Poly-A k-mers were silently discarded.** `A` encodes as all-zero bits, which
  the old code also used as its "ambiguous base" marker; the two collided. Poly-A
  tracts now produce their k-mer correctly.
- **Spurious second tab** after the k-mer column in the matrix is removed.
- **`--delimiter none` with `--count y`** produced unparseable concatenated
  counts; the combination is now rejected.
- **`--index` was read from an uninitialised variable** when omitted from
  `matrix_merge`; it is now validated.
- **k-mer counts wrapped silently** past 65535; they now saturate.
- **The `standard` and `slurm` profiles could never run** — both compiled from a
  path that exists only inside the container image. They now compile from the
  checkout, and so does `slurm_container` (the image is used purely as a
  toolchain, so its source copy can never drift from the code being run).
- Documentation corrected: `--threshold` is a two-sided minor-allele-frequency
  filter on accession occurrence (not a per-accession count filter); core k-mers
  are excluded from the matrix by design; `--cleanup` disables `-resume`.

### Breaking changes

- **Stage 1 output is now `<accession>.kbin` pack files**, not tar archives of
  per-bin files. Intermediate files from older runs cannot be read by this
  version and must be regenerated.
- **The matrix format changed**: the spurious second tab is gone, and poly-A
  k-mers now appear as rows. Parsers depending on the old layout must be updated.
- **Output directory names now embed `k`** (`kmer_count_k<k>/`,
  `matrix_acc<N>_k<k>_...`), so results at different *k* cannot be confused.
- The container image tag referenced by the `slurm_container` profile is
  `v3.0.0`.

## v2.5.1 and earlier

See the git history. v2.x used an in-memory hash-table counting stage, tar-packed
per-bin intermediate files, and a fixed k-mer length of 51.
