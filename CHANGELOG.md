# Changelog

All notable changes to this pipeline are documented here.
Versions follow [semantic versioning](https://semver.org): the major number is
bumped when output or interfaces change in a backwards-incompatible way.

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
- **`bits_to_text`** — expand a `bits` matrix back to text (gzip-aware, optional
  accession-name header). Two interchangeable implementations with identical
  output: `tools/bits_to_text.py` (portable, no build) and a C++ version
  (`src/bits_to_text.cpp`, `make bits_to_text`) that is ~14× faster for large
  matrices. Output is byte-identical to a natively written text matrix.
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
