# Changelog

All notable changes to this pipeline are documented here.
Versions follow [semantic versioning](https://semver.org): the major number is
bumped when output or interfaces change in a backwards-incompatible way.

## v3.7.5

Fixes `-resume` for Stage 2, and renames the count flag. No output change.

### Fixed

- **`-resume` now caches `MATRIX_MERGE`.** Its input was the published pack
  directory (`output_dir/kmer_count_k<k>/`), which is mutated after a run —
  re-publish timestamps, and `--compress_kbin_packs` rewriting each `.kbin` to
  `.kbin.gz` in place — so the task hash changed on every resume and the whole
  merge re-ran (and, once packs were gzipped, could not find them). The directory
  is now passed as a value (a stable path string the task reads directly) rather
  than a hashed `path` input. Correct invalidation is preserved: `k` is in the
  path, the accessions file stays a tracked input, and the Stage-2 params plus
  `min_kmer_count` are folded into the task's script hash.

### Changed

- Renamed **`--count` → `--keep_kmer_counts`** (the bare name was easily confused
  with `--min_kmer_count`). `--count` is kept as a **deprecated alias** (a launch
  warning notes the rename); it will be removed in a future release.

## v3.7.4

Separate publish mode for the Stage 2 matrix, and a codebase-wide comment
cleanup for public readers. No output or matrix format changes.

### Added

- **`kbin_dump`** now detects a gzip-compressed pack (`.kbin.gz`, e.g. left by
  `--compress_kbin_packs`) and exits with an actionable message telling you to
  decompress it first (`pigz -dc …`), instead of failing with a cryptic "bad
  magic". `kbin_dump` seeks by the footer index, so it needs the uncompressed pack.
- **`--matrix_publish_mode`** (default `link`) — controls how the Stage 2 matrix
  reaches `output_dir`, with the same `link`/`copy`/`move` choices as the Stage 1
  packs. The matrix previously always used `copy`. `link` is safe for it: the
  matrix is already gzipped and never re-compressed, so nothing breaks the hard
  link, and it avoids duplicating one of the largest outputs of a big run. A hard
  link survives `--cleanup`. Use `copy` for a fully independent archival copy.

### Changed

- Renamed **`--publish_mode` → `--kmer_count_publish_mode`** so the two publish
  modes name their stages unambiguously. `--publish_mode` is kept as a
  **deprecated alias** (a launch warning notes the rename); it will be removed in
  a future release.
- The same-filesystem launch check for `link` now covers either publish mode.
- Documentation-only comment cleanup across the source and modules: removed
  internal changelog/task references and stale magic numbers, so comments
  describe the current behaviour rather than its development history.

## v3.7.3

Overhauls the two helper tools and adds optional pack compression. No matrix
output changes.

### Added

- **`--compress_kbin_packs`** (default `true`) — after **both** stages finish
  successfully, gzip the published `.kbin` packs, one job per accession, gated on
  all of Stage 2 so a pack is never compressed while it might still be read. The
  gain is modest (~1.5–2×: the 2-bit-packed keys barely compress). Re-running
  Stage 2 needs the packs decompressed first. With `--publish_mode link` +
  `--cleanup false` the uncompressed pack survives in the work dir, so nothing is
  reclaimed until it is removed — a launch warning notes this.

### Changed

- **`kbin_dump`**: removed `--with-bin`; added **`--all_bins DIR`**, which writes
  one file per bin into `DIR/<accession>/<N>.tsv.gz`, decoding across all cores and
  gzipping in-process (no `pigz` pipe). Memory is bounded (~tens of MB) via a
  streamed reader regardless of bin size. `--no-gzip` and `--threads` opt-outs.
- **`bits_to_text`**: `-i/-o` flags; output is now **plain by default** (its gzip
  is single-threaded — `--gzip` opt-in, or pipe to `pigz`); `--encode` infers the
  accession count from the row and **auto-detects the delimiter** (default `none`),
  fixing the common "row has 1 value" error on a default (`--delimiter none`)
  matrix; `--header` removed. `--decode` still needs `-a/-n` to trim byte padding.
- **Removed `tools/bits_to_text.py`** — the C++ binary is authoritative.
- Corrected README: the `--cleanup` default is `false` (was still documented as
  `true` in places), plus the tool sections above.

## v3.7.2

Renames the `--core` flag and corrects its documentation. No output changes, and
`--core` keeps working as a deprecated alias, so upgrading is safe.

### Changed

- **`--core` renamed to `--write_core_kmers`** (a boolean, default `false`). The
  flag writes a per-bin `<bin>_core.txt` listing k-mers present in every accession.
  It is independent of every other flag: it only adds that file and never changes
  the matrix. The old `--core y` / `--core n` still works but prints a deprecation
  warning and will be removed in a later release.

### Fixed

- **Corrected the core-k-mers documentation.** The help, README, and code comments
  previously said core k-mers are "excluded from the matrix," implying this flag
  removes them. It does not: whether all-accession k-mers appear in the matrix is
  decided solely by `--threshold` (they are dropped only at `--threshold >= 1`; at
  the default `0` they remain). The flag only controls the side file.

## v3.7.1

Fixes a container bind for symlinked reads, enriches the launch banner, and tidies
code comments for publication. No output or interface changes.

### Fixed

- **Symlinked reads whose target sits under a symlinked prefix (e.g. `/project`
  itself a symlink to `/lustre2/project`) now resolve inside the container.** The
  `slurm_container` profile bound each `data_dir` symlink's target directory only
  in its fully-resolved (canonical) form, e.g. `/lustre2/project/.../data`. Inside
  the container the symlink is followed via the path it literally stores,
  `/project/.../data`, whose prefix is absent there, so it dangled and Stage 1
  failed with "no R1 file found" even though the file existed. The profile now
  binds the target directory in **both** forms — canonical and as the symlink
  stores it — so the read resolves whichever prefix the link references.
  (Workaround on older versions: repoint the data symlinks at their canonical
  target with `ln -sfn "$(readlink -f "$f")" "$f"`.)

### Changed

- **`matrix_merge` logs its merge plan before merging.** It prints the shard and
  thread count, or the reason it ran single-threaded (`--threads 1`, or an input
  below the sharding threshold), so the effective Stage 2 parallelism is visible
  in the task log up front rather than only inferable from the end summary.
- **The launch banner now reports what a run is actually configured to do.** It
  adds the pipeline version, the active profile, and the work directory (runtime
  facts, not parameters), plus previously-omitted parameters `publish_mode`,
  `kmer_count_budget_gb` (shown as `0 (auto)`), `kmer_count_read_threads`,
  `queue_size`, and the robustness ceilings (`max_retries`,
  `max_memory`/`max_cpus`/`max_time`). Fields are regrouped to mirror `--help`.
- **Code comments genericised for public release.** Removed internal breadcrumbs
  (audit task numbers, a commit hash, a renamed-file note) and the project's
  specific accession count, restating costs as the per-accession rate they
  describe. Comment-only; no behaviour change.

## v3.7.0

Parallelises the Stage 2 merge. Output is byte-identical to previous releases.

### Changed

- **The Stage 2 matrix merge now runs in parallel.** Previously the k-way merge
  was single-threaded and the requested CPUs only fed pigz. It now range-shards
  the bin's key space into `S = 2^B` contiguous slices and merges them
  concurrently. The shard is the top bits of the k-mer key; because the first
  base sits in the key's highest bits and A<C<G<T maps to 0<1<2<3, those bits are
  the sort prefix, so every accession's copy of a k-mer lands in the same shard
  and each shard is an independent, contiguous sub-range of every sorted pack.
  The per-shard outputs are concatenated in shard order, so the result is
  **byte-identical** to a single-threaded run (verified by md5 across encodings,
  the MAF and core filters, and both single-threaded and parallel paths on real
  data). Measured ~4.4× faster at 8 threads on a 62.7 M-k-mer bin.
- **`--matrix_merge_cpus` default 6 → 16.** The CPUs now accelerate the merge
  itself, not only pigz. Returns diminish past ~8–12 threads (pigz also saturates
  around 8), so 16 is a sensible ceiling; higher mainly costs RAM.
- **`--matrix_merge_memory` note updated.** Merge memory is now
  `O(matrix_merge_cpus × accessions)` — each concurrent shard worker keeps a small
  (~36 KB) read buffer per accession — but stays independent of genome size. At
  the default 16 CPUs this is ~7 GB for a 12,600-accession panel, so the 16 GB
  default is unchanged and still has headroom (a shortfall is caught by the retry
  escalation from v3.6.0).

## v3.6.0

Adds failure-recovery for large runs, right-sizes the Stage 2 job, and tidies the
tools and docs. No output format changes.

### Added

- **Automatic retry with resource escalation (SLURM profiles).** A job killed by
  an out-of-memory or timeout signal is retried with *more* memory and wallclock
  each attempt (`attempt ×` the base), so the occasional straggler in a large
  cohort self-heals instead of failing the whole run. New parameters:
  `--max_retries` (default 2), and `--max_memory` / `--max_cpus` / `--max_time`
  (defaults 512 GB / 64 / 24h) which set `resourceLimits`, clamping every request
  — including the escalated retries — to what a node can provide. Set the three
  ceilings to your largest node.

### Changed

- **`--matrix_merge_cpus` default 32 → 6.** The merge itself is single-threaded;
  the requested CPUs only compress that bin's output with pigz, which saturates
  around 8 threads. Requesting 32 left cores idle and, on a shared cluster,
  reduced how many bins could run at once.
- **`--matrix_merge_memory` default 128 GB → 16 GB.** The merge's memory is
  proportional to the number of accessions (~76 KB each), independent of genome
  size, so 16 GB covers well over 100k accessions for any genome; a rare shortfall
  is caught by the retry escalation above.
- **`--help` is reorganised** into labelled groups (input/output, k-mers & matrix,
  Stage 1 resources, Stage 2 resources, scheduler, work dir & publishing).
- **The `tools/` folder is flattened.** The prebuilt `kbin_dump` / `bits_to_text`
  binaries moved from `tools/bin/` up into `tools/`, and the obsolete
  `compare_to_baseline.sh` was removed.
- The README now has a table of contents.

## v3.5.1

Default and documentation fixes.

### Changed

- **`--cleanup` now defaults to `false`.** The work directory is kept by default
  so `-resume` can skip finished accessions out of the box (with the default
  `--publish_mode link` the packs it keeps are hard links, so this costs no extra
  data storage). Pass `--cleanup true` to delete it on success, which disables
  `-resume`.

### Fixed

- The `--help` text showed a stale `370.GB` default for `--kmer_count_memory` and
  `--matrix_merge_memory` (the real default is 128 GB), had a duplicated
  `--encoding bits` line, and omitted `--kmer_count_budget_gb` and
  `--kmer_count_read_threads`. Corrected, and the memory defaults now read from
  the actual parameter values so they cannot drift again.
- Parameter documentation is made generic: cohort-specific figures and a
  site-specific SLURM account example were replaced with neutral wording and a
  round example cohort, so the docs stand on their own.

## v3.5.0

Makes `-resume` usable again without giving up the low storage footprint, exposes
two more knobs, fixes symlinked input under the container, and adds a pack viewer.
No output format changes.

### Added

- **`--publish_mode`** (`link` | `copy` | `move`, default `link`). Controls how
  each Stage 1 pack reaches `output_dir`. `link` hard-links the pack instead of
  moving it, so the pack also stays in the work dir and **`-resume` can skip
  finished accessions**, with **no second copy of the data** (one inode). It
  requires `output_dir` and the work dir on the same filesystem, which is checked
  at launch (a clear error, not a mid-run failure). `copy` works on any layout at
  ~2x storage; `move` is the previous behaviour (lowest footprint, no resume).
- **`--queue_size`** (default 200) — the maximum number of jobs the SLURM
  executor keeps submitted (queued + running) at once.
- **`kbin_dump`** — a standalone tool to inspect a Stage 1 `.kbin` pack or export
  its k-mers as text (`<kmer><TAB><count>`). It reads `k` from the pack's footer,
  so one binary reads a pack of any k. Shipped to `results/bin/` on every run,
  committed as a prebuilt static binary under `tools/bin/`, and attached to each
  tagged release.
- **Prebuilt static binaries** of `kbin_dump` and `bits_to_text` under
  `tools/bin/` and as release assets, so they can be used without running the
  pipeline or compiling.

### Fixed

- **Symlinked FASTQ in `data_dir` now work under the `slurm_container` profile.**
  Singularity binds `data_dir`, but a symlink's target lives outside it and
  dangled inside the container; the real target directories are now bound too.
  (Non-container profiles were unaffected.)

### Changed

- The `mode: 'move'` publish that broke `-resume` is now `--publish_mode` (see
  above), defaulting to `link`.
- Documentation: the benchmarks are rewritten as a self-contained, current-state
  report with a per-release resource comparison across every category.

## v3.4.0

Speeds up Stage 1 input decompression. No format or interface changes; output is
byte-identical to v3.3.0.

### Changed — performance

- **Stage 1 (`KMER_COUNT`) decompresses the two mates concurrently** (one reader
  thread per file, `--kmer_count_read_threads`, default 2). Decompression was the
  single-threaded bottleneck of the read phase once the finalize sort was
  parallelised; reading the mates in parallel roughly halves it on NVMe / parallel
  filesystems. Set `--kmer_count_read_threads 1` to serialise on a single spinning
  disk, where concurrent reads seek-thrash.
- **Gzip inflate uses zlib-ng when available** (~3× faster than stock zlib on
  FASTQ: 23.2 s → 7.4 s per rice mate), selected at compile time and falling back
  to stock zlib so the non-container profiles still build on clusters without it.
  The container image now ships zlib-ng (built from source). On few-core nodes the
  read phase is bounded by the parallel k-mer *encoding* rather than inflate, so
  the two libraries converge there; the zlib-ng win shows on many-core nodes where
  encoding is cheap and the (fixed, two-thread) inflate is the wall.

### Added

- **`--kmer_count_read_threads`** (default `2`).

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
