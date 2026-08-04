# tools

Standalone helper tools that read this pipeline's outputs, usable without running
the pipeline.

Prebuilt, statically linked **x86-64 Linux** binaries (committed here so you can
run them without compiling; the same binaries are attached to each release):

- **`kbin_dump`** — inspect or export a Stage 1 `.kbin` pack: the whole pack to
  stdout, or `--all_bins DIR` to write one file per bin under `DIR/<accession>/`,
  decoded in parallel and gzipped in-process. Reads `k` from the pack's footer,
  so one binary reads a pack of any k.
- **`bits_to_text`** — convert the bit-packed matrix to and from text.

Usage for both is documented in the top-level [README](../README.md). Being
statically linked, the binaries run on an HPC login node or a laptop regardless of
its glibc / libstdc++ version.

## Provenance and rebuilding

The binaries are a convenience snapshot built from `src/` (with `g++ -O2 -static`).
The authoritative source is `src/kbin_dump.cpp` and `src/bits_to_text.cpp`. If you
change those, rebuild with:

```bash
make -C src kbin_dump bits_to_text
```

Every pipeline run also compiles fresh copies into its own `results/bin/`.
