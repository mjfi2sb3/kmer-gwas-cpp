# tools

Standalone helper tools that read this pipeline's outputs, usable without running
the pipeline.

Prebuilt, statically linked **x86-64 Linux** binaries (committed here so you can
run them without compiling; the same binaries are attached to each release):

- **`kbin_dump`** — inspect or export a Stage 1 `.kbin` pack. Reads `k` from the
  pack's footer, so one binary reads a pack of any k.
- **`bits_to_text`** — convert the bit-packed matrix to and from text.

Also here:

- **`bits_to_text.py`** — a portable Python version of `bits_to_text`, with the
  identical interface and output (no build needed; ~14× slower on large matrices).

Usage for all three is documented in the top-level [README](../README.md). Being
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
