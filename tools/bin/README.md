# Prebuilt tools

Ready-to-run, statically linked **x86-64 Linux** binaries of the two standalone
helper tools, committed here so you can use them without running the pipeline or
compiling anything:

- **`kbin_dump`** — inspect or export a Stage 1 `.kbin` pack. Reads `k` from the
  pack's footer, so one binary reads a pack of any k.
- **`bits_to_text`** — convert the bit-packed matrix to and from text.

Being statically linked, they run on an HPC login node or a laptop regardless of
its glibc / libstdc++ version. Usage for both is documented in the top-level
[README](../../README.md).

## Provenance and rebuilding

These are a convenience snapshot built from `src/` (with `g++ -O2 -static`). The
authoritative source is `src/kbin_dump.cpp` and `src/bits_to_text.cpp`. If you
change those, rebuild with:

```bash
make -C src kbin_dump bits_to_text     # writes into src/ (built at your machine's k for kmer_count; these two tools are k-agnostic)
```

The matching static binaries are also attached to each release, and every
pipeline run compiles fresh copies into its own `results/bin/`.
