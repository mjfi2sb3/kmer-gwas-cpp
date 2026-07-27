#!/usr/bin/env python3
"""Convert a bits-encoded k-mer matrix (--encoding bits) to a text matrix.

The bits encoding stores each row as

    <k-mer><TAB><hex>

where <hex> is ceil(N/8) bytes, one bit per accession, LSB-first within each
byte (accession i is bit (i & 7) of byte (i >> 3)). It is compact but only
presence/absence, and not readable by tools that expect delimited columns.

This expands it back to the text form

    <k-mer><TAB><v0><delim><v1>...<v(N-1)>

Because the last hex byte may contain padding bits, the number of accessions is
required (via the same accessions file the pipeline was run with, or -n).

Input and output may be plain or gzip-compressed; gzip is chosen for the output
when its name ends in .gz. Both default to stdin/stdout.

Examples
--------
    bits_to_text.py -a accessions.txt  bin0_matrix.tsv.gz  bin0_text.tsv.gz
    zcat bin0_matrix.tsv.gz | bits_to_text.py -n 12600 --delimiter tab > out.tsv
"""

import argparse
import gzip
import signal
import sys

# Exit quietly when a downstream consumer (head, less, ...) closes the pipe,
# instead of printing a BrokenPipeError traceback.
try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass  # SIGPIPE is not available on Windows


def open_maybe_gz(path, mode):
    """Open a path, or stdin/stdout when path is None, honouring a .gz suffix."""
    if path is None:
        stream = sys.stdin if "r" in mode else sys.stdout
        return stream.buffer if "b" in mode else stream
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def count_accessions(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accessions", metavar="FILE",
                     help="accessions file; the number of non-blank lines is the "
                          "accession count (and column order)")
    src.add_argument("-n", "--num-accessions", type=int, metavar="N",
                     help="accession count, if the accessions file is unavailable")
    p.add_argument("--delimiter", choices=("tab", "space", "none"), default="tab",
                   help="value separator in the output (default: tab)")
    p.add_argument("--header", action="store_true",
                   help="prefix an accession-name header row (requires --accessions)")
    p.add_argument("infile", nargs="?", help="bits matrix (default: stdin)")
    p.add_argument("outfile", nargs="?", help="text matrix (default: stdout)")
    args = p.parse_args(argv)

    if args.accessions:
        names = [l.strip() for l in open(args.accessions) if l.strip()]
        n_acc = len(names)
    else:
        names = None
        n_acc = args.num_accessions
    if n_acc <= 0:
        p.error("accession count must be positive")
    if args.header and names is None:
        p.error("--header requires --accessions (to know the names)")

    delim = {"tab": "\t", "space": " ", "none": ""}[args.delimiter]
    n_bytes = (n_acc + 7) // 8

    fin = open_maybe_gz(args.infile, "rt")
    fout = open_maybe_gz(args.outfile, "wt")

    if args.header:
        fout.write("kmer\t" + delim.join(names) + "\n")

    lineno = 0
    for line in fin:
        lineno += 1
        line = line.rstrip("\n")
        if not line:
            continue
        try:
            kmer, hexstr = line.split("\t", 1)
        except ValueError:
            sys.exit(f"line {lineno}: expected '<kmer><TAB><hex>', got: {line[:60]!r}")
        raw = bytes.fromhex(hexstr)
        if len(raw) != n_bytes:
            sys.exit(f"line {lineno}: row has {len(raw)} bytes but {n_acc} accessions "
                     f"need {n_bytes}. Wrong --accessions/-n, or the file is not "
                     f"bits-encoded.")
        vals = ("1" if (raw[i >> 3] >> (i & 7)) & 1 else "0" for i in range(n_acc))
        fout.write(kmer + "\t" + delim.join(vals) + "\n")

    if args.outfile:
        fout.close()


if __name__ == "__main__":
    main()
