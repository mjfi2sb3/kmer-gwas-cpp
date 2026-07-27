#!/usr/bin/env python3
"""Convert a k-mer matrix between the bit-packed and text forms.

    text:  <k-mer><TAB><v0><delim><v1>...<v(N-1)>     (each v is 0 or 1)
    bits:  <k-mer><TAB><hex>                          (N bits, hex-encoded)

where in the bits form <hex> is ceil(N/8) bytes, one bit per accession, LSB-first
within each byte (accession i is bit (i & 7) of byte (i >> 3)).

Runs both ways:
    --decode  (default)  bits -> text   (expand)
    --encode             text -> bits   (pack; a matrix "compressor")

Because the last hex byte may contain padding bits, the number of accessions is
required (via the same accessions file the pipeline was run with, or -n): it
trims padding when decoding and validates the row width when encoding.

Input and output may be plain or gzip-compressed; gzip is chosen for the output
when its name ends in .gz. Both default to stdin/stdout.

For very large matrices the C++ build (src/bits_to_text.cpp, `make bits_to_text`)
is ~14x faster and produces byte-identical output.

Examples
--------
    bits_to_text.py -a accessions.txt  bin0_bits.tsv.gz  bin0_text.tsv.gz
    bits_to_text.py --encode -a accessions.txt  bin0_text.tsv.gz  bin0_bits.tsv.gz
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
    """Open a path, or stdin/stdout when path is None or '-', honouring .gz."""
    if path is None or path == "-":
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
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--decode", dest="encode", action="store_false", default=False,
                      help="bits -> text (default)")
    mode.add_argument("--encode", dest="encode", action="store_true",
                      help="text -> bits")
    p.add_argument("--delimiter", choices=("tab", "space", "none"), default="tab",
                   help="separator of the TEXT side, read or written (default: tab)")
    p.add_argument("--header", action="store_true",
                   help="(decode only) prefix an accession-name header row; needs --accessions")
    p.add_argument("infile", nargs="?", help="input matrix (default: stdin)")
    p.add_argument("outfile", nargs="?", help="output matrix (default: stdout)")
    args = p.parse_args(argv)

    if args.accessions:
        names = [l.strip() for l in open(args.accessions) if l.strip()]
        n_acc = len(names)
    else:
        names = None
        n_acc = args.num_accessions
    if n_acc <= 0:
        p.error("accession count must be positive")
    if args.header and args.encode:
        p.error("--header applies to --decode only")
    if args.header and names is None:
        p.error("--header requires --accessions (to know the names)")

    delim = {"tab": "\t", "space": " ", "none": ""}[args.delimiter]
    n_bytes = (n_acc + 7) // 8

    fin = open_maybe_gz(args.infile, "rt")
    fout = open_maybe_gz(args.outfile, "wt")

    if args.header and not args.encode:
        fout.write("kmer\t" + delim.join(names) + "\n")

    lineno = 0
    for line in fin:
        lineno += 1
        line = line.rstrip("\n")
        if not line:
            continue
        try:
            kmer, rest = line.split("\t", 1)
        except ValueError:
            sys.exit(f"line {lineno}: expected '<kmer><TAB>...', got: {line[:60]!r}")

        if not args.encode:
            # decode: hex -> 0/1 values
            raw = bytes.fromhex(rest)
            if len(raw) != n_bytes:
                sys.exit(f"line {lineno}: row has {len(raw)} bytes but {n_acc} accessions "
                         f"need {n_bytes}. Wrong --accessions/-n, or the input is not "
                         f"bits-encoded.")
            vals = ("1" if (raw[i >> 3] >> (i & 7)) & 1 else "0" for i in range(n_acc))
            fout.write(kmer + "\t" + delim.join(vals) + "\n")
        else:
            # encode: 0/1 values -> hex. A value is "present" unless it is "0",
            # which collapses a count matrix to presence/absence as well.
            fields = list(rest) if delim == "" else rest.split(delim)
            if len(fields) != n_acc:
                sys.exit(f"line {lineno}: row has {len(fields)} values but {n_acc} "
                         f"accessions were expected. Wrong --accessions/-n or --delimiter.")
            raw = bytearray(n_bytes)
            for i, v in enumerate(fields):
                if v != "0":
                    raw[i >> 3] |= 1 << (i & 7)
            fout.write(kmer + "\t" + raw.hex() + "\n")

    if args.outfile and args.outfile != "-":
        fout.close()


if __name__ == "__main__":
    main()
