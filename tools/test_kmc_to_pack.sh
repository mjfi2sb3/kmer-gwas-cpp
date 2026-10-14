#!/usr/bin/env bash
# ===========================================================================
#  test_kmc_to_pack.sh — reproducible correctness test for kmc_to_pack.
# ===========================================================================
#
#  Proves that converting a KMC database with kmc_to_pack yields a pack whose
#  data (bins and offset table) is BYTE-IDENTICAL to the one kmer_count produces
#  from the same reads, when both are filtered at the same minimum count. The
#  packs' metadata differs by design (source, program), so it is checked
#  separately.
#
#  It builds everything from the repo in a temp dir, generates synthetic reads
#  (including k-mers seen only once), builds KMC databases with the real `kmc`,
#  and checks:
#    - canonical, -b and KMC1-format databases at -ci1, kept as they are
#    - a -ci2 database kept as it is, and a -ci1 database with min_count 2 enforced
#    - -cs1 and `kmc_tools compact` databases (every count 1): same k-mers
#    - an enforced min_count that the counts never reach stops the job
#    - every k-mer is canonical; the metadata records input_mode, min_kmer_count
#      and its source, and the count limit of the database counter
#    - k mismatch, non-KMC file and missing .kmc_suf all stop the job
#
#  Requires on PATH: g++, python3, and KMC (`kmc`, `kmc_tools`), e.g. `module load kmc` on
#  clusters that use environment modules.
#
#  Usage:  tools/test_kmc_to_pack.sh [K] [NUM_BINS]
#          K        odd, 15..63   (default 51)
#          NUM_BINS                (default 16)
# ===========================================================================
set -euo pipefail
unset CDPATH   # a CDPATH setting makes `cd` print the path into captured output

K="${1:-51}"
NUM_BINS="${2:-16}"

# Resolve this script's own real location (symlink-safe), independent of the
# current working directory, so it can be launched from anywhere.
self="${BASH_SOURCE[0]:-$0}"
if command -v readlink >/dev/null 2>&1; then
    rp="$(readlink -f "$self" 2>/dev/null || true)"; [ -n "$rp" ] && self="$rp"
fi
here="$(cd "$(dirname "$self")" && pwd)"

# Pipeline sources: <repo>/src relative to this script by default. Override with
# KMER_GWAS_SRC when running a copy of this script from outside the repo.
SRC="${KMER_GWAS_SRC:-$here/../src}"
if [ ! -f "$SRC/kmc_to_pack.cpp" ]; then
    echo "ERROR: cannot find the pipeline src/ (looked in '$SRC')." >&2
    echo "  Run this script from within the repo, or set" >&2
    echo "  KMER_GWAS_SRC=/path/to/kmer-gwas-cpp/src" >&2
    exit 2
fi
SRC="$(cd "$SRC" && pwd)"
CXXFLAGS="-std=c++17 -O3 -march=native -pthread"
API="$SRC/kmc_api/kmc_file.cpp $SRC/kmc_api/kmer_api.cpp $SRC/kmc_api/mmer.cpp"

need() { command -v "$1" >/dev/null 2>&1 || { echo "MISSING: $1 not on PATH (e.g. 'module load kmc' on clusters that use environment modules)" >&2; exit 2; }; }
need g++; need python3; need kmc; need kmc_tools

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work"
echo "test_kmc_to_pack: K=$K NUM_BINS=$NUM_BINS  workdir=$work"

# -- build the binaries at this K --------------------------------------------
g++ $CXXFLAGS -DKMER_K="$K" -o kmer_count  "$SRC/kmer_count_v3.cpp" -lz
g++ $CXXFLAGS -DKMER_K="$K" -o kmc_to_pack "$SRC/kmc_to_pack.cpp" $API
g++ $CXXFLAGS -o kbin_dump "$SRC/kbin_dump.cpp" -lz

# -- synthetic reads ----------------------------------------------------------
# R1: 400 reads from a random genome, half reverse-complemented (k-mers seen many
# times). R2: 20 reads from a second genome (k-mers seen about once).
python3 - "$K" <<'PY'
import random
random.seed(1)
comp = {'A':'T','C':'G','G':'C','T':'A'}
rc = lambda s: ''.join(comp[c] for c in reversed(s))
def reads(path, genome_len, n, prefix):
    G = ''.join(random.choice('ACGT') for _ in range(genome_len))
    with open(path, 'w') as fh:
        for i in range(n):
            p = random.randint(0, len(G) - 150); s = G[p:p+150]
            if i % 2: s = rc(s)
            fh.write("@%s%d\n%s\n+\n%s\n" % (prefix, i, s, 'I' * 150))
reads('r1.fastq', 3000, 400, 'a')
reads('r2.fastq', 20000, 20, 'b')
reads('hi.fastq', 600, 3000, 'h')   # k-mer counts far above 255
PY
cat r1.fastq r2.fastq > reads.fastq

# -- fastq path (reference), at min_count 1 and 2 ------------------------------
mkdir -p fq1 fq2
( cd fq1 && ../kmer_count ACC "$NUM_BINS" ./ ../r1.fastq ../r2.fastq 1 1 2 ) >/dev/null
( cd fq2 && ../kmer_count ACC "$NUM_BINS" ./ ../r1.fastq ../r2.fastq 1 2 2 ) >/dev/null

# -- KMC databases --------------------------------------------------------------
# KMC's last argument is a working directory that MUST already exist.
mkdir -p tmp
run_kmc() { "$@" >kmc.log 2>&1 || { echo "FAIL: kmc command failed: $*"; cat kmc.log; exit 1; }; }
run_kmc kmc -k"$K" -ci1 -cs65535    reads.fastq db_ci1   tmp
run_kmc kmc -k"$K" -b -ci1 -cs65535 reads.fastq db_b     tmp
run_kmc kmc -k"$K" -ci2 -cs65535    reads.fastq db_ci2   tmp
run_kmc kmc -k"$K" -ci2 -cs1        reads.fastq db_cs1   tmp
run_kmc kmc -k"$K" -ci1 -cs3        reads.fastq db_cs3   tmp
run_kmc kmc_tools transform db_ci1 sort db_kmc1
run_kmc kmc_tools transform db_ci2 compact db_compact
run_kmc kmc -k"$K" -ci1 -cs255   hi.fastq db_hi255 tmp
run_kmc kmc -k"$K" -ci1 -cs65535 hi.fastq db_hi     tmp

# convert <db> <out_dir> [min_count]
convert() { ./kmc_to_pack "$1" "$NUM_BINS" "./$2" ACC "${3:-0}" 0 >/dev/null 2>"$2.log"; }

# Number of bytes of a pack that hold data: the bins plus the offset table.
data_len() {
  python3 - "$1" <<'PY'
import struct, sys
b = open(sys.argv[1], 'rb').read()
index_offset, num_bins = struct.unpack_from('<QI', b, len(b) - 24)
print(index_offset + (num_bins + 1) * 8)
PY
}
# same_data <pack_a> <pack_b>: the data parts are byte-identical
same_data() {
  local la lb; la=$(data_len "$1"); lb=$(data_len "$2")
  [ "$la" = "$lb" ] && cmp -s -n "$la" "$1" "$2"
}
# kmers <pack>: the sorted k-mers of a pack, without counts
kmers() { ./kbin_dump "$1" | cut -f1 | LC_ALL=C sort; }
# counts <pack>: the sorted distinct counts of a pack
counts() { ./kbin_dump "$1" | cut -f2 | sort -u | tr '\n' ' '; }
meta() { ./kbin_dump "$1" --info | sed -n "s/^  $2=//p"; }

fail=0
pass() { echo "PASS: $*"; }
bad()  { echo "FAIL: $*"; fail=1; }

# -- databases kept as they are --------------------------------------------------
for db in ci1 b kmc1; do
  convert "db_$db" "$db"
  same_data fq1/ACC.kbin "$db/ACC.kbin" \
    && pass "$db database, kept as it is == reads at min_count 1 (data byte-identical)" \
    || { bad "$db database differs from reads at min_count 1"; cat "$db.log"; }
done
convert db_ci2 ci2
same_data fq2/ACC.kbin ci2/ACC.kbin \
  && pass "-ci2 database, kept as it is == reads at min_count 2 (data byte-identical)" \
  || bad "-ci2 database differs from reads at min_count 2"

# -- enforced min_count ------------------------------------------------------------
convert db_ci1 ci1_m2 2
same_data fq2/ACC.kbin ci1_m2/ACC.kbin \
  && pass "-ci1 database with min_count 2 enforced == reads at min_count 2" \
  || bad "-ci1 database with min_count 2 differs from reads at min_count 2"
convert db_ci2 ci2_m2 2
grep -q 'filter: skipped' ci2_m2.log && same_data fq2/ACC.kbin ci2_m2/ACC.kbin \
  && pass "-ci2 database with min_count 2: filter skipped, same data" \
  || bad "-ci2 database with min_count 2 was not handled as expected"

# -- databases without real counts (every count 1) ----------------------------------
# kmc -cs1 stores no counter from KMC 3.2 on, and a 1-byte counter before; both
# report count 1. kmc_tools compact also sets every count to 1.
for db in cs1 compact; do
  convert "db_$db" "$db"
  if [ "$(kmers "$db/ACC.kbin")" = "$(kmers fq2/ACC.kbin)" ] && [ "$(counts "$db/ACC.kbin")" = "1 " ]; then
    pass "$db database: same k-mers as reads at min_count 2, every count 1"
  else
    bad "$db database: k-mers or counts not as expected"; cat "$db.log"
  fi
done

# -- an enforced filter the counts never reach must stop the job --------------------
# expect_abort <label> <message-substring> <db_prefix> [min_count]
expect_abort() {
  if ./kmc_to_pack "$3" "$NUM_BINS" ./bad ACC "${4:-0}" 0 >/dev/null 2>bad.log; then
    bad "$1 did NOT abort"
  elif grep -q "$2" bad.log; then
    pass "$1 aborts ($(sed -n 's/^kmc_to_pack: ERROR: //p' bad.log | head -1 | cut -c1-90)...)"
  else
    bad "$1 aborted without the expected message '$2'"; cat bad.log
  fi
}
expect_abort "-cs1 database with min_count 3" "never reach" db_cs1 3
expect_abort "-cs3 database with min_count 5" "never reach" db_cs3 5
convert db_cs3 cs3_m2 2
[ "$(kmers cs3_m2/ACC.kbin)" = "$(kmers fq2/ACC.kbin)" ] \
  && pass "-cs3 database with min_count 2: same k-mers as reads at min_count 2" \
  || bad "-cs3 database with min_count 2: k-mers differ"

# -- invalid input must abort with a clear message ----------------------------------
KM1=$(( K - 2 ))   # still odd
run_kmc kmc -k"$KM1" -ci1 -cs65535 reads.fastq db_wrongk tmp
expect_abort "k mismatch (db k=$KM1, binary k=$K)" "rebuild with -DKMER_K" db_wrongk
printf 'not a kmc database\n' > db_junk.kmc_pre; cp db_junk.kmc_pre db_junk.kmc_suf
expect_abort "non-KMC file" "is not a KMC database file" db_junk
cp db_ci1.kmc_pre db_nosuf.kmc_pre
expect_abort "missing .kmc_suf" "missing file" db_nosuf

# -- metadata --------------------------------------------------------------------------
if [ "$(meta ci1/ACC.kbin input_mode)" = kmc ] && [ "$(meta fq1/ACC.kbin input_mode)" = reads ] \
   && [ "$(meta ci2/ACC.kbin min_kmer_count)" = 2 ] && [ "$(meta ci1_m2/ACC.kbin min_kmer_count)" = 2 ] \
   && [ "$(meta ci1/ACC.kbin min_kmer_count)" = 1 ] && [ "$(meta ci1/ACC.kbin accession)" = ACC ] \
   && [ "$(meta ci1/ACC.kbin canonical_kmers)" = yes ] && [ "$(meta ci1/ACC.kbin count_limit)" = 65535 ] \
   && [ "$(meta ci2/ACC.kbin distinct_kmers_after_filter)" = "$(./kbin_dump ci2/ACC.kbin | wc -l)" ]; then
  pass "metadata records input_mode, accession, min_kmer_count, canonical_kmers, count_limit and distinct k-mers"
else
  bad "metadata not as expected"; ./kbin_dump ci2/ACC.kbin --info
fi

# -- source of the count filter ----------------------------------------------------
( mkdir -p fq_def && cd fq_def && ../kmer_count ACC "$NUM_BINS" ./ ../r1.fastq ../r2.fastq 1 2 2 default ) >/dev/null
if [ "$(meta fq1/ACC.kbin min_kmer_count_source)" = given ] && [ "$(meta fq_def/ACC.kbin min_kmer_count_source)" = default ] \
   && [ "$(meta ci1/ACC.kbin min_kmer_count_source)" = database ] && [ "$(meta ci1_m2/ACC.kbin min_kmer_count_source)" = enforced ] \
   && [ "$(meta ci2_m2/ACC.kbin min_kmer_count_source)" = database ]; then
  pass "min_kmer_count_source: given / default (reads), database / enforced (kmc)"
else
  bad "min_kmer_count_source not as expected"
fi

# -- count limit of the database counter ----------------------------------------------
# -cs255 uses a 1-byte counter: counts stop at 255, and the pack must say so.
convert db_hi255 hi255; convert db_hi hi
want_at=$(./kbin_dump hi/ACC.kbin | awk '$2 >= 255' | wc -l)
got_at=$(meta hi255/ACC.kbin kmers_at_count_limit)
if [ "$(meta hi255/ACC.kbin count_limit)" = 255 ] && [ "$want_at" -gt 0 ] && [ "$got_at" = "$want_at" ] \
   && [ "$(meta hi255/ACC.kbin max_count_per_kmer)" = 255 ] && [ "$(meta hi/ACC.kbin kmers_at_count_limit)" = 0 ] \
   && [ "$(meta hi255/ACC.kbin estimated_kmer_coverage)" != 255 ]; then
  pass "-cs255 database: count_limit 255, kmers_at_count_limit $got_at (as in the -cs65535 database), 255 not taken as coverage"
else
  bad "-cs255 database: count_limit $(meta hi255/ACC.kbin count_limit), kmers_at_count_limit $got_at (want $want_at)"
fi
# -cs1 and compact databases: no counter (limit 1) or a 1-byte counter (limit
# 255), depending on the KMC version; every count is 1, so there is no coverage peak.
for db in cs1 compact; do
  cs=$(grep -oP 'counter_size=\K[0-9]+' "$db.log" | head -1)
  want_limit=$([ "$cs" = 0 ] && echo 1 || echo 255)
  [ "$(meta "$db/ACC.kbin" count_limit)" = "$want_limit" ] && [ "$(meta "$db/ACC.kbin" estimated_kmer_coverage)" = NA ] \
    && pass "$db database (counter_size $cs): count_limit $want_limit, coverage NA" \
    || bad "$db database (counter_size $cs): count_limit $(meta "$db/ACC.kbin" count_limit) (want $want_limit), coverage $(meta "$db/ACC.kbin" estimated_kmer_coverage)"
done

# -- every k-mer in every pack is canonical ---------------------------------------------
for d in fq1 ci1 b kmc1 ci2; do
  ./kbin_dump "$d/ACC.kbin" --all_bins "dump_$d" --no-gzip --threads 1 >/dev/null
done
if python3 - dump_fq1 dump_ci1 dump_b dump_kmc1 dump_ci2 <<'PY'
import glob, sys
comp = str.maketrans('ACGT', 'TGCA')
n = bad = 0
for d in sys.argv[1:]:
    for path in glob.glob(d + '/*/*.tsv'):
        for line in open(path):
            kmer = line.split('\t', 1)[0].strip()
            if not kmer:
                continue
            n += 1
            if kmer > kmer.translate(comp)[::-1]:
                bad += 1
sys.exit(0 if n > 0 and bad == 0 else 1)
PY
then
  pass "all k-mers are canonical (smaller of k-mer and reverse complement)"
else
  bad "non-canonical k-mers found, or no k-mers read"
fi

[ "$fail" = 0 ] && echo "ALL TESTS PASSED (K=$K)" || { echo "SOME TESTS FAILED (K=$K)"; exit 1; }
