#!/usr/bin/env bash
# ===========================================================================
#  check_kmc_db.sh — test kmc_to_pack against a KMC database YOU built.
# ===========================================================================
#
#  Unlike test_kmc_to_pack.sh (which generates its own reads and checks
#  byte-parity against the fastq path), this script takes a KMC database you
#  supply and verifies the conversion using KMC itself as ground truth. No reads
#  are needed.
#
#  Default (cheap, no dump): convert the database, confirm the converter read
#  EVERY k-mer the database holds (in == distinct_kmers), and confirm the pack is
#  structurally valid. Safe on huge production databases.
#
#  --full (exact, expensive): also `kmc_tools dump` the whole database and assert
#  the pack decodes to exactly the same k-mer/count set (counts clamped to 65535,
#  and filtered at min_count when min_count is above the database's -ci). Only for a CANONICAL database, and only when you
#  can afford to dump it — the dump is hundreds of GB for a large genome.
#
#  Requires on PATH: g++ (always); kmc_tools (only for --full), e.g.
#  `module load kmc` on clusters that use environment modules.
#
#  The converted pack is KEPT by default (named <db_basename>.kbin) in the
#  directory you launched from; use --keep <dir> to choose where, or --no-keep to
#  discard it (verify only).
#
#  Usage:  tools/check_kmc_db.sh <kmc_db_prefix> <K> [num_bins] [min_count] [options]
#            <kmc_db_prefix>  path to the database WITHOUT .kmc_pre/.kmc_suf
#            <K>              odd, 15..63, MUST equal the database's k
#            num_bins         default 1500
#            min_count        default 0: keep the database as it is (its own -ci
#                             applies); a value above its -ci is enforced as well
#            --keep <dir>     write <db_basename>.kbin into <dir> (default: launch dir)
#            --no-keep        do not keep the pack (verify only)
#            --full           exact dump comparison (canonical DBs only; expensive)
# ===========================================================================
set -euo pipefail
unset CDPATH   # a CDPATH setting makes `cd` print the path into captured output

DB=""; K=""; NUM_BINS=1500; MINC=0; FULL=0; KEEP=1; KEEPDIR=""
ORIG_PWD="$PWD"                      # where the user launched us (before any cd)
positional=()
while [ $# -gt 0 ]; do
    case "$1" in
        --full)     FULL=1; shift ;;
        --no-keep)  KEEP=0; shift ;;
        --keep)     KEEP=1; KEEPDIR="${2:-}"; [ -n "$KEEPDIR" ] || { echo "ERROR: --keep needs a <dir>" >&2; exit 2; }; shift 2 ;;
        --keep=*)   KEEP=1; KEEPDIR="${1#--keep=}"; shift ;;
        -h|--help)  grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 0 ;;
        -*)         echo "ERROR: unknown option '$1'" >&2; exit 2 ;;
        *)          positional+=("$1"); shift ;;
    esac
done
DB="${positional[0]:-}"; K="${positional[1]:-}"
NUM_BINS="${positional[2]:-1500}"; MINC="${positional[3]:-0}"
[ -n "$DB" ] && [ -n "$K" ] || { grep -E '^#( |$)' "$0" | sed -E 's/^# ?//'; exit 2; }
[ -f "$DB.kmc_pre" ] && [ -f "$DB.kmc_suf" ] || { echo "ERROR: no KMC database at '$DB' (need $DB.kmc_pre and $DB.kmc_suf)" >&2; exit 2; }

# Resolve this script's own location (symlink-safe), independent of CWD.
self="${BASH_SOURCE[0]:-$0}"
if command -v readlink >/dev/null 2>&1; then rp="$(readlink -f "$self" 2>/dev/null || true)"; [ -n "$rp" ] && self="$rp"; fi
here="$(cd "$(dirname "$self")" && pwd)"
SRC="${KMER_GWAS_SRC:-$here/../src}"
[ -f "$SRC/kmc_to_pack.cpp" ] || { echo "ERROR: cannot find pipeline src/ (looked in '$SRC'); set KMER_GWAS_SRC=/path/to/kmer-gwas-cpp/src" >&2; exit 2; }
SRC="$(cd "$SRC" && pwd)"

command -v g++ >/dev/null 2>&1 || { echo "ERROR: g++ not on PATH" >&2; exit 2; }
DBABS="$(cd "$(dirname "$DB")" && pwd)/$(basename "$DB")"

# Resolve the keep directory (default: the launch directory), while PWD is still
# the launch dir. The kept pack survives the temp-dir cleanup.
if [ "$KEEP" = 1 ]; then
    [ -n "$KEEPDIR" ] || KEEPDIR="$ORIG_PWD"
    case "$KEEPDIR" in /*) ;; *) KEEPDIR="$ORIG_PWD/$KEEPDIR" ;; esac
fi

work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT; cd "$work"
CXXFLAGS="-std=c++17 -O3 -march=native -pthread"
API="$SRC/kmc_api/kmc_file.cpp $SRC/kmc_api/kmer_api.cpp $SRC/kmc_api/mmer.cpp"

echo "check_kmc_db: db=$DBABS  K=$K  num_bins=$NUM_BINS  min_count=$MINC  full=$FULL"
g++ $CXXFLAGS -DKMER_K="$K" -o kmc_to_pack "$SRC/kmc_to_pack.cpp" $API
g++ $CXXFLAGS -DKMER_K="$K" -o kbin_dump   "$SRC/kbin_dump.cpp" -lz

# -- convert (converter aborts itself on a k mismatch) ----------------------
./kmc_to_pack "$DBABS" "$NUM_BINS" ./ ACC "$MINC" 0 2>conv.log || { echo "FAIL: converter aborted"; cat conv.log; exit 1; }
cat conv.log

# Keep the produced pack (on by default), named after the database, before the
# temp dir is cleaned up. Done right after conversion so the pack is saved even
# if a later check fails.
if [ "$KEEP" = 1 ]; then
    mkdir -p "$KEEPDIR"
    kept="$KEEPDIR/$(basename "$DBABS").kbin"
    cp -f ACC.kbin "$kept"
    echo "kept pack: $kept"
fi
val() { grep -oP "$1=\\K[A-Za-z0-9]+" conv.log | head -1; }
TOTAL=$(val distinct_kmers); BOTH=$(val both_strands); CI=$(grep -oP 'min_count\(-ci\)=\K[0-9]+' conv.log | head -1)
# The filter the converter applied: min_count only when it is above the database -ci.
if [ "$MINC" -gt "${CI:-1}" ]; then EFF=$MINC; else EFF=1; fi
IN=$(grep -oP 'in=\K[0-9]+' conv.log); OUT=$(grep -oP 'out=\K[0-9]+' conv.log)

fail=0
# 1) the converter must have read every k-mer the database reports
if [ "$IN" = "$TOTAL" ]; then echo "PASS: read all k-mers (in=$IN == distinct_kmers=$TOTAL)";
else echo "FAIL: in=$IN != distinct_kmers=$TOTAL"; fail=1; fi
# 2) output is a sane subset
if [ "${OUT:-0}" -gt 0 ] && [ "${OUT:-0}" -le "${IN:-0}" ]; then echo "PASS: out=$OUT records (<= in, >0)";
else echo "FAIL: out=$OUT not in (0, in]"; fail=1; fi
# 3) the pack is structurally valid and self-describes the right k/bins
info="$(./kbin_dump ACC.kbin --info)" || { echo "FAIL: pack did not open"; exit 1; }
echo "$info" | sed 's/^/    /'
echo "$info" | grep -q "k            : $K" && echo "PASS: pack k == $K" || { echo "FAIL: pack k mismatch"; fail=1; }
tot_pack=$(echo "$info" | grep -oP 'total k-mers : \K[0-9]+')
[ "$tot_pack" = "$OUT" ] && echo "PASS: pack holds out records ($tot_pack)" || { echo "FAIL: pack total $tot_pack != out $OUT"; fail=1; }

# -- optional exact comparison against KMC's own dump -----------------------
if [ "$FULL" = 1 ]; then
    command -v kmc_tools >/dev/null 2>&1 || { echo "FAIL: --full needs kmc_tools on PATH"; exit 1; }
    if [ "$BOTH" != "true" ]; then
        echo "SKIP full compare: database is non-canonical (both_strands=$BOTH); exact check needs re-canonicalisation of the dump, which this script does not do. Structural checks above still hold."
    else
        echo "full compare: dumping the database (this can be very large)..."
        kmc_tools transform "$DBABS" dump /dev/stdout 2>/dev/null \
          | awk -v m="$EFF" '{c=$2+0; if(c>65535)c=65535; if(c>=m) print $1"\t"c}' | LC_ALL=C sort > exp.txt
        ./kbin_dump ACC.kbin | LC_ALL=C sort > got.txt
        if diff -q exp.txt got.txt >/dev/null; then
            echo "PASS: pack decodes to exactly the KMC dump set (counts clamped@65535, filtered@$EFF)"
        else
            echo "FAIL: pack differs from the KMC dump set; first differences:"; diff exp.txt got.txt | head; fail=1
        fi
    fi
fi

[ "$fail" = 0 ] && echo "ALL CHECKS PASSED" || { echo "SOME CHECKS FAILED"; exit 1; }
