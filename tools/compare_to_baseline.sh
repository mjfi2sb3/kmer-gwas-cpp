#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Compare the current pipeline against a baseline git revision, on real data.
#
# Builds BOTH versions from source, runs each over the same accessions, and
# diffs the resulting matrices. Use this before trusting the rewritten pipeline
# on a production cohort.
#
#   tools/compare_to_baseline.sh --data /path/to/fastq --accessions accs.txt
#   tools/compare_to_baseline.sh --data DIR --accessions F --bins 8 --k 51 \
#                                --baseline 3297c96 --workdir /scratch/cmp
#
# TWO DIFFERENCES ARE EXPECTED AND INTENTIONAL — the script accounts for both:
#
#   1. The baseline wrote a spurious second tab after the k-mer column
#      ("KMER\t\t0\t1"). That was a bug; it is normalised away before diffing.
#   2. The baseline silently dropped the poly-A k-mer (AAA...A), because A
#      encodes as 00 and the all-zero key doubled as an "invalid" marker. The
#      current version keeps it, so the new output may have ONE extra row per
#      bin where a poly-A tract exists. Reported separately, not as a failure.
#
# Row ORDER differs by design (the baseline iterated a hash map; the current
# version emits sorted k-mers), so rows are compared as sets.
# ---------------------------------------------------------------------------
set -euo pipefail

DATA=""; ACCS=""; BINS=8; K=51; BASELINE="3297c96"; WORKDIR=""
MINCOUNT=2; THRESHOLD=0; COUNT=n; CORE=n

usage() { sed -n '2,25p' "$0"; exit 1; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data)       DATA="$2"; shift 2 ;;
        --accessions) ACCS="$2"; shift 2 ;;
        --bins)       BINS="$2"; shift 2 ;;
        --k)          K="$2"; shift 2 ;;
        --baseline)   BASELINE="$2"; shift 2 ;;
        --workdir)    WORKDIR="$2"; shift 2 ;;
        --threshold)  THRESHOLD="$2"; shift 2 ;;
        --count)      COUNT="$2"; shift 2 ;;
        --core)       CORE="$2"; shift 2 ;;
        -h|--help)    usage ;;
        *) echo "unknown option: $1" >&2; usage ;;
    esac
done
[[ -n "$DATA" && -n "$ACCS" ]] || usage

REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORKDIR="${WORKDIR:-$(mktemp -d)}"
mkdir -p "$WORKDIR"
NACC=$(grep -cve '^\s*$' "$ACCS")

echo "=========================================================="
echo " baseline   : $BASELINE"
echo " current    : $(git -C "$REPO" rev-parse --short HEAD) + working tree"
echo " accessions : $NACC   bins: $BINS   k: $K"
echo " workdir    : $WORKDIR"
echo "=========================================================="

# -- build the baseline out of git, without touching the working tree --------
echo "[1/5] building baseline binaries from $BASELINE ..."
OLD="$WORKDIR/baseline"; mkdir -p "$OLD/src"
for f in kmer_count_v3.cpp matrix_merge.cpp mmap_io.cpp mmap_io.hpp thread_pool.hpp; do
    git -C "$REPO" show "$BASELINE:src/$f" > "$OLD/src/$f" 2>/dev/null || true
done
( cd "$OLD/src"
  g++ -std=c++17 -O3 -march=native -pthread -o "$OLD/kmer_count"   kmer_count_v3.cpp mmap_io.cpp -lz
  g++ -std=c++17 -O3 -march=native -pthread -o "$OLD/matrix_merge" matrix_merge.cpp  mmap_io.cpp )
# The baseline hardcodes k=51; refuse to compare at any other k.
if [[ "$K" != "51" ]]; then
    echo "  NOTE: the baseline only supports k=51. Re-run with --k 51 to compare," >&2
    echo "        or skip the comparison and validate the new pipeline alone." >&2
    exit 2
fi

echo "[2/5] building current binaries (k=$K) ..."
NEW="$WORKDIR/current"; mkdir -p "$NEW"
( cd "$REPO/src" && make -s KMER_K="$K" BIN="$NEW" >/dev/null )

# -- stage 1 -----------------------------------------------------------------
echo "[3/5] running Stage 1 (both versions) ..."
mkdir -p "$WORKDIR/kc_old" "$WORKDIR/kc_new"
while read -r acc; do
    [[ -z "$acc" ]] && continue
    R1=""; R2=""
    for e in _1.fq _1.fastq _1.fq.gz _1.fastq.gz; do [[ -f "$DATA/$acc$e" ]] && { R1="$DATA/$acc$e"; break; }; done
    for e in _2.fq _2.fastq _2.fq.gz _2.fastq.gz; do [[ -f "$DATA/$acc$e" ]] && { R2="$DATA/$acc$e"; break; }; done
    [[ -n "$R1" && -n "$R2" ]] || { echo "  ERROR: no FASTQ pair for $acc in $DATA" >&2; exit 1; }
    ( cd "$WORKDIR/kc_old" && "$OLD/kmer_count" "$acc" "$BINS" ./ "$R1" "$R2" >/dev/null )
    ( cd "$WORKDIR/kc_new" && "$NEW/kmer_count" "$acc" "$BINS" ./ "$R1" "$R2" 32 "$MINCOUNT" >/dev/null )
    echo "    $acc done"
done < "$ACCS"

# -- stage 2 -----------------------------------------------------------------
echo "[4/5] running Stage 2 (both versions) ..."
mkdir -p "$WORKDIR/mx_old" "$WORKDIR/mx_new"
for ((b=0; b<BINS; b++)); do
    ( cd "$WORKDIR/mx_old" && "$OLD/matrix_merge" --input "$WORKDIR/kc_old/" --accessions "$ACCS" \
        --index "$b" --threshold "$THRESHOLD" --delimiter tab --count "$COUNT" --core "$CORE" \
        --bins "$BINS" --threads 1 >/dev/null )
    ( cd "$WORKDIR/mx_new" && "$NEW/matrix_merge" --input "$WORKDIR/kc_new/" --accessions "$ACCS" \
        --index "$b" --threshold "$THRESHOLD" --delimiter tab --count "$COUNT" --core "$CORE" \
        --bins "$BINS" --threads 1 >/dev/null )
done

# -- compare -----------------------------------------------------------------
echo "[5/5] comparing ..."
cat "$WORKDIR"/mx_old/matrix_*/*_matrix.tsv | sed 's/\t\t/\t/' | sort > "$WORKDIR/old.rows"
cat "$WORKDIR"/mx_new/matrix_*/*_matrix.tsv                    | sort > "$WORKDIR/new.rows"

POLYA=$(printf 'A%.0s' $(seq 1 "$K"))
only_old=$(comm -23 "$WORKDIR/old.rows" "$WORKDIR/new.rows" | wc -l)
only_new=$(comm -13 "$WORKDIR/old.rows" "$WORKDIR/new.rows" | wc -l)
only_new_nonpolya=$(comm -13 "$WORKDIR/old.rows" "$WORKDIR/new.rows" | grep -vc "^$POLYA" || true)
polya_rows=$(comm -13 "$WORKDIR/old.rows" "$WORKDIR/new.rows" | grep -c "^$POLYA" || true)

echo
echo "  rows: baseline $(wc -l < "$WORKDIR/old.rows"), current $(wc -l < "$WORKDIR/new.rows")"
echo "  rows only in baseline            : $only_old"
echo "  rows only in current             : $only_new"
echo "    of which poly-A (expected, #6) : $polya_rows"
echo "    unexplained                    : $only_new_nonpolya"

if [[ "$CORE" == "y" ]]; then
    cat "$WORKDIR"/mx_old/matrix_*/*_core.txt 2>/dev/null | sort > "$WORKDIR/old.core" || true
    cat "$WORKDIR"/mx_new/matrix_*/*_core.txt 2>/dev/null | sort > "$WORKDIR/new.core" || true
    if diff -q "$WORKDIR/old.core" "$WORKDIR/new.core" >/dev/null 2>&1;
       then echo "  core k-mer files                 : identical"
       else echo "  core k-mer files                 : DIFFER"; fi
fi

echo
if [[ "$only_old" -eq 0 && "$only_new_nonpolya" -eq 0 ]]; then
    echo "  *** PASS — output matches the baseline (poly-A difference is intentional) ***"
    echo "  workdir kept at $WORKDIR"
    exit 0
else
    echo "  *** FAIL — unexplained differences; workdir kept at $WORKDIR ***"
    exit 1
fi
