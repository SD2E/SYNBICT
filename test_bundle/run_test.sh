#!/usr/bin/env bash
# Run the two features under test:
#   feature 1  features_to_circuits.py -gn   : annotated SBOL -> logic-gate netlist
#   feature 2  circuit_to_truth_table.py     : netlist -> truth table (via Yosys)
# Annotation (FASTA -> annotated SBOL) runs first so the test starts from raw sequence.
#
# Usage:
#   ./run_test.sh                 quick: 6 Cello circuits          (~30 s)
#   ./run_test.sh --cello         all 66 Cello circuits            (~5 min)
#   ./run_test.sh --md5           all 70 MD5 plasmids              (~10 min)
#   ./run_test.sh --all           both
#   ./run_test.sh --quick out_dir  (any mode takes an optional output dir)
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${SYNBICT_REPO:-$(cd "$HERE/.." && pwd)}"
PY="${SYNBICT_PY:-python}"
YOSYS="${YOSYS:-$(command -v yosys)}"
NS="http://examples.org"

MODE="${1:---quick}"
OUT="${2:-$HERE/out}"
QUICK_SET="0x01 0x07 0x70 0xC9 0xE8 0xEA"

for exe in blastn makeblastdb; do
    command -v "$exe" >/dev/null || { echo "MISSING: $exe not on PATH"; exit 1; }
done
[ -x "$YOSYS" ] || { echo "MISSING: yosys (set \$YOSYS)"; exit 1; }
"$PY" -c "import sbol2, pysam, pandas" 2>/dev/null || { echo "MISSING: python deps (sbol2/pysam/pandas) -- wrong env? set \$SYNBICT_PY"; exit 1; }
export PYTHONPATH="$REPO"

# ---------------------------------------------------------------- one dataset
# $1 = dataset dir (cello|md5), $2 = library, $3 = min target length,
# $4 = extra annotation flags, $5 = output subdir, $6.. = circuit names
#      (no names = every FASTA in the dataset)
run_dataset() {
    local ds="$1" lib="$2" minlen="$3" annoflags="$4" out="$5"; shift 5
    local only=("$@")
    mkdir -p "$out"
    # the annotation step writes its BLAST index as ./test.* -- keep it out of the repo
    ( cd "$out" || exit 1
      echo "== build feature index from $(basename "$lib") =="
      "$PY" -m sequences_to_features -n "$NS" -f "$lib" -bi > build_index.log 2>&1 \
          || { echo "index build FAILED, see $out/build_index.log"; exit 1; }

      local ok=0 fail=0
      for F in "$HERE/$ds"/sequences/*.fasta; do
          local C; C=$(basename "$F" .fasta)
          if [ ${#only[@]} -gt 0 ]; then
              [[ " ${only[*]} " == *" $C "* ]] || continue
          fi
          (
              set -e
              # NOTE: no -nms on purpose; see TESTING.md "Known traps"
              # shellcheck disable=SC2086
              "$PY" -m sequences_to_features -n "$NS" -f "$lib" -t "$F" \
                    -o "$out/${C}_annotated.xml" -blastn -m "$minlen" -np -ni $annoflags
              echo "== FEATURE 1: circuit + gate netlist =="
              "$PY" "$REPO/features_to_circuits/features_to_circuits.py" -n "$NS" -c "$lib" \
                    -t "$out/${C}_annotated.xml" -o "$out/${C}_circuit.xml" -m "$minlen" -gn
              echo "== FEATURE 2: truth table =="
              "$PY" "$REPO/features_to_circuits/circuit_to_truth_table.py" \
                    "$out/${C}_circuit_netlist.json" --yosys "$YOSYS"
          ) > "$out/$C.log" 2>&1 && { ok=$((ok+1)); printf '  %-34s ok\n' "$C"; } \
                                 || { fail=$((fail+1)); printf '  %-34s FAILED (see %s)\n' "$C" "$out/$C.log"; }
      done
      echo "== $ds: $ok ran, $fail failed =="
    )
}

CELLO_LIB="$HERE/cello/library/cello_library.xml"
MD5_LIB="$HERE/md5/library/Eco2C1G5T1_public_library.xml"
# Cello circuits are linear constructs of >=1889 bp, parts >=40 bp.
CELLO_ARGS=(cello "$CELLO_LIB" 1000 "-M 40")
# MD5 plasmids are circular, carry parts down to 14 bp, and need the 90% identity
# threshold the recorded run used: PJR1 is a 62 bp part that matches over 58 bp
# (93.5% coverage-weighted), so the 95% default drops it and the locus is annotated
# as the spurious constitutive Plambda instead. See MD5_validation_record.md.
MD5_ARGS=(md5 "$MD5_LIB" 2000 "-M 14 -cir -pid 90")

case "$MODE" in
    --quick) echo "### Cello, quick subset"
             # shellcheck disable=SC2086
             run_dataset "${CELLO_ARGS[@]}" "$OUT/cello" $QUICK_SET ;;
    --cello) echo "### Cello, all circuits"
             run_dataset "${CELLO_ARGS[@]}" "$OUT/cello" ;;
    --md5)   echo "### MD5 plasmids"
             run_dataset "${MD5_ARGS[@]}" "$OUT/md5" ;;
    --all)   echo "### Cello, all circuits"
             run_dataset "${CELLO_ARGS[@]}" "$OUT/cello"
             echo "### MD5 plasmids"
             run_dataset "${MD5_ARGS[@]}" "$OUT/md5" ;;
    *) echo "usage: $0 [--quick|--cello|--md5|--all] [output_dir]"; exit 2 ;;
esac

echo ""
echo "== now check the results: =="
[ -d "$OUT/cello" ] && echo "   $PY $HERE/check_results.py $OUT/cello --dataset cello"
[ -d "$OUT/md5" ]   && echo "   $PY $HERE/check_results.py $OUT/md5 --dataset md5"
