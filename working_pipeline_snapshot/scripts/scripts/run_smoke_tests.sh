#!/usr/bin/env bash
# ============================================================
# Smoke tests — verify all figure scripts run without error
# Usage (from repo root): bash scripts/run_smoke_tests.sh
# Expected runtime: ~2 min
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

PASS=0; FAIL=0

run_test() {
    local name="$1"; shift
    if "$@" > /dev/null 2>&1; then
        echo "  PASS  $name"
        PASS=$((PASS+1))
    else
        echo "  FAIL  $name"
        FAIL=$((FAIL+1))
    fi
}

echo "=== BPB/MBE smoke tests ==="
echo ""

echo "--- Paper I ---"
run_test "figI1_kids_xi"         python pipeline/figI1_kids_xi.py
run_test "figI2_chi2_bridge"     python pipeline/figI2_chi2_bridge.py
run_test "figI3_power_growth"    python pipeline/figI3_power_growth.py
run_test "figI4_galactic_summary" python pipeline/figI4_galactic_summary.py
run_test "figI5_chi2_landscapes" python pipeline/figI5_chi2_landscapes.py

echo ""
echo "--- Paper II ---"
run_test "fig1_bpb_signal"   python pipeline/fig1_bpb_signal.py results/chains/sparc_active_sample.csv
run_test "fig2_areq_scatter" python pipeline/fig2_areq_scatter.py results/chains/sparc_active_sample.csv
run_test "fig3_validation"   python pipeline/fig3_validation.py results/chains/sparc_active_sample.csv results/chains/regime2_complete.csv
run_test "fig4_chi2_landscapes" python pipeline/fig4_chi2_landscapes.py
run_test "fig5_regime2"      python pipeline/fig5_regime2.py

echo ""
echo "--- Paper V ---"
run_test "plot_corner_paperV" python pipeline/plot_corner_paperV.py \
    --npz results/chains/paperV/kids_kbf_fixed_dynesty.npz \
    --output_dir results/figures

echo ""
echo "--- Paper IV ---"
run_test "figIV_galaxy_bridge" python pipeline/figIV_galaxy_bridge.py
if python -c "import emcee, h5py" 2>/dev/null; then
    run_test "figPaperIV_all" python pipeline/figPaperIV_all.py
else
    echo "  SKIP  figPaperIV_all (emcee/h5py not installed)"
fi

echo ""
echo "======================================="
echo "  Results: ${PASS} passed, ${FAIL} failed"
echo "======================================="

[ "$FAIL" -eq 0 ] && exit 0 || exit 1
