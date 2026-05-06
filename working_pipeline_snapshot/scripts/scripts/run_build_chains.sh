#!/usr/bin/env bash
# ============================================================
# run_build_chains.sh
# Regenerates all derived CSV artefacts from frozen raw inputs.
# Run this if you want to verify the full provenance chain
# from raw SPARC tables to the final per-galaxy summary tables.
#
# Usage (from repo root):
#   bash scripts/run_build_chains.sh
#
# Expected runtime: ~15-20 min (step3 runs 500 permutations)
# Step 5 requires emcee — skipped if not installed.
#
# Dependency graph:
#
#   [frozen raw inputs]
#     MassModels_Lelli2016c.mrt
#     sparc_reference_halo_table_nfw_lcdm_clean.csv
#     per_galaxy_v200_required_amplitude.csv
#     bpb_correction_proxy_grid_z0.csv
#     sparc_all_clean_inner_shape_proxy.csv
#     hybrid_loo_per_galaxy.csv
#     two_stage_areq_per_galaxy.csv
#
#   step1 → mouter_per_galaxy.csv
#   step2 (step1) → loo_predictions.csv
#   step3 (step2) → null_distribution.csv
#   step4 (step2) → delta_chi2_per_galaxy.csv
#   step6 → regime2_c200_scan.csv
#   step7 (step2,4) → sparc_active_sample.csv   ← Paper I/II figures
#   step8 (step7) → chi2_summary.csv             ← Paper I/II figures
#   step9 → regime2_complete.csv                 ← Paper II Figure 5
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

PASS=0; FAIL=0
run_step() {
    local name="$1"; local script="$2"
    echo "  [$name] $script ..."
    if python "$script"; then
        echo "        OK"
        PASS=$((PASS+1))
    else
        echo "        FAILED"
        FAIL=$((FAIL+1))
    fi
}

echo "=== Building derived chain artefacts ==="
echo ""

echo "--- Galaxy branch (Papers I & II) ---"
run_step "1/9" pipeline/step1_compute_mouter.py
run_step "2/9" pipeline/step2_sign_rule_loo.py
echo "  [3/9] pipeline/step3_null_test.py  (~10 min, 500 perms) ..."
if python pipeline/step3_null_test.py; then
    echo "        OK"; PASS=$((PASS+1))
else
    echo "        FAILED"; FAIL=$((FAIL+1))
fi
run_step "4/9" pipeline/step4_delta_chi2.py
run_step "4b/9" pipeline/step4b_type95_baseline.py
run_step "5/9" pipeline/step6_regime2_c200_scan.py
run_step "6/9" pipeline/step7_build_sparc_active_sample.py
run_step "7/9" pipeline/step8_build_chi2_summary.py
run_step "8/9" pipeline/step9_build_regime2_complete.py

echo ""
echo "--- Paper IV chain (requires emcee) ---"
if python -c "import emcee" 2>/dev/null; then
    echo "  [9/9] pipeline/bpb_mcmc_v2.py  (full MCMC ~5h) ..."
    echo "        Skipped in standard build — chain shipped as frozen artefact."
    echo "        To rerun: python pipeline/bpb_mcmc_v2.py"
else
    echo "  [9/9] emcee not installed — Paper IV chain is a frozen release artefact."
    echo "        Install with: pip install emcee h5py"
fi

echo ""
echo "======================================="
echo "  Results: $PASS passed, $FAIL failed"
echo "======================================="

[ "$FAIL" -eq 0 ] && exit 0 || exit 1
