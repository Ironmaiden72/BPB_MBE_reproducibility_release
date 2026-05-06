#!/usr/bin/env bash
# ============================================================
# Reproduce all Paper II figures from frozen derived products
# Usage (from repo root): bash scripts/run_reproduce_paper2.sh
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

echo "=== Paper II reproduction ==="

echo "[0/5] Rebuilding derived CSV inputs..."
python pipeline/step9_build_regime2_complete.py
python pipeline/step7_build_sparc_active_sample.py
python pipeline/step8_build_chi2_summary.py
echo "      -> regime2_complete.csv, sparc_active_sample.csv, chi2_summary.csv"
echo ""

echo "[1/5] Fig 1 — BPB signal (delta_R1 histogram + decision space)..."
python pipeline/fig1_bpb_signal.py \
    results/chains/sparc_active_sample.csv
echo "      -> results/figures/fig1a_signal_hist.pdf"
echo "      -> results/figures/fig1b_decision_space.pdf"

echo "[2/5] Fig 2 — A_req scatter (physical rule vs oracle)..."
python pipeline/fig2_areq_scatter.py \
    results/chains/sparc_active_sample.csv
echo "      -> results/figures/fig2a_areq_physical.pdf"
echo "      -> results/figures/fig2b_areq_oracle.pdf"

echo "[3/5] Fig 3 — statistical validation (null, per-galaxy, regimes, summary)..."
python pipeline/fig3_validation.py \
    results/chains/sparc_active_sample.csv \
    results/chains/regime2_complete.csv
echo "      -> results/figures/fig3a_null.pdf ... fig3d_chi2summary.pdf"

echo "[4/5] Fig 4 — chi2(A) landscapes..."
python pipeline/fig4_chi2_landscapes.py
echo "      -> results/figures/fig4a_IC2574.pdf"
echo "      -> results/figures/fig4b_DDO170.pdf"

echo "[5/5] Fig 5 — Regime 2 quantitative analysis..."
python pipeline/fig5_regime2.py
echo "      -> results/figures/fig5_regime2.pdf"

echo ""
echo "Done. All Paper II figures in: results/figures/"
