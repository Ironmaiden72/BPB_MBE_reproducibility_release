#!/usr/bin/env bash
# ============================================================
# Reproduce all Paper I figures from frozen derived products
# Usage (from repo root): bash scripts/run_reproduce_paper1.sh
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

echo "=== Paper I reproduction ==="

echo "[1/5] Fig I1 — KiDS xi+/xi- vs BPB/LCDM..."
python pipeline/figI1_kids_xi.py
echo "      -> results/figures/figI1_kids_xi.pdf"

echo "[2/5] Fig I2 — chi2(gamma_p) landscape + bridge..."
python pipeline/figI2_chi2_bridge.py
echo "      -> results/figures/figI2_chi2_bridge.pdf"

echo "[3/5] Fig I3 — C_ell^BPB/C_ell^ref + growth ratio..."
python pipeline/figI3_power_growth.py
echo "      -> results/figures/figI3_power_growth.pdf"

echo "[4/5] Fig I4 — galactic summary (two-regime)..."
python pipeline/figI4_galactic_summary.py
echo "      -> results/figures/figI4_galactic_summary.pdf"

echo "[5/5] Fig I5 — chi2(A) landscapes (IC2574, DDO170)..."
python pipeline/figI5_chi2_landscapes.py
echo "      -> results/figures/figI5_chi2_landscapes.pdf"

echo ""
echo "Done. All Paper I figures in: results/figures/"
