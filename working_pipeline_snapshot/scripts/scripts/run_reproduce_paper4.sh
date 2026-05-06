#!/usr/bin/env bash
# ============================================================
# Reproduce all Paper IV figures from the frozen all-phases chain
# Usage (from repo root): bash scripts/run_reproduce_paper4.sh
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

echo "=== Paper IV reproduction ==="

# Sanity check — frozen chain must be present
CHAIN="results/chains/bpb_allphases_chain.h5"
if [ ! -f "$CHAIN" ]; then
    echo "ERROR: $CHAIN not found at repo root."
    echo "       Download from Zenodo (doi:10.5281/zenodo.XXXXXXX) and place at repo root."
    exit 1
fi
echo "Chain found: $CHAIN ($(du -sh $CHAIN | cut -f1))"

echo ""
echo "[1/2] Galaxy-scale bridge figure (no chain required)..."
python pipeline/step9_build_regime2_complete.py
python pipeline/step7_build_sparc_active_sample.py
python pipeline/step8_build_chi2_summary.py
python pipeline/figIV_galaxy_bridge.py
echo "      -> results/figures/figIV_galaxy_bridge.pdf"

# Sanity check — frozen chain must be present
CHAIN="results/chains/bpb_allphases_chain.h5"
if [ ! -f "$CHAIN" ]; then
    echo ""
    echo "NOTE: $CHAIN not found — skipping Paper IV MCMC figures."
    echo "      Download from Zenodo (doi:10.5281/zenodo.XXXXXXX) to reproduce them."
    echo ""
    echo "Done. Available figure: results/figures/figIV_galaxy_bridge.pdf"
    exit 0
fi
echo ""
echo "Chain found: $CHAIN ($(du -sh $CHAIN | cut -f1))"

echo ""
echo "[2/2] All Paper IV MCMC figures (corner, diagnostics, lensing amplitude)..."
if python -c "import emcee, h5py, corner" 2>/dev/null; then
    python pipeline/figPaperIV_all.py
    echo "      -> results/figures/bpb_allphases_corner_color.pdf"
    echo "      -> results/figures/p3_Alens.pdf"
    echo "      -> results/figures/bpb_phase1_diagnostics.pdf"
else
    echo "  SKIP: emcee/h5py/corner not installed."
    echo "        Install with: pip install emcee corner h5py"
    echo "        Then rerun: python pipeline/figPaperIV_all.py"
fi

echo ""
echo "Done. Figures in: results/figures/"
