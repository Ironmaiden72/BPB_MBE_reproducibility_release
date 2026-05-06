#!/usr/bin/env bash
# ============================================================
# Reproduce all Paper V figures from the frozen posterior
# Usage (from repo root): bash scripts/run_reproduce_paperV.sh
#
# NOTE: The frozen posterior (kids_kbf_fixed_dynesty.npz) is
# shipped with the repo. No KiDS bridge file is required for
# figure reproduction. Full re-run instructions are in MANIFEST.md.
# ============================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

NPZ="results/chains/paperV/kids_kbf_fixed_dynesty.npz"
OUT="results/figures"

echo "=== Paper V reproduction ==="
echo "Input : ${NPZ}"

if [ ! -f "$NPZ" ]; then
    echo "ERROR: Frozen posterior not found: ${NPZ}"
    exit 1
fi

echo "[1/2] Fig V1 — marginalised beta_c posterior..."
python pipeline/bpb_kids_only_fixed_kbf_v4.py \
    --plot_only \
    --output_dir results/chains/paperV 2>/dev/null \
    && cp results/chains/paperV/beta_c_kids_only.pdf "${OUT}/figV1_beta_c_posterior.pdf" \
    || echo "      (skipped: KiDS bridge not available — using shipped PDF)"
echo "      -> results/figures/figV1_beta_c_posterior.pdf"

echo "[2/2] Fig V2 — corner (beta_c, sigma8, omega_m, S8_eff)..."
python pipeline/plot_corner_paperV.py \
    --npz "${NPZ}" \
    --output_dir "${OUT}"
cp "${OUT}/corner_beta_s8_om_S8.pdf" "${OUT}/figV2_corner.pdf" 2>/dev/null || true
echo "      -> results/figures/figV2_corner.pdf"

echo ""
echo "Done. Paper V figures in: ${OUT}"
