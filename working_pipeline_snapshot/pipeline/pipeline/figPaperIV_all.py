"""
figPaperIV_all.py
==================
Generate all Paper IV figures from the archived MCMC chain.

Figures produced in results/figures/:
    bpb_allphases_corner_color.pdf   — Figure 1: all-phases posterior corner plot
    bpb_phase1_diagnostics.pdf       — Figure 2: H(z), fσ8(z), activation function
    p3_Alens.pdf                     — Figure 3: CMB lensing amplitude

Input (at repo root):
    bpb_allphases_chain.h5           — archived MCMC chain (92 MB)

Requirements:
    pip install emcee corner h5py matplotlib scipy numpy

Usage (from repo root):
    python pipeline/figPaperIV_all.py
"""

import os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
FIGURES    = os.path.join(REPO_ROOT, "results", "figures")
CHAIN_FILE = os.path.join(REPO_ROOT, "results", "chains", "bpb_allphases_chain.h5")
# Fallback: repo root (legacy location)
if not os.path.exists(CHAIN_FILE):
    CHAIN_FILE = os.path.join(REPO_ROOT, "bpb_allphases_chain.h5")
CHAIN_P3   = os.path.join(REPO_ROOT, "results", "chains", "chain_p3.h5")

os.makedirs(FIGURES, exist_ok=True)
sys.path.insert(0, SCRIPT_DIR)

# Check dependencies
try:
    import emcee, corner, h5py
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install emcee corner h5py")
    sys.exit(1)

# Check chain
if not os.path.exists(CHAIN_FILE):
    raise FileNotFoundError(
        f"Missing MCMC chain: {CHAIN_FILE}\n"
        "This file is required to reproduce Paper IV figures.")

print(f"Using chain: {CHAIN_FILE}")

# ── Figure 1: all-phases corner plot (colour) ──────────────────────────────
print("\n--- Figure 1: corner plot ---")
from make_corner_color import make_corner_color
out_corner = os.path.join(FIGURES, "bpb_allphases_corner_color.pdf")
make_corner_color(CHAIN_FILE, outfile=out_corner)

# ── Figure 2: Phase 1 diagnostics (H(z), fσ8, activation) ────────────────
print("\n--- Figure 2: Phase 1 diagnostics ---")
from bpb_mcmc_v2 import make_diagnostics_plot
import numpy as np

reader   = emcee.backends.HDFBackend(CHAIN_FILE, read_only=True)
n_burn    = max(0, reader.iteration // 10)
chain_flat = reader.get_chain(flat=True, discard=n_burn, thin=3)
log_prob   = reader.get_log_prob(flat=True, discard=n_burn, thin=3)
# MAP = argmax log_prob (exact same point used for Paper IV Fig 2)
theta_map  = chain_flat[np.argmax(log_prob)]
# Save MAP for reproducibility
map_path = os.path.join(REPO_ROOT, "results", "chains", "bpb_allphases_map.npy")
np.save(map_path, theta_map)
print(f"MAP saved: {map_path}")
print(f"MAP params: H0={theta_map[2]:.2f}  beta_c={theta_map[6]:.4f}  sigma8={theta_map[7]:.4f}")
out_diag = os.path.join(FIGURES, "bpb_phase1_diagnostics.pdf")
make_diagnostics_plot(theta_map, out=out_diag)

# ── Figure 3: Phase 3 Alens ───────────────────────────────────────────────
print("\n--- Figure 3: CMB lensing Alens ---")
if os.path.exists(CHAIN_P3):
    from bpb_mcmc_phase3 import make_plots as make_p3_plots
    reader_p3  = emcee.backends.HDFBackend(CHAIN_P3, read_only=True)
    n_burn_p3  = max(0, reader_p3.iteration // 10)
    chain_p3   = reader_p3.get_chain(flat=True, discard=n_burn_p3, thin=3)
    lp_p3      = reader_p3.get_log_prob(flat=True, discard=n_burn_p3, thin=3)
    out_prefix = os.path.join(FIGURES, "p3")
    make_p3_plots(chain_p3, lp_p3, out_prefix=out_prefix)
else:
    print(f"  Phase 3 chain not found: {CHAIN_P3}")
    print(f"  Skipping Figure 3 — copy chain_p3.h5 to results/chains/ to generate.")

print("\nAll Paper IV figures saved to:", FIGURES)
