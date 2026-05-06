"""
figI1_kids_xi.py  --  Figure A1: KiDS-1000 xi+/xi- vs BPB(gp=0.15) and LCDM

Data:
  results/chains/kids1000_xi_data.csv

Output:
  results/figures/figI1_kids_xi.pdf

Notes:
- This script is fully reproducible from input CSV.
- The displayed Δχ² value is the paper value (not recomputed here).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Plot KiDS xi+/xi- comparison")
parser.add_argument("--data", type=str, default=None,
                    help="Path to kids1000_xi_data.csv")
parser.add_argument("--outdir", type=str, default=None,
                    help="Output directory for figures")
args = parser.parse_args()

# ─────────────────────────────────────────────────────────────────────────────
# Paths (root-safe)
# ─────────────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[1]

DATA = Path(args.data) if args.data else ROOT / "results" / "chains" / "kids1000_xi_data.csv"
OUTDIR = Path(args.outdir) if args.outdir else ROOT / "results" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not DATA.exists():
    raise FileNotFoundError(f"Missing input CSV: {DATA}")

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA)

# ─────────────────────────────────────────────────────────────────────────────
# Plot settings
# ─────────────────────────────────────────────────────────────────────────────
C_BPB = "#1f77b4"
C_LCDM = "#ff7f0e"
FS = 8

panels_plus  = [(1,1),(2,3),(5,5)]
panels_minus = [(1,2),(2,4),(4,5)]

fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.0))
fig.suptitle(r"KiDS-1000 $\xi_+$ and $\xi_-$: BPB/MBE vs $\Lambda$CDM", fontsize=9)

# ─────────────────────────────────────────────────────────────────────────────
# xi+
# ─────────────────────────────────────────────────────────────────────────────
for col, (b1,b2) in enumerate(panels_plus):
    ax = axes[0, col]
    sub = df[(df.xi_type=="xip")&(df.bin1==b1)&(df.bin2==b2)].sort_values("theta")

    if sub.empty:
        raise ValueError(f"No data for xi+ bins ({b1},{b2})")

    th = sub["theta"].values
    scale = th**2 * 1e4 # type: ignore

    ax.plot(th, sub["xi_bpb"]*scale, color=C_BPB, lw=1.3,
            label=r"BPB ($\gamma_p=0.15$)")
    ax.plot(th, sub["xi_lcdm"]*scale, color=C_LCDM, lw=1.3, ls="--",
            label=r"$\Lambda$CDM")

    ax.errorbar(th, sub["xi_data"]*scale,
                yerr=sub["xi_err"]*scale,
                fmt="o", ms=4, color="k",
                elinewidth=0.7, capsize=2, zorder=3,
                label="KiDS-1000")

    ax.set_xscale("log")
    ax.set_title(rf"$\xi_+$: bins ({b1},{b2})", fontsize=FS)
    ax.set_ylabel(r"$\theta^2\xi_+\,[\times10^{-4}]$", fontsize=FS-1)
    ax.tick_params(labelsize=FS-1)

    if col == 0:
        ax.legend(fontsize=6, loc="lower left")

# ─────────────────────────────────────────────────────────────────────────────
# xi-
# ─────────────────────────────────────────────────────────────────────────────
for col, (b1,b2) in enumerate(panels_minus):
    ax = axes[1, col]
    sub = df[(df.xi_type=="xim")&(df.bin1==b1)&(df.bin2==b2)].sort_values("theta")

    if sub.empty:
        raise ValueError(f"No data for xi- bins ({b1},{b2})")

    th = sub["theta"].values
    scale = th**2 * 1e4 # type: ignore

    ax.plot(th, sub["xi_bpb"]*scale, color=C_BPB, lw=1.3)
    ax.plot(th, sub["xi_lcdm"]*scale, color=C_LCDM, lw=1.3, ls="--")

    ax.errorbar(th, sub["xi_data"]*scale,
                yerr=sub["xi_err"]*scale,
                fmt="s", ms=4, color="k",
                elinewidth=0.7, capsize=2, zorder=3)

    ax.set_xscale("log")
    ax.set_title(rf"$\xi_-$: bins ({b1},{b2})", fontsize=FS)
    ax.set_xlabel(r"$\theta$ [arcmin]", fontsize=FS-1)
    ax.set_ylabel(r"$\theta^2\xi_-\,[\times10^{-4}]$", fontsize=FS-1)
    ax.tick_params(labelsize=FS-1)

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
out = OUTDIR / "figI1_kids_xi.pdf"
plt.tight_layout(rect=[0,0.01,1,1]) # type: ignore
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out}")