"""
figI3_power_growth.py  --  Figure A3: C_ell^BPB/C_ell^ref + D_BPB(z)/D_ref(z)

Data (in results/chains/):
    bpb_power_spectrum.csv     -- growth ratio D_BPB/D_ref(z)

Output (in results/figures/):
    figI3_power_growth.pdf
    figI3_power_growth.png

Method — Panel (a):
    The C_ell ratio at unit amplitude is computed via a Limber approximation:

        C_ell^BPB / C_ell^ref (ell, bin_i) ≈ R_eff(ell) × [D_BPB(z_eff,i)/D_ref(z_eff,i)]²

    where:
    - z_eff,i is the effective lens redshift for each tomographic bin,
      approximated as 0.55 × z_s,i (lensing kernel peak);
    - z_s,i are the mean source redshifts for KiDS-1000 bins 1-5
      (Heymans et al. 2021);
    - R_eff(ell) = 1 + (R_eff_mean - 1) × exp(-ell/ell_break) captures
      the scale dependence: at large scales (low ell) the BPB primordial
      ratio R_BPB(k) produces maximum enhancement; at small scales
      (high ell) R_BPB → 1. Here R_eff_mean = 1.011 (Paper I §4.1) and
      ell_break = 500 (characteristic scale of BPB enhancement);
    - D_BPB/D_ref(z) is read directly from bpb_power_spectrum.csv.

    This approximation reproduces the 2-10 per cent enhancement shown in
    Paper I Figure A3 and is fully reproducible from stored derived data.

Panel (b): D_BPB(z)/D_ref(z) directly from CSV.

Usage (from repo root):
    python pipeline/figI3_power_growth.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT   = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

PS_CSV = ROOT / "results" / "chains" / "bpb_power_spectrum.csv"
if not PS_CSV.exists():
    raise FileNotFoundError(f"Missing: {PS_CSV}")

# ── Load growth ratio ───────────────────────────────────────────────────────
df_ps = pd.read_csv(PS_CSV)
gr    = df_ps[df_ps["data_type"] == "growth"].dropna(
            subset=["z","D_ratio"]).sort_values("z")
if gr.empty:
    raise ValueError("No growth rows in bpb_power_spectrum.csv")

z_g   = gr["z"].values
D_rat = gr["D_ratio"].values

# ── BPB parameters (Paper I Table A1) ──────────────────────────────────────
R_EFF_MEAN = 1.011   # mean R_BPB over KiDS ell grid (Paper I §4.1)
ELL_BREAK  = 500.0   # characteristic scale of BPB enhancement [arcmin^-1]

# ── KiDS-1000 source redshifts (Heymans et al. 2021, Table 1) ──────────────
# Mean source redshift per tomographic bin
Z_S_BINS = [0.35, 0.55, 0.77, 0.99, 1.27]

# Effective lens redshift: lensing kernel W(z) peaks at ~0.55 × z_s
Z_EFF_FACTOR = 0.55

# ── ell grid for C_ell ratio ────────────────────────────────────────────────
ell_arr = np.logspace(np.log10(80), np.log10(5100), 60)

# ── Style ──────────────────────────────────────────────────────────────────
FS     = 8
COLORS = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

# ── Panel (a): C_ell ratio per diagonal bin ────────────────────────────────
for i, z_s in enumerate(Z_S_BINS):
    z_eff = z_s * Z_EFF_FACTOR
    Dr    = float(np.interp(z_eff, z_g, D_rat))

    # C_ell ratio: growth enhancement × scale-dependent primordial ratio
    R_eff_ell = 1.0 + (R_EFF_MEAN - 1.0) * np.exp(-ell_arr / ELL_BREAK)
    cl_ratio  = R_eff_ell * Dr**2

    ax1.plot(ell_arr, cl_ratio, color=COLORS[i], lw=1.3,
             label=f"bin ({i+1},{i+1})")

ax1.axhline(1.0, color="k", ls="--", lw=0.8)
ax1.set_xscale("log")
ax1.set_xlabel(r"$\ell$", fontsize=FS)
ax1.set_ylabel(r"$C_\ell^\mathrm{BPB}/C_\ell^\mathrm{ref}$", fontsize=FS)
ax1.set_title(r"(a) BPB/reference $C_\ell$ ratio (diagonal bins)", fontsize=FS)
ax1.legend(fontsize=6.5, loc="upper right")
ax1.tick_params(labelsize=7)
ax1.set_xlim(80, 5100)
ax1.set_ylim(0.98, 1.42)

# ── Panel (b): growth ratio ────────────────────────────────────────────────
ax2.plot(z_g, D_rat, color=COLORS[0], lw=1.8, label="BPB growth")
ax2.axhline(1.0, color="#ff7f0e", ls="--", lw=1.3,
            label=r"$\Lambda$CDM ($D_\mathrm{ref}$)")
ax2.axvline(0.5, color="gray", ls=":", lw=0.8)
ax2.set_xlabel("Redshift $z$", fontsize=FS)
ax2.set_ylabel(r"$D_\mathrm{BPB}(z)/D_\mathrm{ref}(z)$", fontsize=FS)
ax2.set_title(r"(b) BPB growth enhancement vs $z$", fontsize=FS)
ax2.legend(fontsize=7, loc="upper left")
ax2.tick_params(labelsize=7)
ax2.set_xlim(0, 3)
ax2.set_ylim(0.98, 1.80)

plt.tight_layout()

out_pdf = OUTDIR / "figI3_power_growth.pdf"
out_png = OUTDIR / "figI3_power_growth.png"
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
print()
for i, z_s in enumerate(Z_S_BINS):
    z_eff = z_s * Z_EFF_FACTOR
    Dr    = float(np.interp(z_eff, z_g, D_rat))
    r_low = (1.0 + (R_EFF_MEAN-1.0)*np.exp(-80/ELL_BREAK)) * Dr**2
    r_hi  = (1.0 + (R_EFF_MEAN-1.0)*np.exp(-5100/ELL_BREAK)) * Dr**2
    print(f"Bin ({i+1},{i+1}): z_eff={z_eff:.2f}  D_ratio={Dr:.3f}  "
          f"C_ell ratio [{r_hi:.3f}→{r_low:.3f}]")
