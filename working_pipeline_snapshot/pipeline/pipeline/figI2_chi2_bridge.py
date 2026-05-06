"""
figI2_chi2_bridge.py  --  Figure A2: chi2(gamma_p) landscape + bridge

Data:
  results/chains/kids1000_chi2_scan.csv

Output:
  results/figures/figI2_chi2_bridge.pdf

Notes:
- Panel (a) is fully data-driven.
- The adopted best-fit is gamma_p=0.15 (physically motivated, not the
  absolute minimum which lies in the excluded undamped-memory regime).
- Panel (b) is a schematic of the BPB cosmological chain.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT   = Path(__file__).resolve().parents[1]
DATA   = ROOT / "results" / "chains" / "kids1000_chi2_scan.csv"
OUTDIR = ROOT / "results" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not DATA.exists():
    raise FileNotFoundError(f"Missing: {DATA}")

df = pd.read_csv(DATA)
if not {"gamma_p","dchi2"}.issubset(df.columns):
    raise ValueError("CSV must contain gamma_p and dchi2 columns")

gp    = df["gamma_p"].values
dchi2 = df["dchi2"].values

# Physically motivated best-fit: gamma_p=0.15 (NOT the absolute minimum)
# The absolute minimum at gamma_p~0.09-0.10 is in the undamped-memory regime
# and is excluded on physical grounds (see Paper I Section 3.2).
GP_ADOPTED  = 0.15
idx_adopted = np.argmin(np.abs(gp - GP_ADOPTED))
dchi2_adopted = dchi2[idx_adopted]

C_BPB = "#1f77b4"
FS    = 8

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2))

# ── Panel (a): chi2 landscape ─────────────────────────────────────────────
ax1.axhline(0, color="gray", ls="--", lw=0.8,
            label=r"$\Lambda$CDM reference ($\Delta\chi^2=0$)")
ax1.fill_between(gp, dchi2, 0, where=(dchi2 < 0),
                 alpha=0.15, color=C_BPB,
                 label=r"$\Delta\chi^2<0$ (BPB better)")
ax1.plot(gp, dchi2, "o-", color=C_BPB, ms=3, lw=1.3,
         label="BPB unified")

# Adopted best-fit marker (physically motivated)
ax1.axvline(GP_ADOPTED, color="k", ls=":", lw=1.0)
ax1.annotate(
    rf"$\gamma_p={GP_ADOPTED:.2f}$" + "\n" +
    rf"$\Delta\chi^2={dchi2_adopted:.2f}$",
    xy=(GP_ADOPTED, dchi2_adopted),
    xytext=(GP_ADOPTED + 0.06, dchi2_adopted - 0.25),
    fontsize=7,
    arrowprops=dict(arrowstyle="->", lw=0.8)
)

# Shade excluded undamped-memory regime (gamma_p < 0.13)
gp_excl = gp[gp < 0.13]
if len(gp_excl) > 0:
    ax1.axvspan(gp.min(), 0.13, alpha=0.07, color="red",
                label=r"Excluded ($\gamma_p\lesssim0.13$)")

ax1.set_xlabel(r"$\gamma_p$ (memory exponent)", fontsize=FS)
ax1.set_ylabel(r"$\Delta\chi^2_\mathrm{BPB} - \chi^2_{\Lambda\mathrm{CDM,best}}$",
               fontsize=FS)
ax1.set_title(r"(a) $\chi^2$ landscape vs $\gamma_p$ (KiDS-1000)", fontsize=FS)
ax1.legend(fontsize=6.5, loc="upper right")
ax1.tick_params(labelsize=7)

# ── Panel (b): cosmological bridge ────────────────────────────────────────
ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.axis("off")
ax2.set_title("(b) Cosmological bridge: KiDS → SPARC", fontsize=FS)

boxes = [
    (0.5, 0.82,
     r"$\mathbf{KiDS\text{-}1000\ \xi_\pm\ fit}$" + "\n"
     + rf"$\gamma_p={GP_ADOPTED:.2f},\ \lambda_p=1.0$" + "\n"
     + r"$H_0=68,\ \Omega_m=0.19,\ \Omega_{m,\mathrm{eff}}=0.352$",
     "#dce9f7"),
    (0.5, 0.52,
     r"$\mathbf{BPB\ transport}$" + "\n"
     + r"$R_\mathrm{BPB}(k)$: mean $= 1.011$" + "\n"
     + r"$R_\mathrm{eff}(k,a)=1+M_p(a)[R_\mathrm{BPB}-1]$" + "\n"
     + r"$\delta R_1(M)=R_\mathrm{eff}(M;A=1)-1$",
     "#ddf0dd"),
    (0.5, 0.16,
     r"$\mathbf{SPARC\ rotation\ curves}$" + "\n"
     + r"$\Delta\chi^2=-4108$ on 38 active gal." + "\n"
     + r"LOO $p=0.00083$, Cohen $d=1.15$",
     "#ddf0dd"),
]

for x, y, txt, col in boxes:
    ax2.text(x, y, txt, ha="center", va="center", fontsize=6.0,
             bbox=dict(boxstyle="round,pad=0.4", fc=col, ec="gray", lw=0.8),
             transform=ax2.transAxes)

for y0, y1 in [(0.64, 0.60), (0.37, 0.33)]:
    ax2.annotate("", xy=(0.5, y1), xytext=(0.5, y0),
                 xycoords="axes fraction", textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0))

plt.tight_layout()
out = OUTDIR / "figI2_chi2_bridge.pdf"
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")
print(f"Adopted best-fit: gamma_p={GP_ADOPTED:.2f}  dchi2={dchi2_adopted:.4f}")
