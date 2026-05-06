"""
fig5_regime2.py
===============
Figure 5 — Regime 2 quantitative analysis (3 panels).
  (a) Violin: c200/cLCDM for active vs inactive-strong
  (b) Per-galaxy Δχ² for the c200^BPB correction (27 galaxies)
  (c) δR1 → αc correlation (25 converged galaxies)

Data sources (frozen in results/chains/):
  results/chains/regime2_complete.csv   — per-galaxy Δχ², αc, δR1 (Regime 2)
  results/chains/sparc_active_sample.csv — c200/cLCDM for active galaxies
  results/chains/regime2_complete.csv    — c200/cLCDM for inactive strong

Output:
  results/figures/fig5_regime2.pdf
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT   = os.path.dirname(SCRIPT_DIR)
CHAINS      = os.path.join(REPO_ROOT, "results", "chains")
FIGURES     = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(FIGURES, exist_ok=True)

# ── Load data from CSVs ───────────────────────────────────────────────────────
r2   = pd.read_csv(os.path.join(CHAINS, "regime2_complete.csv"))
act  = pd.read_csv(os.path.join(CHAINS, "sparc_active_sample.csv"))
inact = pd.read_csv(os.path.join(CHAINS, "regime2_complete.csv"))

# Per-galaxy arrays (sorted by dR1, ascending → most negative first)
r2 = r2.sort_values("delta_R1").reset_index(drop=True)

dR1   = r2["delta_R1"].values
dchi2 = r2["delta_chi2"].values
alpha_c_raw = r2["alpha_c_opt"].values if "alpha_c_opt" in r2.columns else r2["alpha_c"].values
at_boundary = r2["at_boundary"].astype(bool).values.astype(bool) if "at_boundary" in r2.columns \
              else (np.abs(alpha_c_raw) > 18)

converged_mask = ~at_boundary
alpha_c = np.where(converged_mask, alpha_c_raw, np.nan)

conv_idx = np.where(converged_mask)[0]
n_conv   = converged_mask.sum()
n_bnd    = (~converged_mask).sum()

# c200/cLCDM
c_active       = act["c200_lcdm_ratio"].dropna().values
c_inact_strong = inact["c200_lcdm_ratio"].dropna().values

# t-test computed from data
from scipy.stats import ttest_ind as _ttest
_t, _p = _ttest(c_active, c_inact_strong, equal_var=False)
_d = ((c_inact_strong.mean() - c_active.mean()) /
      ((c_active.std(ddof=0)**2 + c_inact_strong.std(ddof=0)**2) / 2)**0.5)
ttest_str = f"$t={_t:.2f}$, $p={_p:.3f}$\nCohen $d={abs(_d):.2f}$"

# Totals
total_dchi2 = dchi2.sum()
n_improved  = (dchi2 < 0).sum()
pbinom_str  = "$p=7.45\\times10^{-9}$"

# Linear fit δR1 → αc (converged only)
dR1_conv = dR1[conv_idx]
ac_conv  = alpha_c[conv_idx]
fit      = np.polyfit(dR1_conv, ac_conv, 1)
r_fit    = np.corrcoef(dR1_conv, ac_conv)[0, 1]
from scipy.stats import pearsonr as _pr
_, p_fit = _pr(dR1_conv, ac_conv)
# Format: p = X.X × 10^{-N} in LaTeX mathtext
_exp  = int(f"{p_fit:.1e}".split("e")[1])
_mant = float(f"{p_fit:.1e}".split("e")[0])
p_fit_str = f"$p={_mant:.1f}\\times10^{{{_exp}}}$"

print(f"Total Δχ² = {total_dchi2:.0f}")
print(f"N improved: {n_improved}/27")
print(f"Converged: {n_conv}, Boundary: {n_bnd}")
print(f"Linear fit: αc = {fit[0]:.2f} × δR1 + {fit[1]:.2f}")
print(f"Correlation r = {r_fit:.3f}")

# ── Colour scheme ─────────────────────────────────────────────────────────────
C_BLUE = "#4C72B0"
C_RED  = "#C44E52"
FS     = 9
rng    = np.random.default_rng(7)

# ── Build figure ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10.5, 3.4))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

# ─── Panel (a): violin c200/cLCDM ────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0])

data_violin = [c_active, c_inact_strong]
labels_v    = [f"Active\n(rising, $n={len(c_active)}$)",
               f"Inactive strong\n(flat, $n={len(c_inact_strong)}$)"]
cols_v      = [C_BLUE, "#DD8452"]

parts = ax1.violinplot(data_violin, positions=[1, 2],
                       showmedians=True, showextrema=False, widths=0.5)
for pc, col in zip(parts["bodies"], cols_v):
    pc.set_facecolor(col); pc.set_alpha(0.65)
parts["cmedians"].set_color("k"); parts["cmedians"].set_linewidth(1.8)

for i, (grp, col) in enumerate(zip(data_violin, cols_v)):
    jx = rng.normal(i + 1, 0.04, len(grp))
    ax1.scatter(jx, grp, c=col, s=6, alpha=0.45, zorder=3)

ax1.axhline(1.0, color="k", ls="--", lw=1.0)
ax1.text(0.50, 0.97,
         ttest_str,
         transform=ax1.transAxes, ha="center", va="top", fontsize=7,
         bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))
ax1.set_xticks([1, 2]); ax1.set_xticklabels(labels_v, fontsize=7)
ax1.set_ylabel(r"$c_{200}/c_{200}^{\Lambda\mathrm{CDM}}$", fontsize=FS)
ax1.set_title("(a) Two-regime structure", fontsize=8)
ax1.tick_params(labelsize=7)

# ─── Panel (b): per-galaxy Δχ² ──────────────────────────────────────────────
ax2 = fig.add_subplot(gs[1])

CLIP          = -500
dchi2_clipped = np.clip(dchi2, CLIP, 0)

for i, (dc, is_conv) in enumerate(zip(dchi2_clipped, converged_mask)):
    ax2.bar(i, dc, width=0.75, color=C_BLUE if is_conv else C_RED, alpha=0.85)

ax2.axhline(0, color="k", lw=0.6)
ax2.set_xlabel(r"Galaxy (sorted by $\delta R_1$)", fontsize=FS)
ax2.set_ylabel(r"$\Delta\chi^2\ (\mathrm{BPB}-\mathrm{ref})$", fontsize=FS)
ax2.text(0.05, 0.05,
         f"{n_improved}/27 improved\n{pbinom_str}\n"
         f"$\\Delta\\chi^2_{{\\rm tot}}={total_dchi2:.0f}$",
         transform=ax2.transAxes, va="bottom", fontsize=7)
leg_h = [Patch(facecolor=C_BLUE, label=f"Converged ({n_conv})"),
         Patch(facecolor=C_RED,  label=f"Boundary ({n_bnd})")]
ax2.legend(handles=leg_h, fontsize=6.5, loc="lower right")
ax2.set_title(r"(b) $\chi^2$ improvement with $c_{200}^\mathrm{BPB}$", fontsize=8)
ax2.set_xlim(-0.5, 26.5); ax2.tick_params(labelsize=7)

# ─── Panel (c): δR1 → αc ─────────────────────────────────────────────────────
ax3 = fig.add_subplot(gs[2])

ax3.scatter(dR1_conv, ac_conv, c=C_BLUE, s=22, zorder=3,
            label=f"Converged ($n={n_conv}$)")
x_line = np.linspace(dR1_conv.min() - 0.003, dR1_conv.max() + 0.003, 100)
ax3.plot(x_line, np.polyval(fit, x_line), color=C_RED, lw=1.5, zorder=2)
ax3.text(0.97, 0.05,
         f"$r={r_fit:.3f}$, {p_fit_str}\n"
         f"$\\alpha_c={fit[0]:.1f}\\,\\delta R_1+{fit[1]:.1f}$",
         transform=ax3.transAxes, va="bottom", ha="right", fontsize=7)
ax3.set_xlabel(r"$\delta R_1$ (BPB correction)", fontsize=FS)
ax3.set_ylabel(r"$\alpha_c$ (optimal $c_{200}$ coupling)", fontsize=FS)
ax3.legend(fontsize=6.5, loc="upper right")
ax3.set_title(r"(c) $\delta R_1\to\alpha_c$ correlation", fontsize=8)
ax3.tick_params(labelsize=7)

fig.suptitle(r"Regime 2: BPB signal projection into $c_{200}$ (virialized halos)",
             fontsize=9, fontweight="bold", y=1.01)

outpath = os.path.join(FIGURES, "fig5_regime2.pdf")
fig.savefig(outpath, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {outpath}")
