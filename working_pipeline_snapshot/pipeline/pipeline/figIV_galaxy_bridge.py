"""
figIV_galaxy_bridge.py  --  Figure for Paper IV (galaxy-scale bridge)

Condensed two-panel figure for Paper IV, adapted from Paper I Fig. I4.
Shows only the summary panels — statistical validation and per-galaxy
detail are in Paper I.

  Left  (a): Two-regime concentration structure  (c200/cLCDM)
  Right (b): Deltachi2 summary — V200-only | Full BPB

Narration (Paper IV context):
  "The galaxy-scale two-regime structure is consistent with the global
   BPB framework calibrated against late-universe probes."

Caption note:
  "Adapted from Paper I Fig. I4. Statistical validation and per-galaxy
   diagnostics are in Paper I."

Sources (results/chains/):
  regime2_complete.csv      <- step9
  sparc_active_sample.csv   <- step7
  chi2_summary.csv          <- step8

Output: results/figures/figIV_galaxy_bridge.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({"text.usetex": False, "mathtext.fontset": "dejavusans",
                     "font.size": 9})

ROOT   = Path(__file__).resolve().parents[1]
CHAINS = ROOT / "results" / "chains"
OUTDIR = ROOT / "results" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

def req(name):
    p = CHAINS / name
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")
    return p

# ── Load ──────────────────────────────────────────────────────────────────────
df_r2  = pd.read_csv(req("regime2_complete.csv"))
df_act = pd.read_csv(req("sparc_active_sample.csv"))
df_sum = pd.read_csv(req("chi2_summary.csv"))

# ── Quantities ────────────────────────────────────────────────────────────────
c_act = df_act["c200_lcdm_ratio"].dropna().to_numpy()
c_ina = df_r2["c200_lcdm_ratio"].dropna().to_numpy()
t_val, _ = ttest_ind(c_act, c_ina, equal_var=False)
d_cohen  = ((c_ina.mean() - c_act.mean()) /
            np.sqrt((c_act.std(ddof=0)**2 + c_ina.std(ddof=0)**2) / 2))

rules = df_sum["rule"].tolist()
vals  = df_sum["minus_dchi2"].tolist()

total_c2    = float(df_r2["delta_chi2"].sum())
oracle_v200 = float(df_sum.loc[df_sum["rule"]=="physical_v200", "minus_dchi2"].values[0])
total_full  = oracle_v200 + abs(total_c2)

# ── Colours ───────────────────────────────────────────────────────────────────
CB  = "#4C72B0"   # blue  — Regime 1 (V200)
CG  = "#55A868"   # green — Full BPB
CGR = "#8E8E8E"   # grey
FS  = 10

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4),
                         gridspec_kw={"wspace": 0.38})

# ── (a) Two-regime violin ─────────────────────────────────────────────────────
ax = axes[0]
groups = [c_act, c_ina]
labels = ["Active\n(rising, n=38)", "Inactive strong\n(flat, n=27)"]
pos    = [1, 2]
cols   = [CB, CG]

vp = ax.violinplot(groups, positions=pos, showmedians=True,
                   showextrema=False, widths=0.55)
for pc, col in zip(vp["bodies"], cols):
    pc.set_facecolor(col); pc.set_alpha(0.65)
vp["cmedians"].set_color("k"); vp["cmedians"].set_linewidth(2.0)

rng = np.random.default_rng(42)
for i, (g, col) in enumerate(zip(groups, cols)):
    ax.scatter(rng.normal(pos[i], 0.05, len(g)), g,
               c=col, s=8, alpha=0.5, zorder=3)

ax.axhline(1.0, color="k", ls="--", lw=1.2)
ax.text(0.5, 0.97,
        f"|t| = {abs(t_val):.2f},  Cohen d = {abs(d_cohen):.2f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", alpha=0.85))

ax.set_xticks(pos)
ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel(r"$c_{200}\,/\,c_{200}^{\Lambda\mathrm{CDM}}$", fontsize=FS)
ax.set_title("(a) Two-regime concentration", fontsize=FS)
ax.tick_params(labelsize=8)

# ΛCDM annotation
ax.text(1.5, 1.04, r"$\Lambda$CDM", ha="center", va="bottom",
        fontsize=7.5, color="k")

# ── (b) Deltachi2 summary ────────────────────────────────────────────────────
ax = axes[1]

LABEL_MAP = {
    "physical_v200": "V200-only\n(Regime 1)",
    "full_bpb":      "Full BPB\n(Regime 1+2)",
}
xlbls = [LABEL_MAP.get(r, r) for r in rules]
bcols  = [CG if "full" in r.lower() else CB for r in rules]

bars = ax.bar(range(len(rules)), vals, color=bcols,
              width=0.50, edgecolor="k", linewidth=0.6)

for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2,
            val + max(vals) * 0.020,
            f"−{val:,}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")

ax.set_xticks(range(len(rules)))
ax.set_xticklabels(xlbls, fontsize=8.5)
ax.set_ylabel(r"$-\Delta\chi^2$", fontsize=FS)
ax.set_ylim(0, max(vals) * 1.22)
ax.set_title(r"(b) $\Delta\chi^2$ summary (118 galaxies)", fontsize=FS)
ax.tick_params(labelsize=8)

leg = [Patch(fc=CB, label="Regime 1 (V200, 38 active)"),
       Patch(fc=CG, label="Full BPB (V200+c\u2082\u2080\u2080, 118 gal.)")]
ax.legend(handles=leg, fontsize=7.5, loc="upper left")

# ── Caption note ──────────────────────────────────────────────────────────────
fig.text(0.5, -0.03,
         "Adapted from Paper I Fig. I4. Statistical validation and per-galaxy "
         "diagnostics are in Paper I.",
         ha="center", va="top", fontsize=7, style="italic", color="grey")

# ── Save ──────────────────────────────────────────────────────────────────────
outpath = OUTDIR / "figIV_galaxy_bridge.pdf"
fig.savefig(outpath, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {outpath}")
print(f"  (a) |t|={abs(t_val):.2f}, d={abs(d_cohen):.2f}")
print(f"  (b) bars: {vals}")
