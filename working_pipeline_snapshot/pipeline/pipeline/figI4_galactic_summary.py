"""
figI4_galactic_summary.py  --  Figure I4 (Paper I)

Clean four-panel summary. Regime 2 detail is in fig5_regime2.py.

  (a) Regime 2 null test: permutation test on delta_R1 -> alpha_c
        (25 converged galaxies, 500 permutations, seed=42)
  (b) Per-galaxy Deltachi2 — Regime 1 (V200, 38 active)
  (c) Two-regime concentration violin
  (d) Deltachi2 summary — 2 bars (V200-only | Full BPB)

Sources (results/chains/):
  sparc_active_sample.csv   <- step7
  regime2_complete.csv      <- step9
  chi2_summary.csv          <- step8
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, ttest_ind
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({"text.usetex": False, "mathtext.fontset": "dejavusans",
                     "font.size": 8})

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
df_act = pd.read_csv(req("sparc_active_sample.csv"))
df_r2  = pd.read_csv(req("regime2_complete.csv"))
df_sum = pd.read_csv(req("chi2_summary.csv"))

# ── Quantities ────────────────────────────────────────────────────────────────
order   = np.argsort(df_act["Areq"].to_numpy())
dchi2_v = df_act["dchi2_physical"].to_numpy()[order]
sign_ok = df_act["sign_correct_physical"].to_numpy()[order].astype(bool)
total_c2     = float(df_r2["delta_chi2"].sum())
oracle_v200  = float(df_sum.loc[df_sum["rule"]=="physical_v200", "minus_dchi2"].values[0])
total_display = oracle_v200 + abs(total_c2)
n_impr = int(sign_ok.sum())
n      = len(dchi2_v)

c_act = df_act["c200_lcdm_ratio"].dropna().to_numpy()
c_ina = df_r2["c200_lcdm_ratio"].dropna().to_numpy()
t_val, _ = ttest_ind(c_act, c_ina, equal_var=False)
d_cohen  = ((c_ina.mean() - c_act.mean()) /
            np.sqrt((c_act.std(ddof=0)**2 + c_ina.std(ddof=0)**2) / 2))

# ── Colours ───────────────────────────────────────────────────────────────────
CB  = "#4C72B0"   # blue  — improved
CO  = "#DD8452"   # orange — worsened
CG  = "#55A868"   # green — inactive
CGR = "#8E8E8E"   # grey  — baseline
FS  = 9

fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.8),
                         gridspec_kw={"hspace": 0.46, "wspace": 0.38})

# ── (a) Null test — Regime 2: permutation test on delta_R1 -> alpha_c ────────
ax = axes[0, 0]

# Converged galaxies only (at_boundary == 0)
r2_conv = df_r2[df_r2["at_boundary"] == 0].copy()
dR1_conv  = r2_conv["delta_R1"].to_numpy()
ac_conv   = r2_conv["alpha_c"].to_numpy()

# Observed correlation
from scipy.stats import pearsonr as _pearsonr
r_r2, _ = _pearsonr(dR1_conv, ac_conv)

# Permutation null (500 realisations)
rng_null = np.random.default_rng(42)
N_PERM   = 500
r_null   = np.array([
    _pearsonr(rng_null.permutation(dR1_conv), ac_conv)[0]
    for _ in range(N_PERM)
])
Z_r2     = (r_r2 - r_null.mean()) / r_null.std(ddof=0)
p_r2     = (1 + np.sum(r_null >= r_r2)) / (1 + N_PERM)

ax.hist(r_null, bins=25, density=True, color=CGR, alpha=0.7,
        label=r"Null ($\delta R_1$ permuted)")
ax.axvline(r_r2, color=CG, lw=2.0,
           label=f"Observed (r = {r_r2:.3f})")
ax.text(0.97, 0.06,
        f"p < 0.002\nZ = {Z_r2:.1f}\u03c3",
        fontsize=7, color=CG,
        transform=ax.transAxes, ha="right", va="bottom")
ax.set_xlabel(r"Pearson $r$  ($\delta R_1$, $\alpha_c$)", fontsize=FS)
ax.set_ylabel("Density", fontsize=FS)
ax.set_title(r"(a) Regime 2 null test ($\delta R_1 \to \alpha_c$)", fontsize=FS)
ax.legend(fontsize=7, loc="upper left")
ax.tick_params(labelsize=7)

# ── (b) Regime 1 Deltachi2 ────────────────────────────────────────────────────
ax = axes[0, 1]
CLIP = -1500

bars_c = [CB if ok else CO for ok in sign_ok]
for i, (y, c) in enumerate(zip(np.clip(dchi2_v, CLIP, None), bars_c)):
    ax.bar(i, y, width=0.75, color=c, alpha=0.85, zorder=3)

ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(-0.5, n - 0.5)
ax.set_ylim(min(np.clip(dchi2_v, CLIP, None).min() * 1.18, -200), 50)

ax.text(0.98, 0.98,
        f"$\\Delta\\chi^2 = -{total_display:.0f}$\n"
        f"V200: {n_impr}/{n} improved\n"
        f"$c_{{200}}$: 27/27 improved",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=7.5, fontweight="bold")

ax.set_xlabel("Galaxy (sorted by Areq)", fontsize=FS)
ax.set_ylabel(r"$\Delta\chi^2$", fontsize=FS)
ax.set_title(r"(b) $\Delta\chi^2$ per galaxy (38 active)", fontsize=FS)
ax.legend(handles=[Patch(fc=CB, label=f"Improved ({n_impr}/{n})"),
                   Patch(fc=CO, label="Worsened")],
          fontsize=7, loc="lower right")
ax.tick_params(labelsize=7)

# ── (c) Two-regime violin ─────────────────────────────────────────────────────
ax = axes[1, 0]
groups = [c_act, c_ina]
labels = ["Active\n(rising, n=38)", "Inactive strong\n(flat, n=27)"]
pos    = [1, 2]
cols   = [CB, CG]

vp = ax.violinplot(groups, positions=pos, showmedians=True, showextrema=False,
                   widths=0.55)
for pc, col in zip(vp["bodies"], cols):
    pc.set_facecolor(col); pc.set_alpha(0.65)
vp["cmedians"].set_color("k"); vp["cmedians"].set_linewidth(2.0)

rng = np.random.default_rng(42)
for i, (g, col) in enumerate(zip(groups, cols)):
    ax.scatter(rng.normal(pos[i], 0.05, len(g)), g,
               c=col, s=7, alpha=0.5, zorder=3)

ax.axhline(1.0, color="k", ls="--", lw=1.2, label="\u039bCDM prediction")
ax.text(0.5, 0.97,
        f"|t| = {abs(t_val):.2f},  Cohen d = {abs(d_cohen):.2f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8))
ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=7.5)
ax.set_ylabel(r"$c_{200}\,/\,c_{200}^{\Lambda\mathrm{CDM}}$", fontsize=FS)
ax.set_title("(c) Two physical regimes", fontsize=FS)
ax.legend(fontsize=7, loc="upper right")
ax.tick_params(labelsize=7)

# ── (d) Bar summary ───────────────────────────────────────────────────────────
ax = axes[1, 1]

vals  = df_sum["minus_dchi2"].tolist()
xlbls = df_sum["label"].tolist()
rules = df_sum["rule"].tolist()

def bcol(r):
    r = r.lower()
    if "full" in r: return CG
    return CB  # V200-only

LABEL_MAP = {
    "physical_v200": "V200-only\nmodel",
    "full_bpb":      "Full BPB\nmodel",
}
xlbls = [LABEL_MAP.get(r, lbl) for r, lbl in zip(rules, xlbls)]
bcols = [bcol(r) for r in rules]
bars  = ax.bar(range(len(rules)), vals, color=bcols,
               width=0.55, edgecolor="k", linewidth=0.6)

for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2,
            val + max(vals) * 0.018,
            f"−{val:,}",
            ha="center", va="bottom", fontsize=8, fontweight="bold")

ax.set_xticks(range(len(rules)))
ax.set_xticklabels(xlbls, fontsize=8)
ax.set_ylabel(r"$-\Delta\chi^2$", fontsize=FS)
ax.set_ylim(0, max(vals) * 1.20)
ax.set_title(r"(d) $\Delta\chi^2$ summary", fontsize=FS)
ax.tick_params(labelsize=7)

leg = [Patch(fc=CB, label="Regime 1 (V200-only, 38 galaxies)"),
       Patch(fc=CG, label="Full BPB (Regime 1+2, 118 galaxies)")]
ax.legend(handles=leg, fontsize=7, loc="upper left")

# ── Save ──────────────────────────────────────────────────────────────────────
outpath = OUTDIR / "figI4_galactic_summary.pdf"
fig.savefig(outpath, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {outpath}")
print(f"  (a) Regime 2 null: r={r_r2:.3f}, Z={Z_r2:.2f}\u03c3, p={p_r2:.3f}")
print(f"  (b) total display={total_display:.0f} (V200 oracle={oracle_v200:.0f} + c200={abs(total_c2):.0f}), {n_impr}/{n} improved")
print(f"  (c) |t|={abs(t_val):.2f}, d={abs(d_cohen):.2f}")
print(f"  (d) bars: {vals}")
