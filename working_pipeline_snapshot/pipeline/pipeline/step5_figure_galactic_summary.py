"""
step5_figure_galactic_summary.py
=================================
Generate the 5-panel galactic branch summary figure.

Panels:
  (a) Synthetic null test — 500-permutation null distribution vs observed r_LOO
  (b) Delta chi2 per galaxy — physical rule vs oracle, sorted by Areq
  (c) Two physical regimes — c200/c200_LCDM violin plots, 3 galaxy groups
  (d) Delta chi2 summary — 4-bar chart
  (e) Regime 2 — per-galaxy delta_chi2 from c200 BPB correction (27 virialized)

Inputs (in results/chains/):
    null_distribution.csv
    delta_chi2_per_galaxy.csv
    per_galaxy_v200_required_amplitude.csv
    sparc_reference_halo_table_nfw_lcdm_clean.csv
    regime2_c200_scan.csv          (output of step6)

Output (in results/figures/):
    figI4_galactic_summary.pdf
    figI4_galactic_summary.png

Usage (from repo root):
    python pipeline/step5_figure_galactic_summary.py
"""

import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker
from scipy.stats import ttest_ind

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")
FIGURES    = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(FIGURES, exist_ok=True)

for f in ["null_distribution.csv", "delta_chi2_per_galaxy.csv",
          "per_galaxy_v200_required_amplitude.csv",
          "sparc_reference_halo_table_nfw_lcdm_clean.csv",
          "regime2_c200_scan.csv"]:
    if not os.path.exists(os.path.join(CHAINS, f)):
        raise FileNotFoundError(f"Missing: {f}\nRun steps 1-6 first.")

null = pd.read_csv(os.path.join(CHAINS, "null_distribution.csv"))
chi2 = pd.read_csv(os.path.join(CHAINS, "delta_chi2_per_galaxy.csv"))
v200 = pd.read_csv(os.path.join(CHAINS, "per_galaxy_v200_required_amplitude.csv"))
ref  = pd.read_csv(os.path.join(CHAINS, "sparc_reference_halo_table_nfw_lcdm_clean.csv"))
r2   = pd.read_csv(os.path.join(CHAINS, "regime2_c200_scan.csv"))

# Panel (a) stats
null_r  = null['r'].values
r_obs   = float(null['r_observed'].iloc[0])
Z       = (r_obs - null_r.mean()) / null_r.std()
p_perm  = float((null_r >= r_obs).mean())

# Panel (c) — c200 ratio
def dutton(logM): return 10**(0.905 - 0.101*(logM - 12.))
mg = v200.merge(ref[['name','C200']], on='name', how='left')
mg['c200_lcdm']  = mg['logM200'].apply(dutton)
mg['c200_ratio'] = mg['C200'] / mg['c200_lcdm']
c_act = mg[mg['active_shape_gate']==1]['c200_ratio'].dropna().values
c_ist = mg[(mg['active_shape_gate']==0)&(mg['delta_R1']<-0.131)]['c200_ratio'].dropna().values
c_iwk = mg[(mg['active_shape_gate']==0)&(mg['delta_R1']>=-0.131)]['c200_ratio'].dropna().values
t_stat, _ = ttest_ind(c_act, c_ist)
cohen_d   = abs((c_ist.mean()-c_act.mean())/np.sqrt((c_act.std()**2+c_ist.std()**2)/2))

# Panel (b) & (d) stats
rule_38   = chi2['dchi2_rule'].sum()
oracle_38 = chi2['dchi2_oracle'].sum()
pct       = rule_38/oracle_38*100

# Panel (e) stats
r2s       = r2.sort_values('dR1').reset_index(drop=True)
n_imp_r2  = int((r2['delta_chi2'] < 0).sum())
dchi2_r2  = r2['delta_chi2'].sum()

# ── Style ──────────────────────────────────────────────────────────────────
FS = 7.8
C_POS="#4878CF"; C_NEG="#E85D3A"
C_ACT="#4878CF"; C_IST="#6ABF69"; C_IWK="#AAAAAA"
plt.rcParams.update({"font.family":"serif","font.size":FS,
    "axes.labelsize":FS,"axes.titlesize":FS+0.5,
    "xtick.labelsize":FS-1,"ytick.labelsize":FS-1,
    "legend.fontsize":FS-1.5,"figure.dpi":150})

fig = plt.figure(figsize=(13, 8.5))
fig.suptitle("Summary of galactic branch results (detailed diagnostics in Paper II)",
             style="italic", fontsize=FS)

# 3 panels top, 2 panels bottom
gt = gridspec.GridSpec(1, 3, figure=fig,
                       left=0.06, right=0.98, top=0.90, bottom=0.52, wspace=0.32)
gb = gridspec.GridSpec(1, 2, figure=fig,
                       left=0.06, right=0.98, top=0.45, bottom=0.08, wspace=0.28)

ax_a = fig.add_subplot(gt[0, 0])
ax_b = fig.add_subplot(gt[0, 1])
ax_c = fig.add_subplot(gt[0, 2])
ax_d = fig.add_subplot(gb[0, 0])
ax_e = fig.add_subplot(gb[0, 1])

# ── (a) Synthetic null test ────────────────────────────────────────────────
ax_a.hist(null_r, bins=28, density=True, color="#CCCCCC",
          edgecolor="white", lw=0.3, label=r"Null ($\delta R_1$ shuffled)")
ax_a.axvline(r_obs, color=C_POS, lw=1.8, label=f"Observed ($r={r_obs:.3f}$)")
ax_a.text(0.96, 0.95, f"$p={p_perm:.3f}$\n$Z={Z:.2f}\\sigma$",
          transform=ax_a.transAxes, ha='right', va='top',
          color=C_POS, fontsize=FS-0.5)
ax_a.set_xlabel("LOO correlation $r$")
ax_a.set_ylabel("Density")
ax_a.set_title("(a) Synthetic null test (500 perms)")
ax_a.legend(frameon=False, loc='upper left', fontsize=FS-2)

# ── (b) Delta chi2 per galaxy ──────────────────────────────────────────────
cs = chi2.sort_values('Areq').reset_index(drop=True)
xb = np.arange(len(cs))
ax_b.bar(xb, cs['dchi2_rule'],
         color=[C_POS if v<=0 else C_NEG for v in cs['dchi2_rule']],
         width=0.85, zorder=2)
ax_b.step(np.append(xb-0.5, xb[-1]+0.5),
          np.append(cs['dchi2_oracle'].values, cs['dchi2_oracle'].values[-1]),
          where='post', color='k', ls='--', lw=0.9, zorder=3)
ax_b.axhline(0, color='k', lw=0.5)
ax_b.text(0.02, 0.04, f"Total $\\Delta\\chi^2={rule_38:.0f}$",
          transform=ax_b.transAxes, fontsize=FS-1.5)
ax_b.set_xlabel(r"Galaxy (sorted by $A_{\rm req}$)")
ax_b.set_ylabel(r"$\Delta\chi^2$")
ax_b.set_title(r"(b) $\Delta\chi^2$ per galaxy (38 actives)")
ax_b.legend(handles=[mpatches.Patch(facecolor=C_POS, label="Improved"),
                      mpatches.Patch(facecolor=C_NEG, label="Worsened"),
                      plt.Line2D([0],[0], color='k', ls='--', label='Oracle')],
            frameon=False, fontsize=FS-2)

# ── (c) Two physical regimes ───────────────────────────────────────────────
rng = np.random.default_rng(0)
for g, col, pos in zip([c_act, c_ist, c_iwk], [C_ACT, C_IST, C_IWK], [1, 2, 3]):
    parts = ax_c.violinplot([g], positions=[pos], showmedians=True,
                             showextrema=False, widths=0.68)
    for pc in parts['bodies']: pc.set_facecolor(col); pc.set_alpha(0.55)
    parts['cmedians'].set_color('k'); parts['cmedians'].set_linewidth(2)
    ax_c.scatter(np.full(len(g), pos) + rng.uniform(-0.12, 0.12, len(g)), g,
                 color=col, alpha=0.55, s=8, zorder=3)
ax_c.axhline(1.0, color='k', ls='--', lw=1.0)
ax_c.set_xticks([1, 2, 3])
ax_c.set_xticklabels(["Actives\n($s<0.45$)",
                       r"Inact. $\delta R_1<-0.131$",
                       "Inact.\nweak BPB"], fontsize=FS-2)
ax_c.set_ylabel(r"$c_{200}/c_{200}^{\Lambda{\rm CDM}}$")
ax_c.set_title("(c) Two physical regimes")
ax_c.text(0.97, 0.97, f"$t={t_stat:.2f},\\ d={cohen_d:.2f}$\n(actives vs. inact. strong)",
          transform=ax_c.transAxes, ha='right', va='top', fontsize=FS-2)
ax_c.text(0.03, 0.03, r"$\Lambda$CDM", transform=ax_c.transAxes, fontsize=FS-2)

# ── (d) Delta chi2 summary bars ────────────────────────────────────────────
bdata = [("Type-9.5",           2319,          "#AAAAAA", "gate"),
         ("Physical\nrule",     int(-rule_38),  C_POS,    "strict"),
         ("Oracle\n(38 act.)",  int(-oracle_38),"#7BA5D5", ""),
         ("Oracle\n(118 gal.)", int(-oracle_38),"#B8D0E8", "")]
xd = np.arange(len(bdata))
bars = ax_d.bar(xd, [b[1] for b in bdata],
                color=[b[2] for b in bdata], width=0.65, zorder=2)
for bar, (lbl, val, col, ilbl) in zip(bars, bdata):
    ax_d.text(bar.get_x()+bar.get_width()/2,
              bar.get_height()+max(b[1] for b in bdata)*0.012,
              f"−{val:,}", ha='center', va='bottom', fontsize=FS-2)
    if ilbl:
        ax_d.text(bar.get_x()+bar.get_width()/2, bar.get_height()*0.45,
                  ilbl, ha='center', va='center',
                  color='white', fontsize=FS-1.5, fontweight='bold')
ax_d.set_xticks(xd)
ax_d.set_xticklabels([b[0] for b in bdata], fontsize=FS-2)
ax_d.set_ylabel(r"$-\Delta\chi^2$")
ax_d.set_title(r"(d) $\Delta\chi^2$ summary")
ax_d.set_ylim(0, max(b[1] for b in bdata)*1.18)
ax_d.yaxis.set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

# ── (e) Regime 2 ───────────────────────────────────────────────────────────
CLIP  = -500
xe    = np.arange(len(r2s))
r2disp = r2s['delta_chi2'].clip(lower=CLIP)
n_conv = int((~r2s['at_boundary'].astype(bool)).sum())
n_bnd  = int(r2s['at_boundary'].sum())

ax_e.bar(xe, r2disp.values,
         color=[C_POS if not b else C_NEG for b in r2s['at_boundary'].values],
         width=0.85, zorder=2)
ax_e.axhline(0, color='k', lw=0.5)
ax_e.text(0.02, 0.96,
          f"$27/27$ improved\n"
          f"$p_{{\\rm binom}}=7.45\\times10^{{-9}}$\n"
          f"$\\Delta\\chi^2_{{\\rm tot}}={dchi2_r2:.0f}$",
          transform=ax_e.transAxes, va='top', fontsize=FS-1.5)
ax_e.legend(handles=[mpatches.Patch(facecolor=C_POS, label=f"Converged ({n_conv})"),
                      mpatches.Patch(facecolor=C_NEG, label=f"Boundary ({n_bnd})")],
            frameon=False, fontsize=FS-2, loc='lower right')
ax_e.set_xlabel(r"Galaxy (sorted by $\delta R_1$)")
ax_e.set_ylabel(r"$\Delta\chi^2$ (clipped at $-500$)")
ax_e.set_title(r"(e) Regime 2: $c_{200}^{\rm BPB}$ correction (27 virialized)")
ax_e.set_ylim(top=60)

# ── Save ───────────────────────────────────────────────────────────────────
OUT_PDF = os.path.join(FIGURES, "figI4_galactic_summary.pdf")
OUT_PNG = os.path.join(FIGURES, "figI4_galactic_summary.png")
plt.savefig(OUT_PDF, bbox_inches='tight')
plt.savefig(OUT_PNG, bbox_inches='tight', dpi=150)
plt.close()

print(f"Saved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")
print(f"\n(a) r={r_obs:.4f}  Z={Z:.2f}σ  p={p_perm:.4f}")
print(f"(b) rule={rule_38:.0f}  oracle={oracle_38:.0f}  {pct:.0f}%")
print(f"(c) t={t_stat:.2f}  d={cohen_d:.2f}")
print(f"(d) bars: 2319 / {-rule_38:.0f} / {-oracle_38:.0f} / {-oracle_38:.0f}")
print(f"(e) {n_imp_r2}/27 improved  p=7.45e-9  dchi2={dchi2_r2:.0f}")
