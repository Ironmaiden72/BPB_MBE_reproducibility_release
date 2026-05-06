"""
fig3_validation.py  --  Paper II Fig. A3 / 4-panel statistical validation

All data read from results/chains/ — no hardcoded paper values.

Sources:
  null_distribution.csv      <- step3  (LOO Regime 1 null test)
  sparc_active_sample.csv    <- step7  (Areq, dchi2, c200_lcdm_ratio)
  regime2_complete.csv       <- step9  (c200_lcdm_ratio)
  chi2_summary.csv           <- step8  (bar chart values)

Outputs (in results/figures/):
  fig3a_null.pdf, fig3b_pergal.pdf, fig3c_regimes.pdf, fig3d_chi2summary.pdf

Usage (from repo root):
    python pipeline/fig3_validation.py \\
        results/chains/sparc_active_sample.csv \\
        results/chains/regime2_complete.csv
"""

import sys, os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, ttest_ind
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")
OUTDIR     = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(OUTDIR, exist_ok=True)

def req(name):
    p = os.path.join(CHAINS, name)
    if not os.path.exists(p): raise FileNotFoundError(f"Missing: {p}")
    return p

csv_act = sys.argv[1] if len(sys.argv) >= 3 else req("sparc_active_sample.csv")
csv_r2  = sys.argv[2] if len(sys.argv) >= 3 else req("regime2_complete.csv")

df_act  = pd.read_csv(csv_act)
df_r2   = pd.read_csv(csv_r2)
df_null = pd.read_csv(req("null_distribution.csv"))
df_sum  = pd.read_csv(req("chi2_summary.csv"))

# Null
null_r = df_null["r"].dropna().to_numpy()
if "r_observed" in df_null.columns:
    r_obs = float(df_null["r_observed"].iloc[0])
else:
    mask = np.isfinite(df_act["Areq"]) & np.isfinite(df_act["Areq_pred_physical"])
    r_obs, _ = pearsonr(df_act.loc[mask,"Areq"].values,
                        df_act.loc[mask,"Areq_pred_physical"].values)
Z_obs  = (r_obs - null_r.mean()) / null_r.std(ddof=0)
p_perm = (1 + np.sum(null_r >= r_obs)) / (1 + len(null_r))

# Per-galaxy
order   = np.argsort(df_act["Areq"].to_numpy())
dchi2_p = df_act["dchi2_physical"].to_numpy()[order]
dchi2_o = df_act["dchi2_oracle"].to_numpy()[order]
sign_ok = df_act["sign_correct_physical"].to_numpy()[order].astype(bool)
total_p = float(dchi2_p.sum()); n = len(dchi2_p)

# Concentrations
c_act = df_act["c200_lcdm_ratio"].dropna().to_numpy()
c_ina = df_r2["c200_lcdm_ratio"].dropna().to_numpy()
t_val, _ = ttest_ind(c_act, c_ina, equal_var=False)
d_cohen  = ((c_ina.mean()-c_act.mean()) /
            np.sqrt((c_act.std(ddof=0)**2+c_ina.std(ddof=0)**2)/2))

# Bar chart
rules = df_sum["rule"].tolist(); vals = df_sum["minus_dchi2"].tolist()
LABEL_MAP = {"physical_v200":"V200-only\nmodel","full_bpb":"Full BPB\nmodel"}
xlbls = [LABEL_MAP.get(r,r) for r in rules]

C_POS="#4C72B0"; C_NEG="#DD8452"; C_GRN="#55A868"; C_GRY="#8E8E8E"; FS=9
plt.rcParams.update({"text.usetex":False,"mathtext.fontset":"dejavusans","font.size":9})
rng = np.random.default_rng(0)

# (a) Null
fig,ax = plt.subplots(figsize=(3.2,2.6))
ax.hist(null_r,bins=25,density=True,color=C_GRY,alpha=0.7,label=r"Null ($\delta R_1$ shuffled)")
ax.axvline(r_obs,color=C_POS,lw=1.8,label=f"Observed ($r={r_obs:.3f}$)")
ax.text(0.97,0.06,f"$p={p_perm:.3f}$\n$Z={Z_obs:.2f}\\sigma$",
        fontsize=7,color=C_POS,transform=ax.transAxes,ha="right",va="bottom")
ax.set_xlabel("LOO correlation $r$",fontsize=FS); ax.set_ylabel("Density",fontsize=FS)
ax.legend(fontsize=6.5,loc="upper left"); ax.tick_params(labelsize=7)
fig.tight_layout()
out=os.path.join(OUTDIR,"fig3a_null.pdf"); fig.savefig(out,dpi=300,bbox_inches="tight"); plt.close(fig)
print(f"Saved: {out}")

# (b) Per-galaxy
CLIP=-1500; cols=[C_POS if ok else C_NEG for ok in sign_ok]
fig,ax = plt.subplots(figsize=(3.8,2.8))
ax.step(np.arange(n),np.clip(dchi2_o,CLIP,None),where="mid",color="k",ls="--",lw=1.0,label="Oracle")
for i,(y,c) in enumerate(zip(np.clip(dchi2_p,CLIP,None),cols)):
    ax.bar(i,y,width=0.7,color=c,alpha=0.8,zorder=3)
ax.axhline(0,color="k",lw=0.6)
ax.set_xlabel(r"Galaxy (sorted by $A_{\rm req}$)",fontsize=FS)
ax.set_ylabel(r"$\Delta\chi^2$",fontsize=FS)
ax.text(0.05,0.05,f"Total $\\Delta\\chi^2={total_p:.0f}$",transform=ax.transAxes,fontsize=7)
leg=[Patch(facecolor=C_POS,label=f"Improved ({sign_ok.sum()}/{n})"),
     Patch(facecolor=C_NEG,label="Worsened"),
     plt.Line2D([0],[0],color="k",ls="--",label="Oracle")]
ax.legend(handles=leg,fontsize=6.5,loc="lower right"); ax.tick_params(labelsize=7)
fig.tight_layout()
out=os.path.join(OUTDIR,"fig3b_pergal.pdf"); fig.savefig(out,dpi=300,bbox_inches="tight"); plt.close(fig)
print(f"Saved: {out}")

# (c) Violin
groups=[c_act,c_ina]; labels=["Actives\n($s<0.45$)",r"Inact. strong"]
cols_v=[C_POS,C_GRN]
fig,ax = plt.subplots(figsize=(3.0,2.9))
parts=ax.violinplot(groups,positions=[1,2],showmedians=True,showextrema=False,widths=0.55)
for pc,col in zip(parts["bodies"],cols_v): pc.set_facecolor(col); pc.set_alpha(0.6)
parts["cmedians"].set_color("k"); parts["cmedians"].set_linewidth(1.5)
for i,(g,col) in enumerate(zip(groups,cols_v)):
    ax.scatter(rng.normal(i+1,0.06,len(g)),g,c=col,s=6,alpha=0.5,zorder=3)
ax.axhline(1.0,color="k",ls="--",lw=1.0)
ax.text(0.5,0.97,f"|t| = {abs(t_val):.2f},  d = {abs(d_cohen):.2f}",
        transform=ax.transAxes,ha="center",va="top",fontsize=6.5,
        bbox=dict(boxstyle="round,pad=0.2",fc="white",alpha=0.8))
ax.set_xticks([1,2]); ax.set_xticklabels(labels,fontsize=7.5)
ax.set_ylabel(r"$c_{200}/c_{200}^{\Lambda{\rm CDM}}$",fontsize=FS); ax.tick_params(labelsize=7)
fig.tight_layout()
out=os.path.join(OUTDIR,"fig3c_regimes.pdf"); fig.savefig(out,dpi=300,bbox_inches="tight"); plt.close(fig)
print(f"Saved: {out}")

# (d) Bar
bcols=[C_GRN if "full" in r.lower() else C_POS for r in rules]
fig,ax = plt.subplots(figsize=(3.0,2.7))
bars=ax.bar(range(len(rules)),vals,color=bcols,width=0.5,edgecolor="k",linewidth=0.5)
for bar,val in zip(bars,vals):
    ax.text(bar.get_x()+bar.get_width()/2,val+max(vals)*0.018,
            f"$-{val:,}$".replace(",","\\,"),ha="center",va="bottom",fontsize=7,fontweight="bold")
ax.set_xticks(range(len(rules))); ax.set_xticklabels(xlbls,fontsize=7.5)
ax.set_ylabel(r"$-\Delta\chi^2$",fontsize=FS); ax.set_ylim(0,max(vals)*1.20); ax.tick_params(labelsize=7)
fig.tight_layout()
out=os.path.join(OUTDIR,"fig3d_chi2summary.pdf"); fig.savefig(out,dpi=300,bbox_inches="tight"); plt.close(fig)
print(f"Saved: {out}")

print(f"\n  (a) r={r_obs:.3f}, Z={Z_obs:.2f}σ, p={p_perm:.3f}")
print(f"  (b) total={total_p:.0f}, {sign_ok.sum()}/{n} improved")
print(f"  (c) |t|={abs(t_val):.2f}, d={abs(d_cohen):.2f}")
print(f"  (d) bars: {vals}")
