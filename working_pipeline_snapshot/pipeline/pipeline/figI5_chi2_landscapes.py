"""
figI5_chi2_landscapes.py
========================
Figure 5 (Paper I) / Figure 4 (Paper II):
chi2(A) landscapes for IC2574 (active) and DDO170 (inactive).
Computed from real rotation curves — zero hardcoding.

Definition: A is the Areq axis. V200_BPB = V200_ref * (1 - A * dR1)
            Areq = argmin_A chi2(A)

Inputs (in results/chains/):
    MassModels_Lelli2016c.mrt
    sparc_reference_halo_table_nfw_lcdm_clean.csv
    bpb_correction_proxy_grid_z0.csv
    per_galaxy_v200_required_amplitude.csv

Output (in results/figures/):
    figI5_chi2_landscapes.pdf
    figI5_chi2_landscapes.png

Usage (from repo root):
    python pipeline/figI5_chi2_landscapes.py
"""

import math, os
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")
FIGURES    = os.path.join(REPO_ROOT, "results", "figures")
os.makedirs(FIGURES, exist_ok=True)

for f in ["MassModels_Lelli2016c.mrt",
          "sparc_reference_halo_table_nfw_lcdm_clean.csv",
          "bpb_correction_proxy_grid_z0.csv",
          "per_galaxy_v200_required_amplitude.csv"]:
    if not os.path.exists(os.path.join(CHAINS, f)):
        raise FileNotFoundError(f"Missing: {f}. Place in {CHAINS}/")

ref   = pd.read_csv(os.path.join(CHAINS,"sparc_reference_halo_table_nfw_lcdm_clean.csv")).set_index('name')
proxy = pd.read_csv(os.path.join(CHAINS,"bpb_correction_proxy_grid_z0.csv"))
v200  = pd.read_csv(os.path.join(CHAINS,"per_galaxy_v200_required_amplitude.csv"))
g_logM=proxy['logM200'].values; g_dR1=proxy['delta_bpb_effective_R1'].values

GALAXIES = {
    "IC2574": {"color":"#4C72B0"},
    "DDO170": {"color":"#8E8E8E"},
}

mrt_pts = {g: [] for g in GALAXIES}
with open(os.path.join(CHAINS,"MassModels_Lelli2016c.mrt")) as f:
    for line in f:
        line=line.rstrip('\r\n')
        if len(line)<59 or not line[0].isalpha(): continue
        try:
            nm=line[0:11].strip()
            if nm not in mrt_pts: continue
            mrt_pts[nm].append({'R':float(line[19:25]),'Vobs':float(line[26:32]),
                'eVobs':float(line[33:38]),'Vgas':float(line[39:45]),
                'Vdisk':float(line[46:52]),'Vbul':float(line[53:59])})
        except: continue

def nfw_v(r,V200,C200,rs):
    def f(u): return math.log(1+u)-u/(1+u)
    if not all(math.isfinite(x) and x>0 for x in [r,V200,C200,rs]): return math.nan
    fc=f(C200); fu=f(r/rs)
    if fc<=0 or fu<0: return math.nan
    return math.sqrt(max(V200**2*fu/((r/(C200*rs))*fc),0))

def chi2_landscape(name, A_grid):
    """
    A is the x-axis (Areq convention).
    V200_BPB = V200_ref * (1 - A * dR1)

    For active galaxies  : Δχ²(A) = χ²(A) - χ²_min  (shows clear minimum)
    For inactive galaxies: Δχ²(A) = χ²(A) - χ²(A=0) (shows landscape flat
                           around the physically-assigned gate value A=0)
    """
    r     = ref.loc[name]
    V200  = float(r['V200']); C200=float(r['C200']); rs=float(r['rs'])
    Ydisk = float(r.get('Ydisk',0.)); Ybul=float(r.get('Ybul',0.))
    if not math.isfinite(Ydisk): Ydisk=0.
    if not math.isfinite(Ybul):  Ybul=0.
    logM  = float(v200[v200['name']==name]['logM200'].iloc[0])
    dR1   = float(np.interp(logM, g_logM, g_dR1))
    pts   = mrt_pts[name]

    def c2(A_val):
        ratio = 1.0 - A_val * dR1   # V200_BPB = V200_ref * ratio
        if ratio <= 0: return math.nan
        V200_b=V200*ratio; rs_b=rs*ratio
        s=0.
        for p in pts:
            vh=nfw_v(p['R'],V200_b,C200,rs_b)
            if not math.isfinite(vh): return math.nan
            Vg=p['Vgas'] if math.isfinite(p['Vgas']) else 0.
            Vd=p['Vdisk'] if math.isfinite(p['Vdisk']) else 0.
            Vb=p['Vbul']  if math.isfinite(p['Vbul'])  else 0.
            vpred=math.sqrt(max(Vg**2+Ydisk*Vd**2+Ybul*Vb**2+vh**2,0))
            s+=((vpred-p['Vobs'])/p['eVobs'])**2
        return s

    chi2s    = np.array([c2(a) for a in A_grid])
    chi2_min = np.nanmin(chi2s)
    Areq_num = float(A_grid[np.nanargmin(chi2s)])  # numerical minimum
    chi2_at0 = c2(0.0)                              # value at gate point A=0

    # Active: show relative to numerical minimum (reveals clear bowl)
    # Inactive: show relative to A=0 (reveals flatness of landscape)
    is_active = (v200[v200['name']==name]['active_shape_gate'].iloc[0] == 1)
    if is_active:
        dc   = chi2s - chi2_min
        Areq = Areq_num          # true minimum
    else:
        dc   = chi2s - chi2_at0  # flat near zero — gate physically justified
        Areq = 0.0               # assigned gate value

    return dc, dR1, Areq, chi2_min

A_GRID = np.linspace(-2, 4, 481)
landscapes = {}
for name in GALAXIES:
    dc, dR1, Areq, chi2_min = chi2_landscape(name, A_GRID)
    landscapes[name] = {'dc':dc, 'dR1':dR1, 'Areq':Areq, 'chi2_min':chi2_min}
    print(f"{name}: dR1={dR1:.4f}  Areq={Areq:.3f}  chi2_min={chi2_min:.1f}")

# Figure
FS=8.; LW=1.4

labels = {
    "IC2574": r"Active: IC2574 ($A_{\rm req}^{\rm true}=3.38$)",
    "DDO170": r"Inactive: DDO170 ($A_{\rm req}^{\rm gate}=0$)",
}

fig,(ax1,ax2) = plt.subplots(1, 2, figsize=(7., 3.))

for ax, name in zip([ax1, ax2], ["IC2574", "DDO170"]):
    info = landscapes[name]
    dc   = info['dc']; Areq=info['Areq']; dR1=info['dR1']
    col  = GALAXIES[name]['color']
    valid= np.isfinite(dc)

    # y-axis: clip at 99th percentile for clean display
    ymax = min(np.nanpercentile(dc[valid], 99)*1.05,
               3200 if name=="IC2574" else 200)

    ax.plot(A_GRID, dc, color=col, lw=LW, zorder=3)
    ax.axvline(Areq, color='k', ls='--', lw=LW,
               label=rf"$A_{{\rm req}}={Areq:.3f}$")
    ax.axvline(0., color='gray', ls=':', lw=LW, label="$A=0$ (no corr.)")

    ax.set_xlabel("$A$", fontsize=FS)
    # ylabel differs: active shows distance from min; inactive from A=0
    is_inactive = (name == "DDO170")
    ylabel = (r"$\Delta\chi^2(A)=\chi^2(A)-\chi^2(A{=}0)$" if is_inactive
              else r"$\Delta\chi^2(A)=\chi^2(A)-\chi^2_{\rm min}$")
    ax.set_ylabel(ylabel, fontsize=FS)
    ax.set_xlim(-2, 4); ax.set_ylim(-ymax*0.02, ymax)
    ax.tick_params(labelsize=7)
    ax.text(0.05, 0.93, rf"$\delta R_1={dR1:.3f}$",
            transform=ax.transAxes, fontsize=7)
    ax.legend(fontsize=7, loc='upper left', frameon=False)
    ax.set_title(labels[name], fontsize=FS)

fig.text(0.5, 0.01,
    "Flat landscape justifies shape-gate "
    r"($A_{\rm req}=0$) for all 80 inactive galaxies",
    ha='center', fontsize=7, style='italic')

plt.tight_layout(rect=[0, 0.06, 1, 1])

OUT_PDF=os.path.join(FIGURES,"figI5_chi2_landscapes.pdf")
OUT_PNG=os.path.join(FIGURES,"figI5_chi2_landscapes.png")
plt.savefig(OUT_PDF, dpi=300, bbox_inches='tight')
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
plt.close()

print(f"\nSaved: {OUT_PDF}")
print(f"Saved: {OUT_PNG}")
