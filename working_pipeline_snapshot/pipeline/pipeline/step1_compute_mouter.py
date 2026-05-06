"""
step1_compute_mouter.py
=======================
Compute the outer rotation-curve residual mouter for each galaxy.

Definition (Paper II §4.2):
    mouter = mean(Vobs - V_pred_ref)  for  R > R_half
where:
    R_half     = median radius of all rotation-curve points for the galaxy
    V_pred_ref = sqrt(Vgas² + Ydisk·Vdisk² + Ybul·Vbul² + V_NFW²)
    V_NFW      = NFW halo velocity at Li+2020 best-fit (V200, c200, rs)

Inputs (in results/chains/ relative to repo root):
    MassModels_Lelli2016c.mrt
    sparc_reference_halo_table_nfw_lcdm_clean.csv

Output (in results/chains/):
    mouter_per_galaxy.csv   [name, mouter, n_inner, n_outer, R_half]

Usage (from anywhere in the repo):
    python pipeline/step1_compute_mouter.py
"""

import math
import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths — always relative to repo root (parent of pipeline/)
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

MRT_PATH = os.path.join(CHAINS, "MassModels_Lelli2016c.mrt")
REF_PATH = os.path.join(CHAINS, "sparc_reference_halo_table_nfw_lcdm_clean.csv")
OUT_PATH = os.path.join(CHAINS, "mouter_per_galaxy.csv")

for p in [MRT_PATH, REF_PATH]:
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"\nMissing input: {p}\n"
            f"Place all input files in: {CHAINS}/")

# ---------------------------------------------------------------------------
# Parse MRT
# ---------------------------------------------------------------------------
print(f"Reading {os.path.basename(MRT_PATH)} ...")
mrt_rows = []
with open(MRT_PATH) as f:
    for line in f:
        line = line.rstrip('\r\n')
        if len(line) < 59 or not line[0].isalpha():
            continue
        try:
            mrt_rows.append({
                'name':  line[0:11].strip(),
                'R':     float(line[19:25]),
                'Vobs':  float(line[26:32]),
                'eVobs': float(line[33:38]),
                'Vgas':  float(line[39:45]),
                'Vdisk': float(line[46:52]),
                'Vbul':  float(line[53:59]),
            })
        except:
            continue

mrt = pd.DataFrame(mrt_rows)
print(f"  {len(mrt)} points, {mrt['name'].nunique()} galaxies")

# ---------------------------------------------------------------------------
# Load NFW reference fits
# ---------------------------------------------------------------------------
print(f"Reading {os.path.basename(REF_PATH)} ...")
ref      = pd.read_csv(REF_PATH)
ref_dict = ref.set_index('name').to_dict('index')
print(f"  {len(ref)} galaxies")

# ---------------------------------------------------------------------------
# NFW velocity
# ---------------------------------------------------------------------------
def nfw_v(r_kpc, V200, C200, rs):
    if not all(math.isfinite(x) and x > 0 for x in [r_kpc, V200, C200, rs]):
        return math.nan
    def f(u):
        return math.log(1.0 + u) - u / (1.0 + u)
    r200 = C200 * rs
    x    = r_kpc / r200
    u    = r_kpc / rs
    fc   = f(C200)
    fu   = f(u)
    if x <= 0 or fc <= 0 or fu < 0:
        return math.nan
    v2 = V200 * V200 * fu / (x * fc)
    return math.sqrt(v2) if v2 >= 0 else math.nan

# ---------------------------------------------------------------------------
# Compute mouter for each galaxy
# ---------------------------------------------------------------------------
results = []
for gal in mrt['name'].unique():
    pts = mrt[mrt['name'] == gal].sort_values('R')
    if len(pts) < 4 or gal not in ref_dict:
        continue

    r     = ref_dict[gal]
    V200  = r.get('V200',  math.nan)
    C200  = r.get('C200',  math.nan)
    rs    = r.get('rs',    math.nan)
    Ydisk = r.get('Ydisk', 0.0)
    Ybul  = r.get('Ybul',  0.0)

    if not all(math.isfinite(x) and x > 0 for x in [V200, C200, rs]):
        continue
    if not math.isfinite(Ydisk): Ydisk = 0.0
    if not math.isfinite(Ybul):  Ybul  = 0.0

    R_half = float(pts['R'].median())
    outer  = pts[pts['R'] > R_half]
    inner  = pts[pts['R'] <= R_half]

    if len(outer) == 0:
        continue

    residuals = []
    for _, p in outer.iterrows():
        Vg = p['Vgas']  if math.isfinite(p['Vgas'])  else 0.0
        Vd = p['Vdisk'] if math.isfinite(p['Vdisk']) else 0.0
        Vb = p['Vbul']  if math.isfinite(p['Vbul'])  else 0.0
        vh = nfw_v(p['R'], V200, C200, rs)
        if not math.isfinite(vh):
            continue
        vbar2 = Vg**2 + Ydisk * Vd**2 + Ybul * Vb**2
        Vpred = math.sqrt(max(vbar2 + vh**2, 0.0))
        residuals.append(p['Vobs'] - Vpred)

    if len(residuals) == 0:
        continue

    results.append({
        'name':    gal,
        'mouter':  round(float(np.mean(residuals)), 4),
        'n_inner': len(inner),
        'n_outer': len(residuals),
        'R_half':  round(R_half, 4),
    })

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
df = pd.DataFrame(results)
df.to_csv(OUT_PATH, index=False)

print(f"\nSaved: {OUT_PATH}")
print(f"  {len(df)} galaxies with mouter")
print(f"  mouter range: [{df['mouter'].min():.2f}, {df['mouter'].max():.2f}] km/s")
print(f"  mouter mean:  {df['mouter'].mean():.2f} km/s")
