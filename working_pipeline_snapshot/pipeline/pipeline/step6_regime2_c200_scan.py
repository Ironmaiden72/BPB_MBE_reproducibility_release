"""
step6_regime2_c200_scan.py
==========================
Quantitative demonstration of Regime 2 (Paper II §7):
BPB signal projects into c200 for virialized inactive galaxies.

Protocol (Paper II §7.3):
    Population: inactive (sinner > 0.45) strong-signal (dR1 < -0.131) galaxies,
                excluding UGC03580 (extreme chi2_ref = 5086, NFW failure).
    For each galaxy:
        1. Fix V200 = V200_ref  (locked by flat rotation curve)
        2. For each alpha_c in [-20, +20]:
               c200_BPB = c200_ref * (1 + alpha_c * dR1)
               rs_BPB   = r200_ref / c200_BPB          (r200 unchanged)
           3. Minimise chi2 over Ydisk in [0.05, 3.0]
           4. Record alpha_c_opt, delta_chi2 = chi2(c200_BPB) - chi2(c200_ref)

    Statistics:
        - n_improved / 27  (pbinom)
        - n_positive_alpha_c / n_converged  (psign)
        - pearsonr(dR1, alpha_c_opt)  for converged galaxies

Inputs (in results/chains/):
    MassModels_Lelli2016c.mrt
    sparc_reference_halo_table_nfw_lcdm_clean.csv
    bpb_correction_proxy_grid_z0.csv
    per_galaxy_v200_required_amplitude.csv     (provides sinner, dR1, active_gate)

Output (in results/chains/):
    regime2_c200_scan.csv   [name, dR1, c200_ref, alpha_c_opt, delta_chi2,
                              at_boundary, converged]

Usage (from repo root):
    python pipeline/step6_regime2_c200_scan.py
"""

import math
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import pearsonr, binomtest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

MRT_PATH   = os.path.join(CHAINS, "MassModels_Lelli2016c.mrt")
REF_PATH   = os.path.join(CHAINS, "sparc_reference_halo_table_nfw_lcdm_clean.csv")
PROXY_PATH = os.path.join(CHAINS, "bpb_correction_proxy_grid_z0.csv")
V200_PATH  = os.path.join(CHAINS, "per_galaxy_v200_required_amplitude.csv")
OUT_PATH   = os.path.join(CHAINS, "regime2_c200_scan.csv")

for p in [MRT_PATH, REF_PATH, PROXY_PATH, V200_PATH]:
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"\nMissing: {p}\n"
            f"Run steps 1-2 first and place inputs in {CHAINS}/")

# ---------------------------------------------------------------------------
# Protocol parameters (Paper II §7.3)
# ---------------------------------------------------------------------------
DR1_THRESHOLD    = -0.131   # strong-signal gate
SINNER_THRESHOLD =  0.45    # inactive gate (sinner > this)
EXCLUDE          = {"UGC03580"}  # extreme chi2_ref outlier
ALPHA_C_GRID     = np.linspace(-20, 20, 801)
YDISK_BOUNDS     = (0.05, 3.0)
BOUNDARY_TOL     = 0.5      # |alpha_c| > 20 - tol → at boundary

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print("Loading inputs ...")
ref   = pd.read_csv(REF_PATH).set_index('name')
proxy = pd.read_csv(PROXY_PATH)
v200  = pd.read_csv(V200_PATH)

g_logM = proxy['logM200'].values
g_dR1  = proxy['delta_bpb_effective_R1'].values

# Select inactive strong-signal population
inactive_strong = v200[
    (v200['active_shape_gate'] == 0) &
    (v200['delta_R1'] < DR1_THRESHOLD) &
    (~v200['name'].isin(EXCLUDE))
].copy()

print(f"Inactive strong-signal galaxies: {len(inactive_strong)} "
      f"(excl. {EXCLUDE})")

# ---------------------------------------------------------------------------
# Parse MRT for this population only
# ---------------------------------------------------------------------------
target_names = set(inactive_strong['name'])
mrt_pts = {}
print(f"Parsing MRT for {len(target_names)} galaxies ...")

with open(MRT_PATH) as f:
    for line in f:
        line = line.rstrip('\r\n')
        if len(line) < 59 or not line[0].isalpha():
            continue
        try:
            name = line[0:11].strip()
            if name not in target_names:
                continue
            if name not in mrt_pts:
                mrt_pts[name] = []
            mrt_pts[name].append({
                'R':     float(line[19:25]),
                'Vobs':  float(line[26:32]),
                'eVobs': float(line[33:38]),
                'Vgas':  float(line[39:45]),
                'Vdisk': float(line[46:52]),
                'Vbul':  float(line[53:59]),
            })
        except:
            continue

print(f"  Rotation curves loaded: {len(mrt_pts)} galaxies")

# ---------------------------------------------------------------------------
# NFW velocity
# ---------------------------------------------------------------------------
def nfw_v(r, V200, C200, rs):
    """NFW circular velocity. rs = r200/C200."""
    def f(u): return math.log(1.0 + u) - u / (1.0 + u)
    if not all(math.isfinite(x) and x > 0 for x in [r, V200, C200, rs]):
        return math.nan
    r200 = C200 * rs
    x    = r / r200
    u    = r / rs
    fc   = f(C200)
    fu   = f(u)
    if x <= 0 or fc <= 0 or fu < 0:
        return math.nan
    return math.sqrt(max(V200**2 * fu / (x * fc), 0.0))

# ---------------------------------------------------------------------------
# Chi2 at fixed (V200, c200_BPB) with Ydisk optimised analytically
# ---------------------------------------------------------------------------
def chi2_at_c200(name, alpha_c, V200, C200_ref, r200_ref,
                 Ybul, dR1, pts):
    """
    c200_BPB = C200_ref * (1 + alpha_c * dR1)
    rs_BPB   = r200_ref / c200_BPB   (r200 fixed, V200 fixed)
    Minimise chi2 over Ydisk in [0.05, 3.0]
    """
    c200_bpb = C200_ref * (1.0 + alpha_c * dR1)
    if c200_bpb <= 0 or not math.isfinite(c200_bpb):
        return math.nan, math.nan

    rs_bpb = r200_ref / c200_bpb
    if rs_bpb <= 0 or not math.isfinite(rs_bpb):
        return math.nan, math.nan

    def chi2_ydisk(Ydisk):
        c2 = 0.0
        for p in pts:
            vh = nfw_v(p['R'], V200, c200_bpb, rs_bpb)
            if not math.isfinite(vh):
                return 1e10
            Vg = p['Vgas']  if math.isfinite(p['Vgas'])  else 0.0
            Vd = p['Vdisk'] if math.isfinite(p['Vdisk']) else 0.0
            Vb = p['Vbul']  if math.isfinite(p['Vbul'])  else 0.0
            vbar2 = Vg**2 + Ydisk * Vd**2 + Ybul * Vb**2
            vpred = math.sqrt(max(vbar2 + vh**2, 0.0))
            c2 += ((vpred - p['Vobs']) / p['eVobs'])**2
        return c2

    res = minimize_scalar(chi2_ydisk, bounds=YDISK_BOUNDS, method='bounded')
    return res.fun, res.x

# ---------------------------------------------------------------------------
# Reference chi2 (alpha_c = 0, Ydisk from Li+2020)
# ---------------------------------------------------------------------------
def chi2_ref_fixed(name, V200, C200, rs, Ydisk, Ybul, pts):
    c2 = 0.0
    for p in pts:
        vh = nfw_v(p['R'], V200, C200, rs)
        if not math.isfinite(vh):
            return math.nan
        Vg = p['Vgas']  if math.isfinite(p['Vgas'])  else 0.0
        Vd = p['Vdisk'] if math.isfinite(p['Vdisk']) else 0.0
        Vb = p['Vbul']  if math.isfinite(p['Vbul'])  else 0.0
        vbar2 = Vg**2 + Ydisk * Vd**2 + Ybul * Vb**2
        vpred = math.sqrt(max(vbar2 + vh**2, 0.0))
        c2 += ((vpred - p['Vobs']) / p['eVobs'])**2
    return c2

# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------
print(f"\nScanning alpha_c over [{ALPHA_C_GRID[0]:.0f}, {ALPHA_C_GRID[-1]:.0f}] "
      f"({len(ALPHA_C_GRID)} points) for {len(inactive_strong)} galaxies ...")

results = []

for _, row in inactive_strong.iterrows():
    name = row['name']
    logM = float(row['logM200'])
    dR1  = float(np.interp(logM, g_logM, g_dR1))

    if name not in ref.index or name not in mrt_pts:
        print(f"  SKIP {name}: missing ref or rotation curve")
        continue

    r    = ref.loc[name]
    V200  = float(r['V200'])
    C200  = float(r['C200'])
    rs    = float(r['rs'])
    Ydisk = float(r.get('Ydisk', 0.0))
    Ybul  = float(r.get('Ybul',  0.0))

    if not all(math.isfinite(x) and x > 0 for x in [V200, C200, rs]):
        print(f"  SKIP {name}: invalid NFW params")
        continue
    if not math.isfinite(Ydisk): Ydisk = 0.0
    if not math.isfinite(Ybul):  Ybul  = 0.0

    r200_ref = C200 * rs  # fixed
    pts      = mrt_pts[name]

    # Reference chi2 at alpha_c = 0
    c2_ref = chi2_ref_fixed(name, V200, C200, rs, Ydisk, Ybul, pts)
    if not math.isfinite(c2_ref):
        print(f"  SKIP {name}: non-finite chi2_ref")
        continue

    # Scan alpha_c
    chi2_scan = np.full(len(ALPHA_C_GRID), np.nan)
    for j, ac in enumerate(ALPHA_C_GRID):
        c2, _ = chi2_at_c200(name, ac, V200, C200, r200_ref, Ybul, dR1, pts)
        chi2_scan[j] = c2

    valid = np.isfinite(chi2_scan)
    if not valid.any():
        print(f"  SKIP {name}: no valid alpha_c")
        continue

    best_j    = int(np.nanargmin(chi2_scan))
    ac_opt    = float(ALPHA_C_GRID[best_j])
    c2_best   = float(chi2_scan[best_j])
    dchi2     = c2_best - c2_ref
    at_bnd    = int(abs(ac_opt) >= (20.0 - BOUNDARY_TOL))
    converged = int(not at_bnd)

    results.append({
        'name':        name,
        'dR1':         round(dR1, 6),
        'c200_ref':    round(C200, 4),
        'chi2_ref':    round(c2_ref, 2),
        'alpha_c_opt': round(ac_opt, 4),
        'chi2_best':   round(c2_best, 2),
        'delta_chi2':  round(dchi2, 2),
        'at_boundary': at_bnd,
        'converged':   converged,
    })
    print(f"  {name:12s}  dR1={dR1:.4f}  alpha_c={ac_opt:6.2f}  "
          f"dchi2={dchi2:8.1f}  {'BOUNDARY' if at_bnd else 'converged'}")

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
df = pd.DataFrame(results)
df.to_csv(OUT_PATH, index=False)

n_total   = len(df)
n_imp     = int((df['delta_chi2'] < 0).sum())
conv      = df[df['converged'] == 1]
n_conv    = len(conv)
n_pos_ac  = int((conv['alpha_c_opt'] > 0).sum())

p_binom   = binomtest(n_imp, n_total, 0.5, alternative='greater').pvalue
p_sign    = binomtest(n_pos_ac, n_conv, 0.5, alternative='greater').pvalue

r_corr, p_corr = pearsonr(conv['dR1'].values, conv['alpha_c_opt'].values)

dchi2_tot = df['delta_chi2'].sum()

print(f"\n{'='*55}")
print(f"REGIME 2 RESULTS")
print(f"{'='*55}")
print(f"  n galaxies     : {n_total}")
print(f"  n improved     : {n_imp}/{n_total}  (pbinom={p_binom:.2e})")
print(f"  n converged    : {n_conv}")
print(f"  n alpha_c > 0  : {n_pos_ac}/{n_conv}  (psign={p_sign:.2e})")
print(f"  r(dR1, alpha_c): {r_corr:.4f}  (p={p_corr:.2e})")
print(f"  total dchi2    : {dchi2_tot:.0f}")
print(f"\n  Paper II targets:")
print(f"    27/27 improved, pbinom=1.5e-8")
print(f"    25/25 alpha_c>0, psign≈0")
print(f"    r=0.981, p≈0")
print(f"    dchi2=-21528")
print(f"\nSaved: {OUT_PATH}")
