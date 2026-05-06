"""
step4b_type95_baseline.py
==========================
Compute the Type-9.5 morphological baseline for the 38 active galaxies.

Sign rule: predict +1 (positive V200 correction) if galaxy Type <= 9.5,
           predict -1 otherwise. This is the morphological baseline
           used in Paper I/II Table 1 (Type-9.5 row).

For each active galaxy:
    sign_pred_type95 = +1 if Type <= 9.5, else -1
    alpha_type95 = -sign_pred_type95 * |Areq|  (same amplitude as physical rule)
    delta_chi2_type95(i) = chi2(alpha_type95) - chi2_ref

Inputs (in results/chains/):
    MassModels_Lelli2016c.mrt
    sparc_reference_halo_table_nfw_lcdm_clean.csv
    bpb_correction_proxy_grid_z0.csv
    per_galaxy_v200_required_amplitude.csv     (Type column used here)
    loo_predictions.csv                         (output of step2)

Output (in results/chains/):
    type95_per_galaxy.csv   [name, Type, sign_type95, sign_true,
                              correct_type95, dchi2_type95]

Usage (from repo root):
    python pipeline/step4b_type95_baseline.py
"""

import math
import os
import numpy as np
import pandas as pd

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
LOO_PATH   = os.path.join(CHAINS, "loo_predictions.csv")
OUT_PATH   = os.path.join(CHAINS, "type95_per_galaxy.csv")

for p in [MRT_PATH, REF_PATH, PROXY_PATH, V200_PATH, LOO_PATH]:
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"\nMissing: {p}\n"
            f"Run steps 1-2 first and place all inputs in {CHAINS}/")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print("Loading inputs ...")
ref   = pd.read_csv(REF_PATH).set_index('name')
proxy = pd.read_csv(PROXY_PATH)
v200  = pd.read_csv(V200_PATH)
loo   = pd.read_csv(LOO_PATH)
type_map = v200.set_index('name')['Type'].to_dict()
loo['Type'] = loo['name'].map(type_map)
loo['sign_type95'] = loo['Type'].apply(lambda t: 1 if t <= 9.5 else -1)

# dR1 interpolation grid
g_logM = proxy['logM200'].values
g_dR1  = proxy['delta_bpb_effective_R1'].values

# logM200 per galaxy from v200
logM_map = dict(zip(v200['name'], v200['logM200']))

# ---------------------------------------------------------------------------
# Parse MRT — load only active galaxies
# ---------------------------------------------------------------------------
active_names = set(loo['name'])
mrt_pts = {}
print(f"Parsing MRT for {len(active_names)} active galaxies ...")
with open(MRT_PATH) as f:
    for line in f:
        line = line.rstrip('\r\n')
        if len(line) < 59 or not line[0].isalpha():
            continue
        try:
            name = line[0:11].strip()
            if name not in active_names:
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

print(f"  Loaded rotation curves for {len(mrt_pts)} galaxies")

# ---------------------------------------------------------------------------
# NFW velocity
# ---------------------------------------------------------------------------
def nfw_v(r, V200, C200, rs):
    def f(u): return math.log(1.0 + u) - u / (1.0 + u)
    if not all(math.isfinite(x) and x > 0 for x in [r, V200, C200, rs]):
        return math.nan
    x  = r / (C200 * rs)
    u  = r / rs
    fc = f(C200)
    fu = f(u)
    if x <= 0 or fc <= 0 or fu < 0:
        return math.nan
    return math.sqrt(max(V200**2 * fu / (x * fc), 0.0))

# ---------------------------------------------------------------------------
# Chi2 at given alpha
# ---------------------------------------------------------------------------
def chi2_at(name, alpha, V200, C200, rs, Ydisk, Ybul, dR1):
    ratio = 1.0 + alpha * dR1
    if ratio <= 0:
        return 1e10
    V200_b = V200 * ratio
    rs_b   = rs   * ratio          # c200 unchanged, rs scales with V200
    c2 = 0.0
    for p in mrt_pts.get(name, []):
        vh = nfw_v(p['R'], V200_b, C200, rs_b)
        if not math.isfinite(vh):
            return 1e10
        vbar2 = (p['Vgas']**2
                 + Ydisk * p['Vdisk']**2
                 + Ybul  * p['Vbul']**2)
        vpred = math.sqrt(max(vbar2 + vh**2, 0.0))
        c2   += ((vpred - p['Vobs']) / p['eVobs'])**2
    return c2

# ---------------------------------------------------------------------------
# Compute delta chi2 for each active galaxy
# ---------------------------------------------------------------------------
print("Computing delta chi2 ...")
results = []

for _, row in loo.iterrows():
    name       = row['name']
    sign_pred  = int(row['sign_type95'])   # type-9.5 morphological rule
    sign_true  = int(row['sign_true'])
    correct    = bool(sign_pred == sign_true)
    Areq       = float(row['Areq'])

    galtype    = float(row['Type'])

    if name not in ref.index or name not in logM_map:
        print(f"  SKIP {name}: not in ref or logM map")
        continue

    r      = ref.loc[name]
    V200   = float(r['V200'])
    C200   = float(r['C200'])
    rs     = float(r['rs'])
    Ydisk  = float(r.get('Ydisk', 0.0))
    Ybul   = float(r.get('Ybul',  0.0))

    if not all(math.isfinite(x) and x > 0 for x in [V200, C200, rs]):
        print(f"  SKIP {name}: invalid NFW params")
        continue

    if not math.isfinite(Ydisk): Ydisk = 0.0
    if not math.isfinite(Ybul):  Ybul  = 0.0

    logM = logM_map[name]
    dR1  = float(np.interp(logM, g_logM, g_dR1)) # type: ignore

    # Amplitude = |Areq_scan| — single definition
    amp = abs(Areq)

    # alpha: gate wrong-sign predictions at 0
    alpha_rule = -(sign_pred * amp) if correct else 0.0   # gate wrong preds at 0

    c2_ref    = chi2_at(name, 0.0,        V200, C200, rs, Ydisk, Ybul, dR1)
    c2_rule   = chi2_at(name, alpha_rule, V200, C200, rs, Ydisk, Ybul, dR1)

    results.append({
        'name':           name,
        'Type':           galtype,
        'sign_type95':    sign_pred,
        'sign_true':      sign_true,
        'correct_type95': int(correct),
        'dchi2_type95':   round(c2_rule - c2_ref, 2),
    })

# ---------------------------------------------------------------------------
# Summary and save
# ---------------------------------------------------------------------------
df = pd.DataFrame(results)
df.to_csv(OUT_PATH, index=False)

total_type95 = df['dchi2_type95'].sum()
n_correct    = int(df['correct_type95'].sum())

print(f"\n{'='*50}")
print(f"TYPE-9.5 BASELINE RESULTS")
print(f"{'='*50}")
print(f"  n galaxies       : {len(df)}")
print(f"  sign correct     : {n_correct}/38")
print(f"  Δχ²_type95 (38)  : {total_type95:.0f}")
print(f"\n  Paper I/II reference value: -2319")
print(f"  Computed value:             {total_type95:.0f}")
print(f"\nSaved: {OUT_PATH}")
