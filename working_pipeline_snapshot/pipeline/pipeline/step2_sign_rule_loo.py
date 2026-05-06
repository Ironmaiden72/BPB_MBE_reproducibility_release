"""
step2_sign_rule_loo.py
======================
3-condition physical sign rule with strict leave-one-out (LOO)
cross-validation on 38 active SPARC galaxies.

Sign rule (Paper II Eq. 6):
    +1  if  dR1 < T_dR1                                    [strong BPB]
    +1  if  sinner < T_s  AND  mouter * (-dR1) > 0         [moderate]
    -1  otherwise

At each of 38 LOO folds, T_dR1 and T_s are refitted on the 37 training
galaxies by maximising balanced accuracy (BA).

Inputs (in results/chains/):
    per_galaxy_v200_required_amplitude.csv   (Areq, sinner, dR1)
    mouter_per_galaxy.csv                    (output of step1)

Output (in results/chains/):
    loo_predictions.csv   [name, dR1, sinner, mouter, Areq,
                           sign_true, sign_pred, correct,
                           T_dR1_fold, T_s_fold, Areq_pred]

Usage (from anywhere in the repo):
    python pipeline/step2_sign_rule_loo.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, binomtest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

V200_PATH = os.path.join(CHAINS, "per_galaxy_v200_required_amplitude.csv")
MO_PATH   = os.path.join(CHAINS, "mouter_per_galaxy.csv")
OUT_PATH  = os.path.join(CHAINS, "loo_predictions.csv")

for p in [V200_PATH, MO_PATH]:
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"\nMissing input: {p}\n"
            f"Run step1 first, and place input files in: {CHAINS}/")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
v200   = pd.read_csv(V200_PATH)
mo     = pd.read_csv(MO_PATH)

active = v200[v200['active_shape_gate'] == 1].copy()
active = active.merge(mo[['name', 'mouter']], on='name', how='left')

n = len(active)
dR1    = active['delta_R1'].values
sinner = active['s_inner'].values
mouter = active['mouter'].values
Areq   = active['Areq_v_best'].values
st     = np.sign(Areq)

print(f"Active galaxies: {n}")
print(f"  Areq > 0: {(Areq > 0).sum()},  Areq < 0: {(Areq < 0).sum()}")

# ---------------------------------------------------------------------------
# Sign rule and LOO
# ---------------------------------------------------------------------------
T1_GRID = np.linspace(-0.25, -0.05, 41)
T2_GRID = np.linspace(0.20,  0.45,  26)

def apply_rule(dR1, si, mo, T1, T2):
    return np.where(dR1 < T1, 1,
           np.where((si < T2) & (mo * (-dR1) > 0), 1, -1))

def ba(sp, st):
    tp = np.sum((st == 1)  & (sp == 1))
    tn = np.sum((st == -1) & (sp == -1))
    fp = np.sum((st == -1) & (sp == 1))
    fn = np.sum((st == 1)  & (sp == -1))
    s  = tp/(tp+fn) if (tp+fn) > 0 else 0.
    c  = tn/(tn+fp) if (tn+fp) > 0 else 0.
    return 0.5*(s+c)

def best_T(dR1, si, mo, st):
    best_ba, best = -1., (-0.131, 0.320)
    for t1 in T1_GRID:
        for t2 in T2_GRID:
            b = ba(apply_rule(dR1, si, mo, t1, t2), st)
            if b > best_ba:
                best_ba, best = b, (t1, t2)
    return best

print(f"\nRunning LOO ({n} folds) ...")
sign_pred  = np.zeros(n, dtype=int)
T_dR1_fold = np.zeros(n)
T_s_fold   = np.zeros(n)

for i in range(n):
    mask = np.arange(n) != i
    T1, T2 = best_T(dR1[mask], sinner[mask], mouter[mask], st[mask])
    sign_pred[i]  = int(apply_rule(dR1[i:i+1], sinner[i:i+1],
                                   mouter[i:i+1], T1, T2)[0])
    T_dR1_fold[i] = T1
    T_s_fold[i]   = T2

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
correct   = (sign_pred == st)
n_correct = int(correct.sum())

tp = int(((st == 1)  & (sign_pred == 1)).sum())
tn = int(((st == -1) & (sign_pred == -1)).sum())
fp = int(((st == -1) & (sign_pred == 1)).sum())
fn = int(((st == 1)  & (sign_pred == -1)).sum())
BA_loo  = 0.5*(tp/(tp+fn) + tn/(tn+fp))
p_binom = binomtest(n_correct, n, 0.5, alternative='greater').pvalue

Areq_pred = sign_pred * np.abs(Areq)
r_loo, _  = pearsonr(Areq, Areq_pred)

print(f"\n{'='*50}")
print(f"LOO RESULTS")
print(f"{'='*50}")
print(f"  sign correct : {n_correct}/{n}")
print(f"  p_binom      : {p_binom:.5f}")
print(f"  BA_LOO       : {BA_loo:.4f}")
print(f"  r_LOO        : {r_loo:.4f}")
print(f"  T_dR1 mean   : {T_dR1_fold.mean():.4f} ± {T_dR1_fold.std():.4f}")
print(f"  T_s   mean   : {T_s_fold.mean():.4f}  ± {T_s_fold.std():.4f}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
pd.DataFrame({
    'name':       active['name'].values,
    'dR1':        dR1,
    'sinner':     sinner,
    'mouter':     mouter,
    'Areq':       Areq,
    'sign_true':  st.astype(int),
    'sign_pred':  sign_pred,
    'correct':    correct.astype(int),
    'T_dR1_fold': T_dR1_fold,
    'T_s_fold':   T_s_fold,
    'Areq_pred':  Areq_pred,
}).to_csv(OUT_PATH, index=False)

print(f"\nSaved: {OUT_PATH}")
