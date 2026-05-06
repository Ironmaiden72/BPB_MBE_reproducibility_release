"""
step3_null_test.py
==================
Synthetic null test: permute delta_R1 across the 38 active galaxies
500 times and compute the LOO correlation at each permutation.

This tests whether the BPB signal in delta_R1 is specifically responsible
for the sign prediction performance, vs. any structural feature of the
galaxy sample.

Method (Paper II §5.3):
    For each of 500 permutations:
      1. Shuffle dR1 values across galaxies (sinner, mouter, Areq unchanged)
      2. Rerun full LOO with threshold refitting at each fold
      3. Compute r = pearsonr(Areq_true, sign_pred × |Areq_true|)
    The observed r_LOO (from loo_predictions.csv) is compared to this null.

Inputs (in results/chains/):
    loo_predictions.csv    (output of step2)

Output (in results/chains/):
    null_distribution.csv  [r, r_observed]
        r          : 500 permutation LOO correlations
        r_observed : observed LOO correlation (same value in every row)

Usage (from repo root):
    python pipeline/step3_null_test.py
    python pipeline/step3_null_test.py --n_perms 500 --seed 42
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--n_perms", default=500, type=int,
                    help="Number of permutations (default: 500)")
parser.add_argument("--seed",    default=42,  type=int,
                    help="Random seed (default: 42)")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

LOO_PATH = os.path.join(CHAINS, "loo_predictions.csv")
OUT_PATH = os.path.join(CHAINS, "null_distribution.csv")

if not os.path.exists(LOO_PATH):
    raise FileNotFoundError(
        f"\nMissing input: {LOO_PATH}\n"
        f"Run step2 first.")

# ---------------------------------------------------------------------------
# Load LOO predictions
# ---------------------------------------------------------------------------
loo = pd.read_csv(LOO_PATH)

dR1    = loo['dR1'].values.copy()
sinner = loo['sinner'].values
mouter = loo['mouter'].values
Areq   = loo['Areq'].values
st     = loo['sign_true'].values
n      = len(loo)

# Observed r_LOO (from step2, single definition throughout)
r_observed = float(pearsonr(Areq, loo['Areq_pred'].values)[0]) # type: ignore
print(f"r_observed (from step2) = {r_observed:.4f}")

# ---------------------------------------------------------------------------
# Sign rule (identical to step2)
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

def loo_r(dR1_perm):
    """Full LOO correlation with permuted dR1."""
    sp = np.zeros(n, dtype=int)
    for i in range(n):
        mask = np.arange(n) != i
        T1, T2 = best_T(dR1_perm[mask], sinner[mask], mouter[mask], st[mask])
        sp[i] = int(apply_rule(dR1_perm[i:i+1], sinner[i:i+1],
                               mouter[i:i+1], T1, T2)[0])
    return pearsonr(Areq, sp * np.abs(Areq))[0]

# ---------------------------------------------------------------------------
# 500 permutations
# ---------------------------------------------------------------------------
print(f"\nRunning {args.n_perms} permutations (seed={args.seed}) ...")
rng    = np.random.default_rng(args.seed)
null_r = np.zeros(args.n_perms)

for k in range(args.n_perms):
    null_r[k] = loo_r(rng.permutation(dR1)) # type: ignore
    if (k + 1) % 50 == 0:
        print(f"  {k+1:3d}/{args.n_perms}  "
              f"null mean={null_r[:k+1].mean():.3f}  "
              f"std={null_r[:k+1].std():.3f}")

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
null_mean = null_r.mean()
null_std  = null_r.std()
Z         = (r_observed - null_mean) / null_std
p_perm    = (null_r >= r_observed).mean()

print(f"\n{'='*50}")
print(f"NULL TEST RESULTS")
print(f"{'='*50}")
print(f"  null mean    : {null_mean:.4f}")
print(f"  null std     : {null_std:.4f}")
print(f"  r_observed   : {r_observed:.4f}")
print(f"  Z            : {Z:.2f}σ")
print(f"  p_perm       : {p_perm:.4f}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
pd.DataFrame({
    "r":          null_r,
    "r_observed": r_observed,
}).to_csv(OUT_PATH, index=False)

print(f"\nSaved: {OUT_PATH}  ({args.n_perms} rows)")
