"""
build_null_distribution.py
--------------------------
Reconstruct null_distribution.csv by permuting delta_R1 across the 38 active
galaxies 500 times and computing the LOO correlation at each permutation.

This reproduces Paper II §5.3 / Paper I §5.3 synthetic null test.

Method
------
For each of 500 permutations:
  1. Shuffle dR1 values across galaxies (preserving all other columns intact).
  2. Re-apply the 3-condition physical sign rule with LOO threshold refitting
     (T_dR1, T_s re-optimised on 37 training galaxies at each fold).
  3. Compute the predicted Areq as:  sign_pred × |Areq_true|
     (amplitude bottleneck is the sign; using |Areq_true| as amplitude proxy
     gives the upper-bound correlation under that sign — identical to the
     procedure used in the paper's amplitude model).
  4. Record the Pearson r between Areq_pred and Areq_true.

The observed LOO correlation (r_observed) is also computed from the real data
(non-permuted) and written as a separate column.

Inputs
------
  sparc_active_sample.csv   (columns: dR1, sinner, Areq, mouter)

Outputs
-------
  null_distribution.csv     (columns: r, r_observed)
    - 500 rows, one per permutation
    - r_observed is the same value in every row (for easy reading by figI4)

Usage
-----
  python build_null_distribution.py
  python build_null_distribution.py --data /path/to/data/ --n_perms 500 --seed 42
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--data",    default=None, help="Path to data directory")
parser.add_argument("--n_perms", default=500,  type=int)
parser.add_argument("--seed",    default=42,   type=int)
args = parser.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))

def _find_data_dir(override=None):
    if override:
        return os.path.abspath(override)
    for d in [os.path.join(HERE, "data"),
              os.path.join(HERE, "..", "data")]:
        if os.path.isdir(d):
            return os.path.abspath(d)
    raise FileNotFoundError("Data directory not found. Use --data.")

DATA_DIR = _find_data_dir(args.data)

df = pd.read_csv(os.path.join(DATA_DIR, "sparc_active_sample.csv"))
for col in ("dR1", "sinner", "Areq", "mouter"):
    if col not in df.columns:
        raise ValueError(f"sparc_active_sample.csv missing column: {col}")

dR1_real  = df["dR1"].values.copy()
sinner    = df["sinner"].values
Areq_true = df["Areq"].values
mouter    = df["mouter"].values
n         = len(df)

# ---------------------------------------------------------------------------
# Sign rule — identical to Paper II Eq.(6)
# ---------------------------------------------------------------------------
T_DR1_SCAN = np.linspace(-0.25, -0.05, 41)   # search grid for T_dR1
T_S_SCAN   = np.linspace(0.20,  0.45,  26)   # search grid for T_s

def apply_rule(dR1_row, sinner_row, mouter_row, T_dr1, T_s):
    signs = np.where(
        dR1_row < T_dr1, +1,
        np.where((sinner_row < T_s) & (mouter_row * (-dR1_row) > 0), +1, -1)
    )
    return signs

def balanced_accuracy(sign_pred, sign_true):
    tp = np.sum((sign_pred == +1) & (sign_true == +1))
    tn = np.sum((sign_pred == -1) & (sign_true == -1))
    fp = np.sum((sign_pred == +1) & (sign_true == -1))
    fn = np.sum((sign_pred == -1) & (sign_true == +1))
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.
    return 0.5 * (sens + spec)

def best_thresholds(dR1_train, sinner_train, mouter_train, Areq_train):
    """Find (T_dR1, T_s) maximising BA on training set."""
    sign_true = np.sign(Areq_train)
    best_ba, best_T = -1, (-0.131, 0.320)
    for t1 in T_DR1_SCAN:
        for t2 in T_S_SCAN:
            sp = apply_rule(dR1_train, sinner_train, mouter_train, t1, t2)
            ba = balanced_accuracy(sp, sign_true)
            if ba > best_ba:
                best_ba, best_T = ba, (t1, t2)
    return best_T

def loo_correlation(dR1_perm):
    """
    Run full LOO with threshold refitting on the given dR1 vector.
    Returns Pearson r between predicted and true Areq.
    """
    Areq_pred = np.zeros(n)
    sign_true = np.sign(Areq_true)

    for i in range(n):
        mask = np.arange(n) != i
        T_dr1, T_s = best_thresholds(
            dR1_perm[mask], sinner[mask], mouter[mask], Areq_true[mask]
        )
        sign_pred_i = apply_rule(
            dR1_perm[i:i+1], sinner[i:i+1], mouter[i:i+1], T_dr1, T_s
        )[0]
        # Amplitude proxy: sign × |Areq_true|  (sign is the bottleneck)
        Areq_pred[i] = sign_pred_i * abs(Areq_true[i])

    r, _ = pearsonr(Areq_true, Areq_pred)
    return r

# ---------------------------------------------------------------------------
# Observed LOO correlation (real, non-permuted data)
# ---------------------------------------------------------------------------
print("Computing observed LOO correlation (real data)...")
r_observed = loo_correlation(dR1_real)
print(f"  r_observed = {r_observed:.4f}")

# ---------------------------------------------------------------------------
# 500 permutations
# ---------------------------------------------------------------------------
rng = np.random.default_rng(args.seed)
null_r = np.zeros(args.n_perms)

print(f"\nRunning {args.n_perms} permutations (seed={args.seed})...")
for k in range(args.n_perms):
    dR1_perm = rng.permutation(dR1_real) # type: ignore
    null_r[k] = loo_correlation(dR1_perm)
    if (k + 1) % 50 == 0:
        print(f"  {k+1}/{args.n_perms}  null_r so far: "
              f"mean={null_r[:k+1].mean():.3f}  std={null_r[:k+1].std():.3f}")

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
Z     = (r_observed - null_r.mean()) / null_r.std()
p_val = (null_r >= r_observed).mean() # type: ignore
print(f"\nResults:")
print(f"  null mean = {null_r.mean():.4f}  std = {null_r.std():.4f}")
print(f"  r_observed = {r_observed:.4f}")
print(f"  Z = {Z:.2f}σ   p = {p_val:.4f}")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = pd.DataFrame({"r": null_r, "r_observed": r_observed})
outpath = os.path.join(DATA_DIR, "null_distribution.csv")
out.to_csv(outpath, index=False)
print(f"\nSaved: {outpath}  ({args.n_perms} rows)")
