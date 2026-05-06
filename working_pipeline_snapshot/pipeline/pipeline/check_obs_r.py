#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.cwd()
f = ROOT / "results" / "chains" / "hybrid_loo_per_galaxy.csv"
df = pd.read_csv(f)

x = df["A_true"].to_numpy(dtype=float)
y = df["A_pred_hybrid_loo"].to_numpy(dtype=float)

m = np.isfinite(x) & np.isfinite(y)
x = x[m]
y = y[m]

def pearson_r(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float((a @ b) / np.sqrt((a @ a) * (b @ b)))

r_obs = pearson_r(x, y)
print("observed_r =", r_obs)