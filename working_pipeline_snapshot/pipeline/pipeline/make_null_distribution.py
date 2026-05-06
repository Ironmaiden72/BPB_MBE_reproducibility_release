#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

try:
    from scipy.stats import norm
except ImportError:
    raise SystemExit("Installe scipy: pip install scipy")

ROOT = Path.cwd()
CHAINS = ROOT / "results" / "chains"
INFILE = CHAINS / "hybrid_loo_per_galaxy.csv"
OUTFILE = CHAINS / "null_distribution.csv"

if not INFILE.exists():
    raise FileNotFoundError(f"Fichier introuvable: {INFILE}")

df = pd.read_csv(INFILE)

# Adapte ici si les noms exacts diffèrent
x_col = "dR1"
y_col = "Areq_loo"

missing = [c for c in [x_col, y_col] if c not in df.columns]
if missing:
    raise ValueError(
        f"Colonnes manquantes {missing}. Colonnes dispo: {list(df.columns)}"
    )

x = df[x_col].to_numpy(dtype=float)
y = df[y_col].to_numpy(dtype=float)

mask = np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]

if len(x) < 3:
    raise ValueError("Pas assez de points valides pour calculer une corrélation.")

def pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    da = np.sqrt(np.sum(a * a))
    db = np.sqrt(np.sum(b * b))
    if da == 0 or db == 0:
        return np.nan
    return float(np.sum(a * b) / (da * db))

r_obs = pearson_r(x, y)

rng = np.random.default_rng(0)
nperm = 500
r_vals = np.empty(nperm, dtype=float)

for i in range(nperm):
    x_perm = rng.permutation(x)
    r_vals[i] = pearson_r(x_perm, y)

# test unilatéral: proportion des permutations >= r observé
p_val = (1 + np.sum(r_vals >= r_obs)) / (1 + nperm)
z_val = float(norm.isf(p_val))

out = pd.DataFrame({
    "r_value": r_vals,
    "observed_r": np.full(nperm, r_obs),
    "p_value": np.full(nperm, p_val),
    "z_value": np.full(nperm, z_val),
})

OUTFILE.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTFILE, index=False)

print(f"OK -> {OUTFILE}")
print(f"observed_r = {r_obs:.6f}")
print(f"p_value    = {p_val:.6f}")
print(f"z_value    = {z_val:.6f}")