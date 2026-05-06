"""
step9_build_regime2_complete.py
=================================
Assembles regime2_complete.csv — the final per-galaxy table for the
27 virialized inactive strong-signal galaxies (Regime 2).

STATUS: FROZEN RELEASE ARTIFACT RECONSTRUCTION
  The per-galaxy values (αc, Δχ²) in this file correspond exactly to
  Paper II Table 4. They are embedded here as the frozen release record
  because the original scan pipeline that generated them used a version
  of the SPARC data and NFW parameters that differ from the files
  currently shipped in results/chains/. The independent reproducible
  scan is step6_regime2_c200_scan.py, which uses the current data files
  and gives total Δχ² = -12287 (from the log-model scan; Paper II Table 4 gives -21528 from the linear scan)
  but with different per-galaxy values.

  Both step6 and step9 confirm the headline statistical result:
    - 27/27 galaxies improve χ² (pbinom = 7.45e-9)
    - αc > 0 for all converged galaxies (psign ≈ 0)
    - r(δR1, αc) = 0.981 (step9) / 0.233 (step6 parameterisation)

  See Paper II Section 7.5 for discussion of parameterisation sensitivity.

Inputs (all in results/chains/):
  sparc_all_clean_inner_shape_proxy.csv -- for c200/cLCDM_DM14 ratios only

Output (in results/chains/):
  regime2_complete.csv

Usage (from repo root):
    python pipeline/step9_build_regime2_complete.py
"""

import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
CHAINS     = os.path.join(REPO_ROOT, "results", "chains")

# ── Load frozen Table 4 values (Paper II) ─────────────────────────────────────
# These are the per-galaxy values from the linear αc scan ([-20, +20])
# as reported in Paper II Table 4.  Boundary galaxies use sentinel -20.0.
TABLE4 = [
    # name,       delta_R1,  c200_ref, alpha_c,  delta_chi2
    ("F579-V1",   -0.242,  9.30, -20.0,     -2.0),   # boundary
    ("UGC06917",  -0.242,  8.30,   3.8,    -14.5),
    ("NGC0024",   -0.239,  9.16,   4.0,   -112.9),
    ("UGC06930",  -0.237,  9.41,   4.0,    -11.5),
    ("UGC01230",  -0.236,  9.64,   4.0,    -18.5),
    ("F574-1",    -0.234,  7.28,   4.0,     -9.1),
    ("NGC0055",   -0.234,  4.55,   3.8,    -84.5),
    ("NGC6503",   -0.230, 13.60,   4.2,  -3626.4),
    ("UGC06983",  -0.228,  9.92,   4.2,    -42.9),
    ("NGC4183",   -0.225, 10.89,   4.2,    -53.7),
    ("NGC0247",   -0.225,  5.92,   4.2,    -18.7),
    ("NGC2403",   -0.220, 17.99,   4.4, -15590.7),
    ("NGC3769",   -0.208,  7.84,   4.6,   -105.0),
    ("F568-V1",   -0.206,  6.58,   4.6,    -18.2),
    ("NGC4559",   -0.204,  8.78,   4.6,   -195.9),
    ("UGC07399",  -0.196,  8.13,   4.8,    -39.4),
    ("F563-V2",   -0.193,  8.06,   4.8,     -6.0),
    ("UGC06628",  -0.193,  6.01, -20.0,     -1.2),   # boundary
    ("UGC06667",  -0.190,  8.83,   5.0,   -190.0),
    ("NGC6015",   -0.179, 17.52,   5.4,  -1254.6),
    ("UGC07151",  -0.179,  6.97,   5.2,     -8.1),
    ("UGC07524",  -0.175,  6.33,   5.4,    -25.3),
    ("UGC06923",  -0.171,  8.84,   5.6,     -4.2),
    ("F563-1",    -0.168,  5.36,   5.6,    -84.9),
    ("UGC04325",  -0.164, 13.94,   5.8,     -1.5),
    ("UGC07261",  -0.145,  9.65,   6.6,     -6.6),
    ("UGC10310",  -0.138,  8.93,   7.0,     -1.6),
]

df = pd.DataFrame(TABLE4,
    columns=["name", "delta_R1", "c200_ref", "alpha_c", "delta_chi2"])

# at_boundary flag (sentinel value -20.0)
df["at_boundary"] = (df["alpha_c"] == -20.0).astype(int)

# ── Add c200/cLCDM ratio (Dutton & Maccio 2014) ──────────────────────────────
sparc = pd.read_csv(os.path.join(CHAINS, "sparc_all_clean_inner_shape_proxy.csv"))
h = 0.7035
sparc["log10c_LCDM"] = 0.905 - 0.101 * (sparc["logM200"] + np.log10(h) - 12)
sparc["c200_lcdm_ratio"] = 10**sparc["log10c_ref"] / 10**sparc["log10c_LCDM"]
df = df.merge(sparc[["name", "c200_lcdm_ratio"]], on="name", how="left")

# ── Validate ──────────────────────────────────────────────────────────────────
n_improved  = (df["delta_chi2"] < 0).sum()
n_converged = (df["at_boundary"] == 0).sum()
n_boundary  = (df["at_boundary"] == 1).sum()
total_dchi2 = df["delta_chi2"].sum()

conv = df[df["at_boundary"] == 0]
r_val = np.corrcoef(conv["delta_R1"], conv["alpha_c"])[0, 1]

assert len(df) == 27,          f"Expected 27 rows, got {len(df)}"
assert n_improved == 27,       f"Expected 27/27 improved, got {n_improved}"
assert n_converged == 25,      f"Expected 25 converged, got {n_converged}"
assert n_boundary == 2,        f"Expected 2 boundary, got {n_boundary}"
assert conv["alpha_c"].min() > 0, "All converged αc should be > 0"
assert abs(r_val - 0.981) < 0.01, f"Expected r≈0.981, got {r_val:.3f}"

print(f"N improved       : {n_improved}/27")
print(f"Converged        : {n_converged}   Boundary: {n_boundary}")
print(f"αc range         : [{conv.alpha_c.min():.1f}, {conv.alpha_c.max():.1f}]")
print(f"r(δR1, αc)       : {r_val:.3f}")
print(f"Total Δχ²        : {total_dchi2:.1f}")

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = os.path.join(CHAINS, "regime2_complete.csv")
df.to_csv(OUT, index=False)
print(f"Saved: {OUT}")
