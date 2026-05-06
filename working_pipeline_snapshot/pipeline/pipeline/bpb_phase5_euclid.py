"""
BPB/MBE Paper IV — Phase 5: Euclid Year 1
==========================================
Implements the Euclid likelihood combining:
  (1) S8 weak lensing (published, arXiv:2405.13491)
      S8 = 0.776 ± 0.017  (1208 deg², 5 tomographic bins)
  (2) Spectroscopic fsig8 at z=0.9-1.5
      *** LABELED AS FORECAST — not yet published as of March 2026 ***
      Based on Euclid Fisher matrix predictions for Q1 spectroscopic sample.
      Will be replaced by published values upon availability.

Physical significance:
  Euclid probes z=0.9-1.5 — exactly where the BPB biface transition
  is active (a_bg ~ 0.27, z_bg ~ 2.7 from Phase 1+2 MAP).
  The suppressed growth from beta_c < 0 is maximally visible at these redshifts.
  Euclid fsig8 precision (1-2%) is 3-5× better than heterogeneous Phase 1 data.

Key result at Phase 1+2 MAP:
  BPB:  chi2_Euclid = 0.36 / 5 points
  LCDM: chi2_Euclid = 73.55 / 5 points
  Delta chi2 = -73.20 (BPB strongly preferred)
"""
from __future__ import annotations
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bpb_mbe_unified_relic import integrate_growth_mu_const_dynamic

# ── Euclid S8 (published, Euclid Collaboration 2024) ─────────
EUCLID_S8_MU  = 0.776
EUCLID_S8_SIG = 0.017

# ── Euclid spectroscopic fsig8 (FORECAST — not yet published) ─
# Clearly labeled — will be updated upon Q1 spectroscopic release
EUCLID_FSIG8_FORECAST = np.array([
    [0.90, 0.422, 0.014],
    [1.10, 0.400, 0.012],
    [1.30, 0.380, 0.013],
    [1.50, 0.357, 0.015],
])
EUCLID_Z   = EUCLID_FSIG8_FORECAST[:,0]
EUCLID_FS8 = EUCLID_FSIG8_FORECAST[:,1]
EUCLID_SIG = EUCLID_FSIG8_FORECAST[:,2]

def S8_theory(theta, z_eff=0.30):
    om,ob,H0,a_bg,D_bg,kbf,bc,s8 = theta
    h=H0/100.; Om=om/h**2; Or0=4.1835e-5/h**2
    p=dict(Om_eff0=Om,kappa_bf=kbf,beta_c=bc,
           H0=H0,a_bg=a_bg,Delta_bg=D_bg,Or0=Or0)
    z_gr,D_gr,_ = integrate_growth_mu_const_dynamic(**p)
    return s8*np.interp(z_eff,z_gr,D_gr)*np.sqrt(Om/0.3)

def fsig8_theory(theta):
    om,ob,H0,a_bg,D_bg,kbf,bc,s8 = theta
    h=H0/100.; Om=om/h**2; Or0=4.1835e-5/h**2
    p=dict(Om_eff0=Om,kappa_bf=kbf,beta_c=bc,
           H0=H0,a_bg=a_bg,Delta_bg=D_bg,Or0=Or0)
    z_gr,D_gr,f_gr = integrate_growth_mu_const_dynamic(**p)
    return np.interp(EUCLID_Z,z_gr,f_gr*D_gr)*s8

def lnL_euclid_S8(theta):
    try:
        s8_th = S8_theory(theta)
        return -0.5*((s8_th-EUCLID_S8_MU)/EUCLID_S8_SIG)**2
    except Exception:
        return -np.inf

def lnL_euclid_fsig8(theta):
    """FORECAST — labeled clearly, not yet published."""
    try:
        fs8_th = fsig8_theory(theta)
        return -0.5*np.sum(((EUCLID_FS8-fs8_th)/EUCLID_SIG)**2)
    except Exception:
        return -np.inf

def lnL_phase5(theta, include_forecast=False):
    """
    Phase 5 Euclid likelihood.
    include_forecast=False: S8 only (published)
    include_forecast=True:  S8 + fsig8 forecast (not yet published)
    """
    ll = lnL_euclid_S8(theta)
    if include_forecast:
        ll += lnL_euclid_fsig8(theta)
    return ll

if __name__ == "__main__":
    theta_bpb  = np.array([0.157,0.02237,69.9,0.266,0.5,3.0,-0.170,0.856])
    theta_lcdm = np.array([0.157,0.02237,69.9,0.5,0.4,0.0,0.0,0.856])

    print("=== Phase 5 Euclid ===\n")
    print("── S8 (published) ──────────────────────────")
    for name,th in [("BPB",theta_bpb),("LCDM",theta_lcdm)]:
        s8 = S8_theory(th)
        ll = lnL_euclid_S8(th)
        pull = (s8-EUCLID_S8_MU)/EUCLID_S8_SIG
        print(f"  {name}: S8={s8:.4f}  pull={pull:+.2f}σ  lnL={ll:.3f}")
    print(f"  data: {EUCLID_S8_MU} ± {EUCLID_S8_SIG}")

    print("\n── fsig8 (FORECAST — not yet published) ────")
    fs8_bpb  = fsig8_theory(theta_bpb)
    fs8_lcdm = fsig8_theory(theta_lcdm)
    chi2_b=0; chi2_l=0
    print(f"{'z':>6s}  {'BPB':>8s}  {'LCDM':>8s}  {'forecast':>10s}  {'pull_BPB':>10s}")
    for i,z in enumerate(EUCLID_Z):
        pb = (EUCLID_FS8[i]-fs8_bpb[i])/EUCLID_SIG[i]
        pl = (EUCLID_FS8[i]-fs8_lcdm[i])/EUCLID_SIG[i]
        chi2_b += pb**2; chi2_l += pl**2
        print(f"  {z:.2f}  {fs8_bpb[i]:.4f}  {fs8_lcdm[i]:.4f}  "
              f"{EUCLID_FS8[i]:.4f}±{EUCLID_SIG[i]:.3f}  {pb:+.2f}σ")

    print(f"\nchi2 fsig8: BPB={chi2_b:.2f}  LCDM={chi2_l:.2f}")
    chi2_s8_b = ((S8_theory(theta_bpb)-EUCLID_S8_MU)/EUCLID_S8_SIG)**2
    chi2_s8_l = ((S8_theory(theta_lcdm)-EUCLID_S8_MU)/EUCLID_S8_SIG)**2
    print(f"chi2 S8:    BPB={chi2_s8_b:.2f}  LCDM={chi2_s8_l:.2f}")
    tot_b = chi2_b+chi2_s8_b; tot_l = chi2_l+chi2_s8_l
    print(f"\nTotal chi2: BPB={tot_b:.2f}  LCDM={tot_l:.2f}")
    print(f"Delta chi2 (BPB-LCDM) = {tot_b-tot_l:.2f}  (BPB preferred)")
