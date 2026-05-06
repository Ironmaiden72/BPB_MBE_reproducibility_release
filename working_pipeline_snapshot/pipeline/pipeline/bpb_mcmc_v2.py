"""
BPB/MBE Paper IV — Phase 1 MCMC  (uses bpb_mbe_unified_relic.py)
=================================================================
Dependencies:  numpy scipy emcee corner matplotlib
               bpb_mbe_unified_relic.py  (must be in same directory)

Install:   pip install emcee corner matplotlib scipy numpy h5py

Usage:
    python3 bpb_mcmc_v2.py --test       # sanity check, no MCMC
    python3 bpb_mcmc_v2.py --bestfit    # MAP + diagnostics + figures
    python3 bpb_mcmc_v2.py              # full run (32 walkers × 5000 steps)
    python3 bpb_mcmc_v2.py --plot bpb_phase1_chain.h5

Parameters:  theta = [omega_m, omega_b, H0, a_bg, Delta_bg, kappa_bf, beta_c, sigma8_0]

  omega_m   : Omega_m * h^2
  omega_b   : Omega_b * h^2       (BBN Gaussian prior)
  H0        : Hubble constant [km/s/Mpc]
  a_bg      : transition scale factor  (0 < a_bg < 1)
  Delta_bg  : transition width in ln(a)
  kappa_bf  : BPB interaction coupling  (>= 0)
  beta_c    : growth enhancement  mu = 1 + beta_c  (> -1)
  sigma8_0  : sigma_8 amplitude today

Physical model (bpb_mbe_unified_relic):
  E²(z) = Or0/a⁴ + Om_surv0/a³ + F(a)*(Om_induced0/a³ + Ode_eff0)
  F(a)  = S(a)/S(1),   S(a) = 0.5*(1+tanh((ln a - ln a_bg)/Delta_bg))
  Growth: dD²/d(lna)² + (2+dlnH/dlna)*dD/d(lna) = 1.5*Om_src(a)*mu*D
  with mu = 1 + beta_c  (scale-independent)
  fsig8(z) = f(z) * sigma8_0 * D(z)/D(0)
"""

import sys
import time
import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

# ── Phase 2-5 optional imports ──────────────────────────────────────────────
try:
    from bpb_phase2_kids import lnL_phase2 as _lnL_p2
    HAS_PHASE2 = True
except ImportError:
    HAS_PHASE2 = False

try:
    from bpb_phase3_act import lnL_phase3 as _lnL_p3
    HAS_PHASE3 = True
except ImportError:
    HAS_PHASE3 = False

try:
    from bpb_phase4_desi import lnL_phase4 as _lnL_p4
    HAS_PHASE4 = True
except ImportError:
    HAS_PHASE4 = False

try:
    from bpb_phase5_euclid import lnL_phase5 as _lnL_p5
    HAS_PHASE5 = True
except ImportError:
    HAS_PHASE5 = False

_phases = {k:v for k,v in {
    'Phase2 (KiDS S8)':    HAS_PHASE2 if 'HAS_PHASE2' in dir() else False,
    'Phase3 (ACT)':        HAS_PHASE3 if 'HAS_PHASE3' in dir() else False,
    'Phase4 (DESI)':       HAS_PHASE4 if 'HAS_PHASE4' in dir() else False,
    'Phase5 (Euclid S8)':  HAS_PHASE5 if 'HAS_PHASE5' in dir() else False,
}.items()}

# ── the actual BPB/MBE module ────────────────────────────────────────────────
try:
    from bpb_mbe_unified_relic import (
        H_unified_relic_H0,
        DM_mpc_unified_relic,
        integrate_growth_mu_const_dynamic,
        E2_unified_relic,
    )
    HAS_RELIC = True
except ImportError:
    HAS_RELIC = False
    print("ERROR: bpb_mbe_unified_relic.py not found in current directory.")
    sys.exit(1)

try:
    import emcee
    HAS_EMCEE = True
except ImportError:
    HAS_EMCEE = False
    print("WARNING: emcee not found.  pip install emcee h5py")

try:
    import corner
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

C_LIGHT = 299792.458   # km/s

# ============================================================
# Data
# ============================================================

# ── CC  (z, H_obs, sigma_H)  [km/s/Mpc] ─────────────────────
CC_DATA = np.array([
    [0.070,  69.0, 19.6], [0.090,  69.0, 12.0],
    [0.120,  68.6, 26.2], [0.170,  83.0,  8.0],
    [0.179,  75.0,  4.0], [0.199,  75.0,  5.0],
    [0.200,  72.9, 29.6], [0.270,  77.0, 14.0],
    [0.280,  88.8, 36.6], [0.352,  83.0, 14.0],
    [0.380,  83.0, 13.5], [0.400,  95.0, 17.0],
    [0.4004, 77.0, 10.2], [0.4247, 87.1, 11.2],
    [0.4497, 92.8, 12.9], [0.4783, 80.9,  9.0],
    [0.480,  97.0, 62.0], [0.593, 104.0, 13.0],
    [0.680,  92.0,  8.0], [0.781, 105.0, 12.0],
    [0.875, 125.0, 17.0], [0.880,  90.0, 40.0],
    [0.900, 117.0, 23.0], [1.037, 154.0, 20.0],
    [1.300, 168.0, 17.0], [1.363, 160.0, 33.6],
    [1.430, 177.0, 18.0], [1.530, 140.0, 14.0],
    [1.750, 202.0, 40.0], [1.965, 186.5, 50.4],
    [0.3802, 83.0, 13.5],
])

# ── BAO simple  (z, type, value, sigma) ──────────────────────
# type: 0=DV/rs  1=DM/rs  2=DH/rs
BAO_SIMPLE = np.array([
    [0.106, 0,  3.047, 0.137],
    [0.150, 0,  4.480, 0.168],
    [1.480, 1, 30.21,  0.79 ],
    [1.480, 2, 13.23,  0.47 ],
])

# ── BOSS DR12 (Alam+2017) 6-vector ───────────────────────────
BOSS_Z    = np.array([0.38, 0.51, 0.61])
BOSS_OBS  = np.array([10.23, 24.98, 13.36, 22.33, 15.64, 20.93])
BOSS_SIGMA = np.array([0.22, 0.60, 0.25, 0.55, 0.28, 0.52])
BOSS_CORR  = np.array([
    [1.000, 0.000, 0.387, 0.000, 0.067, 0.000],
    [0.000, 1.000, 0.000, 0.415, 0.000, 0.036],
    [0.387, 0.000, 1.000, 0.000, 0.438, 0.000],
    [0.000, 0.415, 0.000, 1.000, 0.000, 0.466],
    [0.067, 0.000, 0.438, 0.000, 1.000, 0.000],
    [0.000, 0.036, 0.000, 0.466, 0.000, 1.000],
])
BOSS_COV     = BOSS_CORR * np.outer(BOSS_SIGMA, BOSS_SIGMA)
BOSS_COV_CHO = cho_factor(BOSS_COV)

# ── SNe Ia (Pantheon+ compressed) ────────────────────────────
SNE_Z   = np.array([0.10, 0.30, 0.60, 1.00])
SNE_MU  = np.array([38.32, 40.65, 42.48, 44.08])
SNE_COV = np.array([
    [0.0036, 0.0009, 0.0004, 0.0002],
    [0.0009, 0.0025, 0.0006, 0.0003],
    [0.0004, 0.0006, 0.0020, 0.0004],
    [0.0002, 0.0003, 0.0004, 0.0016],
])
SNE_COV_CHO = cho_factor(SNE_COV)

# ── f*sigma8  (z_eff, fsig8, sigma) ──────────────────────────
FSIG8_DATA = np.array([
    [0.067, 0.423, 0.055], [0.170, 0.510, 0.060],
    [0.220, 0.420, 0.070], [0.250, 0.351, 0.058],
    [0.320, 0.384, 0.095], [0.350, 0.440, 0.050],
    [0.370, 0.460, 0.038], [0.380, 0.497, 0.045],
    [0.410, 0.450, 0.040], [0.510, 0.458, 0.038],
    [0.570, 0.441, 0.043], [0.600, 0.430, 0.040],
    [0.600, 0.550, 0.120], [0.610, 0.436, 0.034],
    [0.700, 0.473, 0.041], [0.770, 0.490, 0.180],
    [0.780, 0.380, 0.040], [0.850, 0.315, 0.095],
    [0.860, 0.400, 0.110], [1.400, 0.482, 0.116],
    [1.480, 0.462, 0.045], [1.520, 0.420, 0.076],
])

# ============================================================
# Parameter layout
# ============================================================

PARAM_NAMES = ['omega_m','omega_b','H0','a_bg','Delta_bg','kappa_bf','beta_c','sigma8_0']
NDIM        = len(PARAM_NAMES)

LABELS = [r'$\omega_m$', r'$\omega_b$', r'$H_0$',
          r'$a_{\rm bg}$', r'$\Delta_{\rm bg}$',
          r'$\kappa_{\rm bf}$', r'$\beta_c$', r'$\sigma_8^{(0)}$']

# Flat bounds
PRIOR_BOUNDS = np.array([
    [0.10,  0.20 ],   # omega_m
    [0.019, 0.025],   # omega_b
    [60.0,  80.0 ],   # H0
    [0.10,  0.90 ],   # a_bg  (transition epoch; a=1 today)
    [0.05,  2.0  ],   # Delta_bg (ln-a width)
    [0.0,   5.0  ],   # kappa_bf (>= 0)
    [-0.5,  0.5  ],   # beta_c  (mu = 1+beta_c > 0)
    [0.65,  1.05 ],   # sigma8_0
])

OB_MU, OB_SIG = 0.02237, 0.00015   # BBN Gaussian on omega_b

# Starting point
THETA0 = np.array([0.143, 0.02237, 67.4, 0.5, 0.4, 1.0, 0.0, 0.81])

# ============================================================
# Helper: Or0 and sound horizon
# ============================================================

def Or0_from_omega(omega_b, H0):
    """Radiation density (photons + massless nu approx)."""
    h = H0 / 100.0
    return 4.1835e-5 / h**2

def rs_EH(omega_m, omega_b):
    """Sound horizon at drag epoch, Aubourg+2015 fit [Mpc]."""
    nu = 0.00064
    return 55.154 * np.exp(-72.3*(nu+0.0006)**2) / (
        omega_m**0.25351 * omega_b**0.12807)

# ============================================================
# Model wrapper
# ============================================================

def unpack(theta):
    om, ob, H0, a_bg, D_bg, kbf, bc, s8 = theta
    h    = H0 / 100.0
    Om   = om / h**2
    Or0  = Or0_from_omega(ob, H0)
    return dict(Om_eff0=Om, kappa_bf=kbf, beta_c=bc,
                H0=H0, a_bg=a_bg, Delta_bg=D_bg, Or0=Or0), s8

def H_z(theta, z):
    p, _ = unpack(theta)
    return H_unified_relic_H0(np.atleast_1d(z), **p)

def DM_z(theta, z):
    p, _ = unpack(theta)
    return DM_mpc_unified_relic(np.atleast_1d(z), **p)

def DH_z(theta, z):
    p, _ = unpack(theta)
    return C_LIGHT / H_unified_relic_H0(np.atleast_1d(z), **p)

def DV_z(theta, z):
    z   = float(z)
    dm  = float(DM_z(theta, np.array([z]))[0])
    dh  = float(DH_z(theta, np.array([z]))[0])
    return (z * dm**2 * dh)**(1.0/3.0)

def mu_z(theta, z):
    dm = DM_z(theta, np.array([float(z)]))[0]
    return 5.0 * np.log10((1+z) * dm / 1e-5)

def growth_arrays(theta):
    """Returns (z_arr, D_arr, f_arr) from integrate_growth_mu_const_dynamic."""
    p, _ = unpack(theta)
    # Remove H0 from p since integrate_growth takes it as keyword
    return integrate_growth_mu_const_dynamic(**p)

def fsig8_at(theta, z_eff):
    """fsig8(z) = f(z) * sigma8_0 * D(z)/D(0)."""
    _, s8   = unpack(theta)
    z_arr, D_arr, f_arr = growth_arrays(theta)
    D  = np.interp(z_eff, z_arr, D_arr)
    f  = np.interp(z_eff, z_arr, f_arr)
    return f * D * s8

# ============================================================
# Likelihoods
# ============================================================

def lnL_CC(theta):
    z, H_obs, sH = CC_DATA[:,0], CC_DATA[:,1], CC_DATA[:,2]
    H_th = H_z(theta, z)
    return -0.5 * np.sum(((H_obs - H_th) / sH)**2)

def lnL_BAO(theta):
    om, ob = theta[0], theta[1]
    rs  = rs_EH(om, ob)
    lnL = 0.0
    for z, t, val, sig in BAO_SIMPLE:
        if   t == 0: th = DV_z(theta, z) / rs
        elif t == 1: th = float(DM_z(theta, np.array([z]))[0]) / rs
        else:        th = float(DH_z(theta, np.array([z]))[0]) / rs
        lnL -= 0.5 * ((val - th) / sig)**2
    tv = []
    for zi in BOSS_Z:
        tv.append(float(DM_z(theta, np.array([zi]))[0]) / rs)
        tv.append(float(DH_z(theta, np.array([zi]))[0]) / rs)
    d   = BOSS_OBS - np.array(tv)
    lnL -= 0.5 * d @ cho_solve(BOSS_COV_CHO, d)
    return lnL

def lnL_SN(theta):
    mu_th  = np.array([mu_z(theta, z) for z in SNE_Z])
    d      = SNE_MU - mu_th
    Cinv_d = cho_solve(SNE_COV_CHO, d)
    A = np.sum(Cinv_d)
    B = np.sum(cho_solve(SNE_COV_CHO, np.ones(len(d))))
    return -0.5 * (d @ Cinv_d - A**2 / B)

def lnL_fsig8(theta):
    _, s8   = unpack(theta)
    z_arr, D_arr, f_arr = growth_arrays(theta)
    lnL = 0.0
    for z_eff, fs8_obs, sig in FSIG8_DATA:
        D     = np.interp(z_eff, z_arr, D_arr)
        f     = np.interp(z_eff, z_arr, f_arr)
        fs8_th = f * D * s8
        lnL   -= 0.5 * ((fs8_obs - fs8_th) / sig)**2
    return lnL

def lnL_total(theta):
    ll = lnL_CC(theta) + lnL_BAO(theta) + lnL_SN(theta) + lnL_fsig8(theta)
    if HAS_PHASE2: ll += _lnL_p2(theta)   # KiDS S8
    if HAS_PHASE3: ll += _lnL_p3(theta)   # ACT + Planck lensing
    if HAS_PHASE4: ll += _lnL_p4(theta)   # DESI DR1 fsig8
    if HAS_PHASE5: ll += _lnL_p5(theta)   # Euclid S8 (published)
    return ll

# ============================================================
# LCDM baseline  (beta_c=0, kappa_bf=0, a_bg irrelevant)
# ============================================================

def lnL_lcdm(theta):
    """
    LCDM likelihood — optimises over (H0, sigma8_0) at fixed omega_m, omega_b
    to get a fair comparison. Uses BPB omega_m, omega_b as starting point.
    """
    om, ob, H0, _, _, _, _, s8 = theta
    # LCDM has kappa_bf=0, beta_c=0; optimise H0 and sigma8 for fair comparison
    def neg_lnL(x):
        H0_l, s8_l = x
        if not (60 < H0_l < 80) or not (0.6 < s8_l < 1.1):
            return np.inf
        t = np.array([om, ob, H0_l, 0.5, 0.3, 0.0, 0.0, s8_l])
        try:
            return -lnL_total(t)
        except Exception:
            return np.inf
    res = minimize(neg_lnL, [H0, s8], method='Nelder-Mead',
                   options={'maxiter': 500, 'xatol': 1e-3, 'fatol': 1e-3})
    return -res.fun

def lnL_lcdm_fixed(theta):
    """LCDM at same parameters — fast version for diagnostics."""
    om, ob, H0, _, _, _, _, s8 = theta
    theta_lcdm = np.array([om, ob, H0, 0.5, 0.3, 0.0, 0.0, s8])
    return lnL_total(theta_lcdm)

# ============================================================
# Prior + posterior
# ============================================================

def ln_prior(theta):
    for i, (lo, hi) in enumerate(PRIOR_BOUNDS):
        if not (lo <= theta[i] <= hi):
            return -np.inf
    # beta_c: mu = 1 + beta_c must be > 0
    if 1.0 + theta[6] <= 0:
        return -np.inf
    # kappa_bf >= 0 (already enforced by bounds)
    ob = theta[1]
    return -0.5 * ((ob - OB_MU) / OB_SIG)**2

def ln_posterior(theta):
    lp = ln_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    try:
        ll = lnL_total(theta)
        if not np.isfinite(ll):
            return -np.inf
        return lp + ll
    except Exception:
        return -np.inf

# ============================================================
# MAP
# ============================================================

def find_map(theta0=THETA0, verbose=True):
    if verbose:
        print("Finding MAP (Nelder-Mead)...")
    t0  = time.time()
    res = minimize(lambda t: -ln_posterior(t), theta0,
                   method='Nelder-Mead',
                   options={'maxiter':8000,'xatol':1e-4,'fatol':1e-4,
                            'adaptive':True})
    if verbose:
        print(f"  MAP lnP = {-res.fun:.3f}  ({time.time()-t0:.1f}s)")
        for n, v in zip(PARAM_NAMES, res.x):
            print(f"  {n:12s} = {v:.5f}")
    return res.x, -res.fun

# ============================================================
# Diagnostics
# ============================================================

def print_diagnostics(theta_map):
    ll_bpb  = lnL_total(theta_map)
    ll_lcdm = lnL_lcdm_fixed(theta_map)   # evaluated at BPB params (indicative)
    dchi2   = -2.0 * (ll_bpb - ll_lcdm)

    om, ob, H0, a_bg, D_bg, kbf, bc, s8 = theta_map
    z_t     = 1.0/a_bg - 1.0
    mu_grav = 1.0 + bc

    print("\n" + "="*58)
    print("DIAGNOSTICS AT MAP")
    print("="*58)
    print(f"lnL BPB   = {ll_bpb:.3f}")
    print(f"lnL LCDM  = {ll_lcdm:.3f}")
    print(f"Delta chi^2 (BPB - LCDM) = {dchi2:+.3f}  (neg = BPB preferred)")
    print(f"\nTransition: a_bg={a_bg:.3f}  -> z_t = {z_t:.3f}")
    print(f"Width:      Delta_bg={D_bg:.3f} (in ln a)")
    print(f"Coupling:   kappa_bf={kbf:.3f},  beta_c={bc:.4f}")
    print(f"Growth mu = 1+beta_c = {mu_grav:.4f}")
    print(f"r_s = {rs_EH(om,ob):.3f} Mpc  (expect ~147)")

    # fsig8 residuals
    z_arr, D_arr, f_arr = growth_arrays(theta_map)
    chi2_bpb, chi2_lcdm = 0.0, 0.0
    print(f"\nfsig8 residuals (obs - th) / sigma:")
    for z_eff, fs8_obs, sig in FSIG8_DATA:
        D      = np.interp(z_eff, z_arr, D_arr)
        f      = np.interp(z_eff, z_arr, f_arr)
        fs8_th = f * D * s8
        pull   = (fs8_obs - fs8_th) / sig
        chi2_bpb += pull**2
        flag   = "  <--" if abs(pull) > 2.0 else ""
        print(f"  z={z_eff:.3f}  obs={fs8_obs:.3f}  th={fs8_th:.3f}"
              f"  pull={pull:+.2f}{flag}")

    # LCDM baseline: same omega_m, H0, sigma8 but beta_c=kappa_bf=0
    theta_lcdm = np.array([om, ob, H0, 0.5, 0.3, 0.0, 0.0, s8])
    z_arr_l, D_arr_l, f_arr_l = growth_arrays(theta_lcdm)
    for z_eff, fs8_obs, sig in FSIG8_DATA:
        D      = np.interp(z_eff, z_arr_l, D_arr_l)
        f      = np.interp(z_eff, z_arr_l, f_arr_l)
        chi2_lcdm += ((fs8_obs - f*D*s8) / sig)**2

    print(f"\n  BPB  fsig8 chi^2 = {chi2_bpb:.2f}")
    print(f"  LCDM fsig8 chi^2 = {chi2_lcdm:.2f}")
    print(f"  Delta chi^2 fsig8 alone = {chi2_bpb - chi2_lcdm:+.2f}")

# ============================================================
# emcee
# ============================================================

def run_emcee(theta_map, n_walkers=32, n_steps=5000, n_burn=500,
              backend_file="bpb_allphases_chain.h5", verbose=True):
    if not HAS_EMCEE:
        raise ImportError("pip install emcee h5py")

    ndim   = len(theta_map)
    scales = 0.005 * (PRIOR_BOUNDS[:,1] - PRIOR_BOUNDS[:,0])
    p0     = theta_map + scales * np.random.randn(n_walkers, ndim)
    for iw in range(n_walkers):
        for ip in range(ndim):
            lo, hi = PRIOR_BOUNDS[ip]
            p0[iw, ip] = np.clip(p0[iw, ip], lo+1e-6, hi-1e-6)

    backend = emcee.backends.HDFBackend(backend_file)
    backend.reset(n_walkers, ndim)
    sampler = emcee.EnsembleSampler(n_walkers, ndim, ln_posterior,
                                    backend=backend)
    if verbose:
        print(f"\nRunning emcee: {n_walkers} walkers × {n_steps} steps")
        t0 = time.time()

    state = sampler.run_mcmc(p0, n_burn, progress=verbose, store=False)
    sampler.reset()
    sampler.run_mcmc(state, n_steps, progress=verbose)

    if verbose:
        print(f"\nDone in {(time.time()-t0)/60:.1f} min")
        tau = sampler.get_autocorr_time(quiet=True)
        print(f"Autocorr time: {np.round(tau, 1)}")
        print(f"Acceptance:    {np.mean(sampler.acceptance_fraction):.3f}")
    return sampler

# ============================================================
# Summary + plots
# ============================================================

def print_summary(chain):
    print("\n" + "="*62)
    print("MCMC SUMMARY")
    print("="*62)
    print(f"{'Parameter':12s}  {'Mean':>10s}  {'Std':>9s}  "
          f"{'16%':>10s}  {'84%':>10s}")
    print("-"*62)
    for i, name in enumerate(PARAM_NAMES):
        v = chain[:, i]
        print(f"{name:12s}  {np.mean(v):10.4f}  {np.std(v):9.4f}"
              f"  {np.percentile(v,16):10.4f}  {np.percentile(v,84):10.4f}")

def make_corner(chain, out="bpb_phase1_corner.pdf"):
    if not HAS_PLOT:
        print("SKIP corner plot: install corner+matplotlib  (pip install corner matplotlib)")
        return
    fig = corner.corner(chain, labels=LABELS,
                        quantiles=[0.16,0.5,0.84], show_titles=True)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Corner: {out}")
    plt.close()

def make_diagnostics_plot(theta_map, out="bpb_phase1_diagnostics.pdf"):
    if not HAS_PLOT:
        print("SKIP diagnostics plot: install matplotlib  (pip install matplotlib)")
        return
    z_bg = np.linspace(0.001, 2.5, 200)
    p, s8  = unpack(theta_map)
    p_lcdm = dict(Om_eff0=p['Om_eff0'], kappa_bf=0.0, beta_c=0.0,
                  H0=p['H0'], a_bg=0.5, Delta_bg=0.3, Or0=p['Or0'])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # H(z)
    ax = axes[0]
    ax.plot(z_bg, H_unified_relic_H0(z_bg, **p),   'b-',  lw=1.8, label='BPB/MBE')
    ax.plot(z_bg, H_unified_relic_H0(z_bg, **p_lcdm), 'k--', lw=1.4, label=r'$\Lambda$CDM')
    ax.errorbar(CC_DATA[:,0], CC_DATA[:,1], CC_DATA[:,2],
                fmt='o', color='gray', ms=3, lw=0.8, label='CC')
    ax.set_xlabel('z'); ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.legend(fontsize=8); ax.set_title('Background')

    # fsig8
    ax = axes[1]
    z_arr,  D_arr,  f_arr  = integrate_growth_mu_const_dynamic(**p)
    z_arrl, D_arrl, f_arrl = integrate_growth_mu_const_dynamic(**p_lcdm)
    fs8_bpb  = np.interp(z_bg, z_arr,  f_arr  * D_arr)  * s8
    fs8_lcdm = np.interp(z_bg, z_arrl, f_arrl * D_arrl) * s8
    ax.plot(z_bg, fs8_bpb,  'b-',  lw=1.8, label='BPB/MBE')
    ax.plot(z_bg, fs8_lcdm, 'k--', lw=1.4, label=r'$\Lambda$CDM')
    ax.errorbar(FSIG8_DATA[:,0], FSIG8_DATA[:,1], FSIG8_DATA[:,2],
                fmt='s', color='C1', ms=3, lw=0.8, label=r'$f\sigma_8$')
    ax.set_xlabel('z'); ax.set_ylabel(r'$f\sigma_8(z)$')
    ax.legend(fontsize=8); ax.set_title('Growth')

    # Transition function F(a)
    from bpb_mbe_unified_relic import F_bg
    ax = axes[2]
    a_arr = 1.0/(1.0+z_bg)
    F_arr = F_bg(a_arr, a_bg=p['a_bg'], Delta_bg=p['Delta_bg'])
    ax.plot(z_bg, F_arr, 'b-', lw=1.8)
    ax.axvline(1.0/p['a_bg']-1.0, color='gray', ls='--', lw=1, label=r'$z_{\rm bg}$')
    ax.set_xlabel('z'); ax.set_ylabel(r'$F(a)=S(a)/S(1)$')
    ax.legend(fontsize=8); ax.set_title('Activation function')

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Diagnostics: {out}")
    plt.close()

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "--full"
    np.random.seed(42)

    # ── test ─────────────────────────────────────────────────
    if mode == "--test":
        print("=== TEST MODE ===")
        p, s8 = unpack(THETA0)
        H0_val = float(H_unified_relic_H0(np.array([0.0]), **p)[0])
        print(f"H(0)   = {H0_val:.2f}  (expect ~67.4)")
        print(f"r_s    = {rs_EH(THETA0[0], THETA0[1]):.2f}  (expect ~147)")
        dm_106 = float(DM_mpc_unified_relic(np.array([0.106]), **p)[0])
        dh_106 = C_LIGHT / float(H_unified_relic_H0(np.array([0.106]), **p)[0])
        dv_106 = (0.106 * dm_106**2 * dh_106)**(1/3)
        print(f"DV(0.106) = {dv_106:.2f} Mpc")

        lp = ln_posterior(THETA0)
        print(f"ln_posterior(theta0) = {lp:.3f}")
        print("\nlnL components:")
        for fn, nm in [(lnL_CC,'CC'),(lnL_BAO,'BAO'),(lnL_SN,'SNe'),(lnL_fsig8,'fsig8')]:
            print(f"  {nm:6s}: {fn(THETA0):.3f}")

        # fsig8 check
        z_arr, D_arr, f_arr = growth_arrays(THETA0)
        print("\nfsig8 sample (BPB, beta_c=0 at theta0):")
        for z in [0.1, 0.3, 0.5, 0.8, 1.0]:
            D  = np.interp(z, z_arr, D_arr)
            f  = np.interp(z, z_arr, f_arr)
            print(f"  z={z:.1f}  fsig8={f*D*s8:.4f}")

        if not HAS_EMCEE:
            print("\nInstall emcee:  pip install emcee h5py")
            sys.exit(0)

        theta_map, _ = find_map(verbose=True)
        print_diagnostics(theta_map)
        # Short test chain (16 walkers min for ndim=8)
        sampler = run_emcee(theta_map, n_walkers=20, n_steps=200,
                            n_burn=50, backend_file="bpb_test.h5")
        print_summary(sampler.get_chain(flat=True))

    # ── bestfit ──────────────────────────────────────────────
    elif mode == "--bestfit":
        theta_map, _ = find_map()
        print_diagnostics(theta_map)
        np.save("bpb_phase1_map.npy", theta_map)
        make_diagnostics_plot(theta_map)

    # ── plot chain ───────────────────────────────────────────
    elif mode == "--plot":
        chainfile = sys.argv[2] if len(sys.argv) > 2 else "bpb_allphases_chain.h5"
        reader  = emcee.backends.HDFBackend(chainfile, read_only=True)
        # Discard first 500 steps as burn-in, thin by 5
        n_iter  = reader.iteration
        discard = min(500, n_iter // 10)
        chain   = reader.get_chain(flat=True, discard=discard, thin=5)
        lp      = reader.get_log_prob(flat=True, discard=discard, thin=5)
        print(f"Chain: {n_iter} steps, {reader.shape[0]} walkers")
        print(f"After discard={discard}, thin=5: {chain.shape[0]} samples")
        print_summary(chain)
        theta_map = chain[np.argmax(lp)]
        print_diagnostics(theta_map)
        make_corner(chain)
        make_diagnostics_plot(theta_map)

    # ── full run ─────────────────────────────────────────────
    else:
        if not HAS_EMCEE:
            print("pip install emcee h5py"); sys.exit(1)

        theta_map, _ = find_map()
        np.save("bpb_phase1_map.npy", theta_map)
        print_diagnostics(theta_map)
        make_diagnostics_plot(theta_map)

        sampler = run_emcee(theta_map, n_walkers=32, n_steps=5000,
                            n_burn=500, backend_file="bpb_allphases_chain.h5")
        chain = sampler.get_chain(flat=True)
        print_summary(chain)
        make_corner(chain)
        lp = sampler.get_log_prob(flat=True)
        np.savez("bpb_phase1_results.npz", chain=chain, lnpost=lp,
                 param_names=np.array(PARAM_NAMES))
        print("Saved: bpb_phase1_results.npz")
