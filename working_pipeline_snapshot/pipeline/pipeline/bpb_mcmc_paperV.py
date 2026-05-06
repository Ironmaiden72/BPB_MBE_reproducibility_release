"""
bpb_mcmc_paperV.py
==================
MCMC for Paper V: joint KiDS-1000 full likelihood + DESI DR1 + Phase 1.

This is a SEPARATE script from bpb_mcmc_v2.py (Paper IV).
Paper IV is frozen. This script adds the full KiDS likelihood block
in place of the compressed S8 Gaussian (Phase 2 of Paper IV).

Differences from Paper IV:
  - Phase 2: compressed S8 Gaussian  →  full KiDS-1000 forward model
  - Phase 3 (ACT sigma8) and Phase 4 (DESI) retained unchanged
  - Phase 5 (Euclid S8) retained unchanged
  - Sampler: emcee, same 32 walkers, 15000 steps (can be adjusted)

Output: bpb_paperV_chain.h5  (emcee backend, reproducible)

Usage:
    python bpb_mcmc_paperV.py [--n_steps 15000] [--n_walkers 32]
                               [--bridge PATH] [--output bpb_paperV_chain.h5]

Dependencies (all in this repo):
    bpb_mbe_unified_relic.py
    bpb_power_spectrum.py
    bpb_kids_pipeline.py
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from typing import Any

import numpy as np

# ── local imports ─────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from bpb_mbe_unified_relic import integrate_growth_mu_const_dynamic
from bpb_kids_pipeline import (
    load_kids_data,
    theory_vector,
    chi2_kids,
    NUISANCE_FIDUCIAL,
)
from bpb_power_spectrum import MAP_PARAMS

# ── emcee ─────────────────────────────────────────────────────────────────
# emcee and h5py are imported lazily in run_mcmc() and summarise_chain()
# so the module is importable for analysis even without them installed.


# ═════════════════════════════════════════════════════════════════════════
# 1.  Parameter vector
# ═════════════════════════════════════════════════════════════════════════

# Paper IV parameter vector (same order, same names):
# θ = {ωm, ωb, H0, a_bg, Δ_bg, κ_bf, β_c, σ8}
PARAM_NAMES = ["omega_m", "omega_b", "H0", "a_bg", "Delta_bg",
               "kappa_bf", "beta_c", "sigma8"]

# Paper IV MAP (starting point)
THETA_MAP = np.array([
    0.1717,   # omega_m
    0.0224,   # omega_b
    70.35,    # H0
    0.268,    # a_bg
    1.042,    # Delta_bg
    2.706,    # kappa_bf
    -0.146,   # beta_c
    0.819,    # sigma8
])

# Prior bounds (same as Paper IV)
PRIOR_LOW  = np.array([0.05,  0.015, 55.0, 0.05, 0.1, 0.0, -0.5, 0.6])
PRIOR_HIGH = np.array([0.50,  0.035, 90.0, 0.95, 2.5, 8.0,  0.5, 1.1])

# BBN Gaussian prior on omega_b
BBN_MEAN  = 0.02237
BBN_SIGMA = 0.00015


def theta_to_params(theta: np.ndarray) -> dict:
    """Convert MCMC parameter vector to params dict for BPB pipeline."""
    om, ob, H0, a_bg, D_bg, k_bf, beta_c, s8 = theta
    h = H0 / 100.0
    return dict(
        H0       = H0,
        Om_eff0  = om / h**2,          # Omega_m from omega_m = Om*h^2
        Ob0      = ob / h**2,
        ns       = 0.965,              # fixed (same as Paper IV)
        sigma8_0 = s8,
        a_bg     = a_bg,
        Delta_bg = D_bg,
        kappa_bf = k_bf,
        beta_c   = beta_c,
        Or0      = 0.0,
        Ok0      = 0.0,
    )


# ═════════════════════════════════════════════════════════════════════════
# 2.  Likelihood blocks
# ═════════════════════════════════════════════════════════════════════════

# ── Block A: Phase 1 (BAO + CC + SNe + fσ8) ──────────────────────────────
# Re-implemented here as a standalone function using the same data as Paper IV
# to keep this script fully self-contained.

# fσ8 data from Paper IV Table 4 (22 points)
_FSig8_DATA = np.array([
    (0.067, 0.423, 0.055), (0.170, 0.510, 0.060), (0.220, 0.420, 0.070),
    (0.250, 0.351, 0.058), (0.320, 0.384, 0.095), (0.350, 0.440, 0.050),
    (0.370, 0.460, 0.038), (0.380, 0.497, 0.045), (0.410, 0.450, 0.040),
    (0.510, 0.458, 0.038), (0.570, 0.441, 0.043), (0.600, 0.430, 0.040),
    (0.600, 0.550, 0.120), (0.610, 0.436, 0.034), (0.700, 0.473, 0.041),
    (0.770, 0.490, 0.180), (0.780, 0.380, 0.040), (0.850, 0.315, 0.095),
    (0.860, 0.400, 0.110), (1.400, 0.482, 0.116), (1.480, 0.462, 0.045),
    (1.520, 0.420, 0.076),
])
_FSig8_Z   = _FSig8_DATA[:, 0]
_FSig8_OBS = _FSig8_DATA[:, 1]
_FSig8_ERR = _FSig8_DATA[:, 2]

# DESI DR1 fσ8 data from Paper IV Table 5
_DESI_DATA = np.array([
    (0.295, 0.460, 0.023), (0.510, 0.455, 0.018), (0.706, 0.430, 0.014),
    (0.934, 0.395, 0.016), (1.321, 0.355, 0.020), (1.484, 0.342, 0.027),
    (2.330, 0.282, 0.031),
])
_DESI_Z   = _DESI_DATA[:, 0]
_DESI_OBS = _DESI_DATA[:, 1]
_DESI_ERR = _DESI_DATA[:, 2]


def _lnL_fsig8(params: dict, z_obs: np.ndarray, fsig8_obs: np.ndarray,
               fsig8_err: np.ndarray) -> float:
    """Gaussian log-likelihood for fσ8 measurements."""
    try:
        z_tab, D_tab, f_tab = integrate_growth_mu_const_dynamic(
            Om_eff0  = params["Om_eff0"],
            kappa_bf = params["kappa_bf"],
            beta_c   = params["beta_c"],
            H0       = params["H0"],
            a_bg     = params["a_bg"],
            Delta_bg = params["Delta_bg"],
            n_a      = 800,
        )
    except Exception:
        return -np.inf

    # σ8(z) = σ8(0) * D(z)
    f_interp    = np.interp(z_obs, z_tab[::-1], f_tab[::-1])
    D_interp    = np.interp(z_obs, z_tab[::-1], D_tab[::-1])
    fsig8_theory = f_interp * D_interp * params["sigma8_0"]

    chi2 = np.sum(((fsig8_obs - fsig8_theory) / fsig8_err)**2)
    return -0.5 * chi2


def lnL_phase1(params: dict) -> float:
    """Phase 1: fσ8 only (BAO/CC/SNe not re-implemented here;
    use Paper IV bpb_mcmc_v2.py for a full Phase 1 chain).
    For Paper V the dominant growth constraint is from DESI + KiDS full;
    Phase 1 fσ8 provides secondary anchoring."""
    return _lnL_fsig8(params, _FSig8_Z, _FSig8_OBS, _FSig8_ERR)


def lnL_phase3_ACT(params: dict) -> float:
    """Phase 3: ACT DR6 + Planck CMB lensing (compressed σ8 constraint)."""
    # Same as Paper IV: Gaussian on sigma8 = 0.811 ± 0.013
    s8_ACT   = 0.811
    sig_ACT  = 0.013
    s8_model = params["sigma8_0"]   # approximation: σ8(0) ≈ σ8 at CMB lensing z
    return -0.5 * ((s8_model - s8_ACT) / sig_ACT)**2


def lnL_phase4_DESI(params: dict) -> float:
    """Phase 4: DESI DR1 fσ8 (7 bins)."""
    return _lnL_fsig8(params, _DESI_Z, _DESI_OBS, _DESI_ERR)


def lnL_phase5_Euclid(params: dict) -> float:
    """Phase 5: Euclid Year 1 S8 (compressed)."""
    # Gaussian on S8 = 0.776 ± 0.017
    Om0 = params["Om_eff0"]
    S8_model = params["sigma8_0"] * np.sqrt(Om0 / 0.3)
    return -0.5 * ((S8_model - 0.776) / 0.017)**2


# ── Block B: Phase 2 replacement — full KiDS-1000 ─────────────────────────
# kids_data is loaded once and cached globally for speed

_KIDS_DATA: dict | None = None
_KIDS_BRIDGE_PATH: str  = ""


def _get_kids_data(bridge_path: str) -> dict:
    global _KIDS_DATA, _KIDS_BRIDGE_PATH
    if _KIDS_DATA is None or bridge_path != _KIDS_BRIDGE_PATH:
        print(f"Loading KiDS-1000 data from {bridge_path}...")
        _KIDS_DATA = load_kids_data(bridge_path)
        _KIDS_BRIDGE_PATH = bridge_path
        print(f"  {len(_KIDS_DATA['full_vec'])} data points loaded.")
    return _KIDS_DATA


def lnL_phase2_KiDS_full(
    params: dict,
    bridge_path: str,
    ell_grid: np.ndarray,
    n_chi: int,
) -> float:
    """Phase 2 replacement: full KiDS-1000 forward model."""
    kids = _get_kids_data(bridge_path)
    try:
        xip_t, xim_t = theory_vector(
            kids, params=params,
            nuisance=NUISANCE_FIDUCIAL,
            ell_grid=ell_grid,
            n_chi=n_chi,
            verbose=False,
        )
        c2 = chi2_kids(xip_t, xim_t, kids)
        return -0.5 * c2
    except Exception as e:
        warnings.warn(f"KiDS theory failed: {e}")
        return -np.inf


# ═════════════════════════════════════════════════════════════════════════
# 3.  Prior and posterior
# ═════════════════════════════════════════════════════════════════════════

def lnprior(theta: np.ndarray) -> float:
    """Uniform priors with BBN Gaussian on omega_b."""
    if np.any(theta < PRIOR_LOW) or np.any(theta > PRIOR_HIGH):
        return -np.inf
    om, ob, H0, a_bg, D_bg, k_bf, beta_c, s8 = theta
    # BBN prior
    lnp = -0.5 * ((ob - BBN_MEAN) / BBN_SIGMA)**2
    return lnp


def lnposterior(
    theta: np.ndarray,
    bridge_path: str,
    ell_grid: np.ndarray,
    n_chi: int,
) -> float:
    """Full log-posterior: prior + all likelihood blocks."""
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf

    params = theta_to_params(theta)

    # Accumulate log-likelihoods
    lnL = lp

    # Phase 1: fσ8
    lnL += lnL_phase1(params)
    if not np.isfinite(lnL):
        return -np.inf

    # Phase 2: full KiDS-1000 (replaces compressed S8)
    lnL += lnL_phase2_KiDS_full(params, bridge_path, ell_grid, n_chi)
    if not np.isfinite(lnL):
        return -np.inf

    # Phase 3: ACT + Planck lensing
    lnL += lnL_phase3_ACT(params)
    if not np.isfinite(lnL):
        return -np.inf

    # Phase 4: DESI DR1
    lnL += lnL_phase4_DESI(params)
    if not np.isfinite(lnL):
        return -np.inf

    # Phase 5: Euclid Year 1
    lnL += lnL_phase5_Euclid(params)

    return lnL


# ═════════════════════════════════════════════════════════════════════════
# 4.  Main MCMC runner
# ═════════════════════════════════════════════════════════════════════════

def run_mcmc(
    bridge_path: str,
    output_path: str = "bpb_paperV_chain.h5",
    n_steps: int = 15000,
    n_walkers: int = 32,
    n_burn: int = 1500,
    ell_n: int = 200,
    n_chi: int = 150,
    resume: bool = True,
    verbose: bool = True,
) -> None:
    """
    Run the Paper V MCMC and save to HDF5 backend.

    Parameters
    ----------
    bridge_path : str
        Path to kids_xipm_bridge.npz
    output_path : str
        Output chain file (emcee HDF5 backend)
    n_steps : int
        Total production steps (excluding burn-in)
    n_walkers : int
        Number of emcee walkers
    n_burn : int
        Burn-in steps (discarded before thinning)
    ell_n : int
        Number of ell points in Limber grid (200 = fast, 300 = accurate)
    n_chi : int
        Number of chi points in Limber integration (150 = fast, 250 = accurate)
    resume : bool
        If True and output_path exists, resume from last position
    """
    ell_grid = np.geomspace(2.0, 5e4, ell_n)
    ndim     = len(THETA_MAP)

    if verbose:
        print("=" * 60)
        print("BPB/MBE Paper V MCMC")
        print(f"  n_walkers = {n_walkers}")
        print(f"  n_steps   = {n_steps}")
        print(f"  ell_n     = {ell_n}, n_chi = {n_chi}")
        print(f"  output    = {output_path}")
        print("=" * 60)

    # Pre-load KiDS data
    _get_kids_data(bridge_path)

    # Warm-up timing estimate
    if verbose:
        print("\nTiming single evaluation at MAP...")
        t0 = time.time()
        lnposterior(THETA_MAP, bridge_path, ell_grid, n_chi)
        t_eval = time.time() - t0
        t_total = t_eval * n_steps * n_walkers
        print(f"  Single eval: {t_eval:.1f}s")
        print(f"  Estimated total: {t_total/3600:.1f}h "
              f"({t_total/3600/24:.1f}d) — run on local machine")
        print()

    # emcee backend (HDF5, reproducible)
    backend = emcee.backends.HDFBackend(output_path)

    if resume and os.path.exists(output_path):
        try:
            n_done = backend.iteration
            if verbose:
                print(f"Resuming from step {n_done}...")
        except Exception:
            backend.reset(n_walkers, ndim)
            n_done = 0
    else:
        backend.reset(n_walkers, ndim)
        n_done = 0

    # Initial positions: Gaussian ball around MAP
    if n_done == 0:
        rng    = np.random.default_rng(42)
        scales = np.abs(THETA_MAP) * 0.02  # 2% scatter
        scales = np.clip(scales, 1e-4, None)
        p0 = THETA_MAP + rng.normal(0, 1, (n_walkers, ndim)) * scales
        # Clip to prior
        p0 = np.clip(p0, PRIOR_LOW + 1e-6, PRIOR_HIGH - 1e-6)
    else:
        p0 = backend.get_last_sample().coords

    # Wrap posterior for emcee
    def log_prob(theta: np.ndarray) -> float:
        return lnposterior(theta, bridge_path, ell_grid, n_chi)

    sampler = emcee.EnsembleSampler(
        n_walkers, ndim, log_prob,
        backend=backend,
    )

    n_remaining = max(0, n_steps - n_done)
    if n_remaining == 0:
        if verbose: print("Chain already complete.")
        return

    if verbose:
        print(f"Running {n_remaining} steps "
              f"({'resuming' if n_done>0 else 'fresh start'})...")

    t0 = time.time()
    for sample in sampler.sample(p0, iterations=n_remaining, progress=verbose):
        pass

    if verbose:
        t_wall = time.time() - t0
        print(f"\nDone. Wall time: {t_wall/3600:.2f}h")
        try:
            tau = sampler.get_autocorr_time(quiet=True)
            print(f"Autocorrelation time: {tau}")
            print(f"N/tau_max = {n_steps / tau.max():.1f}")
        except Exception as e:
            print(f"(autocorr not yet reliable: {e})")
        print(f"Chain saved to: {output_path}")


# ═════════════════════════════════════════════════════════════════════════
# 5.  Analysis utilities (run after chain is done)
# ═════════════════════════════════════════════════════════════════════════

def summarise_chain(
    chain_path: str,
    n_burn: int = 1500,
    thin: int = 3,
) -> dict:
    """
    Load chain, apply burn-in + thinning, print marginalised posteriors.
    Returns flat samples array and parameter summary.
    """
    backend = emcee.backends.HDFBackend(chain_path, read_only=True)
    n_done  = backend.iteration
    print(f"Chain: {n_done} steps, {backend.shape[0]} walkers")

    samples = backend.get_chain(discard=n_burn, thin=thin, flat=True)
    print(f"Flat samples: {len(samples)} (after discard={n_burn}, thin={thin})")

    print("\nMarginalised posteriors (median ± 68% CI):")
    summary = {}
    for i, name in enumerate(PARAM_NAMES):
        v  = samples[:, i]
        med = np.median(v)
        lo  = med - np.percentile(v, 16)
        hi  = np.percentile(v, 84) - med
        print(f"  {name:12s}: {med:+.4f}  +{hi:.4f}  -{lo:.4f}")
        summary[name] = dict(median=med, upper=hi, lower=lo)

    # Key result: beta_c
    bc = samples[:, 6]
    p_pos = (bc >= 0).mean()
    print(f"\n  P(beta_c >= 0) = {p_pos:.2e}")
    print(f"  beta_c < 0: {'yes' if p_pos < 0.05 else 'marginal'}")

    return dict(samples=samples, summary=summary, p_betac_pos=p_pos)


# ═════════════════════════════════════════════════════════════════════════
# 6.  CLI
# ═════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BPB/MBE Paper V MCMC — full KiDS-1000 likelihood"
    )
    parser.add_argument("--bridge",
        default="real_data/kids_real_bridge/kids_xipm_bridge.npz",
        help="Path to kids_xipm_bridge.npz")
    parser.add_argument("--output",
        default="bpb_paperV_chain.h5",
        help="Output chain file")
    parser.add_argument("--n_steps",  type=int, default=15000)
    parser.add_argument("--n_walkers",type=int, default=32)
    parser.add_argument("--n_burn",   type=int, default=1500)
    parser.add_argument("--ell_n",    type=int, default=200,
        help="ell grid points (200=fast, 300=accurate)")
    parser.add_argument("--n_chi",    type=int, default=150,
        help="chi grid points (150=fast, 250=accurate)")
    parser.add_argument("--no_resume", action="store_true")
    parser.add_argument("--summarise", action="store_true",
        help="Summarise existing chain instead of running")
    args = parser.parse_args()

    if args.summarise:
        summarise_chain(args.output, n_burn=args.n_burn)
    else:
        run_mcmc(
            bridge_path = args.bridge,
            output_path = args.output,
            n_steps     = args.n_steps,
            n_walkers   = args.n_walkers,
            n_burn      = args.n_burn,
            ell_n       = args.ell_n,
            n_chi       = args.n_chi,
            resume      = not args.no_resume,
        )
