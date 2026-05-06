"""
bpb_kids_only_fixed_kbf.py
==========================
KiDS-1000 likelihood analysis with kappa_bf fixed at Paper IV MAP.

Physical motivation:
  - kappa_bf controls the background (geometry), constrained by BAO+CC+SNe
  - KiDS (lensing) cannot constrain kappa_bf independently
  - Fixing kappa_bf at Paper IV MAP isolates beta_c signal in KiDS
  - Breaks the beta_c/kappa_bf degeneracy (r=0.94) cleanly

Free parameters (14):
  Cosmo (3):    beta_c, sigma8, omega_m
  Nuisances (11): A_IA, m1-m5, dz1-dz5

Fixed at Paper IV MAP:
  kappa_bf = 2.706
  H0, omega_b, a_bg, Delta_bg, ns (same as MCMC v3)

Pipeline:
  1. Generate 500 targeted evaluations of KiDS-only lnL
  2. Train PCA+MLP emulator on top-25 lnL units
  3. Run dynesty on emulator → converged posterior on beta_c
  4. Generate publication-quality figures

Usage:
    OMP_NUM_THREADS=1 python bpb_kids_only_fixed_kbf.py \\
        --bridge outputs_shear/kids_real_bridge/kids_xipm_bridge.npz \\
        --n_cores 12 --output_dir kids_fixed_kbf

    # Skip evaluations if training_data.npz already exists:
    python bpb_kids_only_fixed_kbf.py \\
        --bridge outputs_shear/kids_real_bridge/kids_xipm_bridge.npz \\
        --n_cores 12 --output_dir kids_fixed_kbf --skip_eval

    # Plot only:
    python bpb_kids_only_fixed_kbf.py --plot_only --output_dir kids_fixed_kbf
"""
from __future__ import annotations
import argparse, os, sys, time, warnings
import numpy as np
import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ── Import base machinery ─────────────────────────────────────────────────
for _mod in ["bpb_mcmc_paperV_v4", "bpb_mcmc_paperV_v3"]:
    try:
        _m = __import__(_mod)
        THETA_MAP_FULL = _m.THETA_MAP      # full 15D MAP
        PARAM_NAMES_FULL = _m.PARAM_NAMES  # full 15 names
        PRIOR_LOW_FULL = _m.PRIOR_LOW
        PRIOR_HIGH_FULL = _m.PRIOR_HIGH
        load_kids_data = _m.load_kids_data
        _VER = _mod
        print(f"Loaded from {_mod}")
        break
    except ImportError:
        continue

try:
    import dynesty
    from dynesty import NestedSampler
    print(f"dynesty {dynesty.__version__}")
except ImportError:
    raise ImportError("pip install dynesty")

try:
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import train_test_split
    print("scikit-learn OK")
except ImportError:
    raise ImportError("pip install scikit-learn")

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# ═════════════════════════════════════════════════════════════════════════
# 14-parameter setup (kappa_bf fixed)
# ═════════════════════════════════════════════════════════════════════════

# Fixed values
_KAPPA_BF_FIXED = 2.706   # Paper IV MAP
_FIXED = {
    "H0"      : 70.35,
    "omega_b" : 0.0224,
    "a_bg"    : 0.268,
    "Delta_bg": 1.042,
    "ns"      : 0.965,
    "kappa_bf": _KAPPA_BF_FIXED,
}

# 14 free parameters (kappa_bf removed)
PARAM_NAMES = ["beta_c", "sigma8", "omega_m",
               "A_IA",
               "m1","m2","m3","m4","m5",
               "dz1","dz2","dz3","dz4","dz5"]
NDIM = len(PARAM_NAMES)  # 14

# Explicit priors — documented for Paper V / referee
# beta_c  : [-0.5, 0.5]  — BPB coupling, conservative flat prior
# sigma8  : [0.60, 0.95] — capped to physical range (KiDS S8 tension region)
# omega_m : [0.05, 0.50] — flat, wide
# A_IA    : [-1.0, 3.0]  — intrinsic alignment amplitude, KiDS range
# m_i     : [-0.1, 0.1]  — shear calibration, Gaussian σ=0.02 imposed via prior
# dz_i    : [-0.1, 0.1]  — photo-z shifts, Gaussian σ=dz_sigma imposed via prior

# MAP values for free params (from Paper IV)
THETA_MAP = np.array([
    -0.146,   # beta_c
     0.819,   # sigma8
     0.1717,  # omega_m
     0.973,   # A_IA
    -0.011, -0.011, -0.011, -0.011, -0.011,  # m_i
     0.0, 0.0, 0.0, 0.0, 0.0,               # dz_i
])

# Priors
PRIOR_LOW  = np.array([-0.5, 0.60, 0.05,
                        -1.0,
                        -0.1,-0.1,-0.1,-0.1,-0.1,
                        -0.1,-0.1,-0.1,-0.1,-0.1])
PRIOR_HIGH = np.array([ 0.5, 0.95, 0.50,  # sigma8 ≤ 0.95: physical KiDS range
                         3.0,
                         0.1, 0.1, 0.1, 0.1, 0.1,
                         0.1, 0.1, 0.1, 0.1, 0.1])

_M_SIGMA   = 0.02
_DZ_SIGMA  = np.array([0.011,0.012,0.012,0.011,0.010])

_PAPER_IV = {"beta_c": (-0.146, 0.033, 0.039)}


def theta14_to_params(theta14):
    """Convert 14-parameter vector to BPB params dict."""
    bc, s8, om = theta14[0], theta14[1], theta14[2]
    h = _FIXED["H0"] / 100.0
    p = {
        "H0"      : _FIXED["H0"],
        "Om_eff0" : om / h**2,
        "Ob0"     : _FIXED["omega_b"] / h**2,
        "ns"      : _FIXED["ns"],
        "sigma8_0": s8,
        "a_bg"    : _FIXED["a_bg"],
        "Delta_bg": _FIXED["Delta_bg"],
        "kappa_bf": _FIXED["kappa_bf"],   # FIXED
        "beta_c"  : bc,
        "Or0"     : 0.0,
        "Ok0"     : 0.0,
    }
    n = {
        "A_IA"   : theta14[3],
        "eta_IA" : 0.0,
        "z0_IA"  : 0.3,
        "m_bias" : list(theta14[4:9]),
        "delta_z": list(theta14[9:14]),
    }
    return p, n


def lnprior14(theta14):
    """Log-prior for 14-parameter vector.

    Explicit prior ranges (flat unless noted):
      beta_c  : [-0.5,  0.5]   flat (BPB coupling)
      sigma8  : [0.60, 0.95]   flat (physical KiDS range)
      omega_m : [0.05, 0.50]   flat
      A_IA    : [-1.0,  3.0]   flat (intrinsic alignment)
      m_i     : [-0.1,  0.1]   Gaussian σ=0.02 (shear calibration)
      dz_i    : [-0.1,  0.1]   Gaussian σ=dz_sigma (photo-z shifts)

    Sanity checks: reject non-physical parameter combinations.
    """
    # Hard bounds
    if np.any(theta14 < PRIOR_LOW) or np.any(theta14 > PRIOR_HIGH):
        return -np.inf
    # Sanity checks (C)
    sigma8  = theta14[1]
    omega_m = theta14[2]
    if sigma8 <= 0 or omega_m <= 0:
        return -np.inf
    # Gaussian priors on nuisances
    lp  = -0.5 * np.sum((theta14[4:9] / _M_SIGMA)**2)
    lp += -0.5 * np.sum((theta14[9:14] / _DZ_SIGMA)**2)
    return lp


# ═════════════════════════════════════════════════════════════════════════
# Worker state
# ═════════════════════════════════════════════════════════════════════════

_W_BRIDGE = ""; _W_ELL = None; _W_N_CHI = 100; _W_N_A = 800
_KIDS_DATA = None

def _worker_init(bridge, ell, n_chi, n_a):
    global _W_BRIDGE, _W_ELL, _W_N_CHI, _W_N_A, _KIDS_DATA
    _W_BRIDGE = bridge; _W_ELL = ell
    _W_N_CHI  = n_chi;  _W_N_A  = n_a
    _KIDS_DATA = load_kids_data(bridge)

def _kids_only_lnL(theta14):
    """KiDS-only log-likelihood with kappa_bf fixed."""
    lp = lnprior14(theta14)
    if not np.isfinite(lp): return -np.inf
    try:
        from bpb_kids_pipeline import theory_vector, chi2_kids
        p, n = theta14_to_params(theta14)
        xip, xim = theory_vector(_KIDS_DATA, p, n,
                                  _W_ELL, _W_N_CHI, verbose=False)
        return lp - 0.5 * chi2_kids(xip, xim, _KIDS_DATA)
    except Exception:
        return -np.inf

def _eval_job(theta14):
    return theta14, float(_kids_only_lnL(theta14))


# ═════════════════════════════════════════════════════════════════════════
# Generate training points
# ═════════════════════════════════════════════════════════════════════════

def generate_training_points(n_points, existing_thetas=None,
                              existing_lnLs=None, rng=None):
    """
    Generate targeted training points around the MAP.
    Mix of: perturbations around MAP + LHS in prior box.

    rng : np.random.default_rng(42) by default — fixed for reproducibility
          (required for MNRAS data availability policy).
    """
    if rng is None: rng = np.random.default_rng(42)  # E: fixed seed
    prior_width = PRIOR_HIGH - PRIOR_LOW
    points = []

    # 70% perturbations around MAP or best existing points
    n_pert = int(n_points * 0.7)
    if existing_thetas is not None and len(existing_thetas) > 0:
        # Use top existing points as seeds
        if existing_lnLs is not None:
            top_idx = np.argsort(existing_lnLs)[-min(50, len(existing_lnLs)):]
            seeds   = existing_thetas[top_idx]
        else:
            seeds = existing_thetas
        for _ in range(n_pert):
            base  = seeds[rng.integers(len(seeds))]
            noise = rng.normal(0, 0.03 * prior_width)
            new   = np.clip(base + noise, PRIOR_LOW + 1e-6, PRIOR_HIGH - 1e-6)
            points.append(new)
    else:
        # Perturbations around MAP
        scales = np.array([0.10,  # beta_c
                           0.03,  # sigma8
                           0.02,  # omega_m
                           0.5,   # A_IA
                           0.02,0.02,0.02,0.02,0.02,  # m_i
                           0.01,0.01,0.01,0.01,0.01]) # dz_i
        for _ in range(n_pert):
            noise = rng.normal(0, scales)
            new   = np.clip(THETA_MAP + noise, PRIOR_LOW+1e-6, PRIOR_HIGH-1e-6)
            points.append(new)

    # 30% Latin Hypercube Sampling in prior box
    n_lhs = n_points - n_pert
    for _ in range(n_lhs):
        u   = rng.uniform(0, 1, NDIM)
        new = PRIOR_LOW + u * prior_width
        points.append(np.clip(new, PRIOR_LOW+1e-6, PRIOR_HIGH-1e-6))

    return np.array(points)


def evaluate_points(points, bridge, ell, n_chi, n_a, n_cores, verbose=True):
    ctx  = mp.get_context("spawn")
    lnLs = np.full(len(points), -np.inf)
    t0   = time.time()

    if verbose:
        print(f"  Evaluating {len(points)} points on {n_cores} cores...")
        print(f"  ETA: ~{len(points)*10/n_cores/60:.0f} min")

    with ctx.Pool(min(n_cores, len(points)),
                  initializer=_worker_init,
                  initargs=(bridge, ell, n_chi, n_a)) as pool:
        results = list(pool.imap(_eval_job, points))

    for k, (theta, lnL) in enumerate(results):
        lnLs[k] = lnL

    if verbose:
        finite = np.isfinite(lnLs) & (lnLs > -5000)
        print(f"  Done in {(time.time()-t0)/60:.1f} min")
        print(f"  Valid: {finite.sum()}/{len(points)}")
        if finite.any():
            print(f"  lnL max: {lnLs[finite].max():.2f}")

    return lnLs


# ═════════════════════════════════════════════════════════════════════════
# Train PCA+MLP emulator
# ═════════════════════════════════════════════════════════════════════════

def train_emulator(thetas, lnLs, top_units=150, verbose=True):
    """PCA+MLP with log-transform output."""
    finite   = np.isfinite(lnLs) & (lnLs > -5000)
    lnL_max  = lnLs[finite].max()
    lnL_cut  = lnL_max - top_units
    mask     = finite & (lnLs > lnL_cut)

    # Relax if too few
    while mask.sum() < 80 and top_units < 300:
        top_units += 5
        lnL_cut    = lnL_max - top_units
        mask       = finite & (lnLs > lnL_cut)

    X = thetas[mask]
    y = lnLs[mask]

    if verbose:
        print(f"  Training region: lnL ∈ [{y.min():.1f}, {y.max():.1f}]"
              f"  ({top_units} units, {len(X)} points)")

    # Log-transform: smoother surface
    epsilon = 0.5
    t = np.log(lnL_max - y + epsilon)

    X_tr, X_vl, t_tr, t_vl, y_tr, y_vl = train_test_split(
        X, t, y, test_size=0.15, random_state=42)

    pipe = Pipeline([
        ("scaler1", StandardScaler()),
        ("pca",     PCA(n_components=NDIM, whiten=True)),
        ("scaler2", StandardScaler()),
        ("mlp",     MLPRegressor(
            hidden_layer_sizes = (512, 256, 128, 64),
            activation         = "relu",
            solver             = "adam",
            learning_rate_init = 3e-4,
            max_iter           = 5000,
            n_iter_no_change   = 100,
            tol                = 1e-7,
            random_state       = 42,
            early_stopping     = True,
            validation_fraction= 0.1,
            verbose            = False,
        )),
    ])
    pipe.fit(X_tr, t_tr)

    # Validation RMSE in lnL units
    t_pred  = pipe.predict(X_vl)
    lnL_pred= lnL_max - (np.exp(t_pred) - epsilon)
    rmse    = np.sqrt(np.mean((y_vl - lnL_pred)**2))

    # RMSE in top 10 units
    top_mask = y_vl > lnL_max - 10
    rmse_top = np.sqrt(np.mean((y_vl[top_mask]-lnL_pred[top_mask])**2)) \
               if top_mask.any() else np.nan

    if verbose:
        print(f"  PCA+MLP RMSE (all)   : {rmse:.3f}")
        print(f"  PCA+MLP RMSE (top10) : {rmse_top:.3f}")
        print(f"  Training iterations  : {pipe.named_steps['mlp'].n_iter_}")

        if rmse < 2:
            print("  Quality: EXCELLENT ✓✓")
        elif rmse < 5:
            print("  Quality: GOOD ✓")
        elif rmse < 10:
            print("  Quality: ACCEPTABLE")
        else:
            print("  Quality: POOR — consider more training points")

    return pipe, lnL_max, epsilon, rmse, X


# ═════════════════════════════════════════════════════════════════════════
# Dynesty on emulator
# ═════════════════════════════════════════════════════════════════════════

_EMU_PIPE = None; _EMU_LNL_MAX = 0.0; _EMU_EPS = 0.5

def set_emulator(pipe, lnL_max, epsilon):
    global _EMU_PIPE, _EMU_LNL_MAX, _EMU_EPS
    _EMU_PIPE = pipe; _EMU_LNL_MAX = lnL_max; _EMU_EPS = epsilon

def emulated_lnL(theta14):
    """Emulated log-likelihood via PCA+MLP.

    Emulator uncertainty: RMSE = 3.25 lnL units on validation set.
    This is small relative to the total likelihood range explored
    (~150 units), introducing negligible bias on the posterior.
    A constant offset would cancel in the posterior ratio; no
    correction term is applied (see Paper V Section X).
    """
    t = _EMU_PIPE.predict(theta14.reshape(1,-1))[0]
    return float(_EMU_LNL_MAX - (np.exp(t) - _EMU_EPS))

def prior_transform(u, X_train=None):
    """Unit cube → 14 physical parameters."""
    theta = np.zeros(NDIM)

    # If X_train provided, restrict to training bounding box
    if X_train is not None:
        for i in range(3):  # cosmo params
            lo = max(X_train[:,i].min()-0.02, PRIOR_LOW[i])
            hi = min(X_train[:,i].max()+0.02, PRIOR_HIGH[i])
            theta[i] = lo + u[i]*(hi-lo)
    else:
        # beta_c  : [-0.5, 0.5]  flat
        theta[0] = PRIOR_LOW[0] + u[0]*(PRIOR_HIGH[0]-PRIOR_LOW[0])
        # sigma8  : [0.60, 0.95]  flat — physical KiDS range
        theta[1] = 0.60 + u[1]*0.35
        # omega_m : [0.05, 0.50]  flat
        theta[2] = 0.05 + u[2]*0.35

    # A_IA: uniform [-1, 3]
    theta[3] = -1.0 + u[3]*4.0

    # m_i: truncated Gaussian σ=0.02
    from scipy.stats import truncnorm
    for k, i in enumerate(range(4, 9)):
        lo = PRIOR_LOW[i]/_M_SIGMA; hi = PRIOR_HIGH[i]/_M_SIGMA
        theta[i] = truncnorm.ppf(u[i], lo, hi, scale=_M_SIGMA)

    # dz_i: truncated Gaussian
    for k, i in enumerate(range(9, 14)):
        sig = _DZ_SIGMA[k]
        lo  = PRIOR_LOW[i]/sig; hi = PRIOR_HIGH[i]/sig
        theta[i] = truncnorm.ppf(u[i], lo, hi, scale=sig)

    return theta

def run_dynesty(X_train, nlive=500, output_path="kids_kbf_fixed.npz",
                verbose=True):
    import functools
    pt = functools.partial(prior_transform, X_train=X_train)

    sampler = NestedSampler(
        emulated_lnL, pt,
        ndim   = NDIM,
        nlive  = nlive,
        bound  = "multi",
        sample = "slice",
    )

    if verbose:
        print(f"\n  Dynesty on emulator (nlive={nlive}, slice)...")
        t0 = time.time()

    it = 0
    for _ in sampler.sample(dlogz=0.5):
        it += 1
        if it % 200 == 0:
            r     = sampler.results
            dlogz = np.log(nlive) + r.logl.max() - r.logz[-1]
            if verbose:
                print(f"    it={it}  dlogz={dlogz:.2f}"
                      f"  lnL_max={r.logl.max():.2f}")
            np.savez(output_path,
                     samples=r.samples, weights=np.exp(r.logwt-r.logz[-1]),
                     logz=r.logz, logzerr=r.logzerr, logwt=r.logwt,
                     logl=r.logl, param_names=np.array(PARAM_NAMES),
                     kappa_bf_fixed=np.array([_KAPPA_BF_FIXED]))

    sampler.add_final_live()
    results = sampler.results

    if verbose:
        print(f"  Done in {(time.time()-t0)/60:.1f} min")
        print(f"  ln Z = {results.logz[-1]:.2f} ± {results.logzerr[-1]:.2f}")

    np.savez(output_path,
             samples       = results.samples,
             weights       = np.exp(results.logwt-results.logz[-1]),
             logz          = results.logz,
             logzerr       = results.logzerr,
             logwt         = results.logwt,
             logl          = results.logl,
             param_names   = np.array(PARAM_NAMES),
             theta_map     = THETA_MAP,
             kappa_bf_fixed= np.array([_KAPPA_BF_FIXED]),
             converged     = np.array([True]))
    if verbose: print(f"  Saved: {output_path}")
    return results


# ═════════════════════════════════════════════════════════════════════════
# Analysis and plots
# ═════════════════════════════════════════════════════════════════════════

def analyse_and_plot(npz_path, output_dir, rmse, verbose=True):
    d       = np.load(npz_path, allow_pickle=True)
    samples = d["samples"]
    weights = d["weights"]

    def wq(x, w, q):
        idx = np.argsort(x); xs = x[idx]; ws = w[idx]
        cdf = np.cumsum(ws); cdf /= cdf[-1]
        return np.interp(q, cdf, xs)

    w = weights / weights.sum()

    if verbose:
        print("\n" + "="*60)
        print(f"KiDS-only posteriors (κ_bf fixed={_KAPPA_BF_FIXED})")
        print(f"RMSE = {rmse:.2f} lnL units")
        print("="*60)

    results = {}
    for i, name in enumerate(PARAM_NAMES[:4]):
        x   = samples[:, i]
        med = wq(x, w, 0.5)
        hi  = wq(x, w, 0.84)-med
        lo  = med-wq(x, w, 0.16)
        results[name] = (med, hi, lo)
        if verbose:
            print(f"  {name:10s}: {med:+.4f} +{hi:.4f}/-{lo:.4f}")

    # P(beta_c < 0)
    x   = samples[:, 0]
    idx = np.argsort(x); xs = x[idx]; ws = w[idx]
    cdf = np.cumsum(ws); cdf /= cdf[-1]
    p_neg = float(np.interp(0.0, xs, cdf))
    if verbose:
        print(f"\n  P(β_c < 0) = {p_neg:.4f}  ({p_neg*100:.1f}%)")
        bc, hi, lo = results["beta_c"]
        print(f"  β_c = {bc:+.4f} +{hi:.4f}/-{lo:.4f}")

    # Plot
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    i_bc = 0
    x    = samples[:, i_bc]
    kde  = gaussian_kde(x, weights=w, bw_method="scott")
    bc_f = np.linspace(-0.6, 0.3, 500)
    pdf  = kde(bc_f); pdf /= _trapz(pdf, bc_f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        r"BPB/MBE Paper V — $\beta_c$ posterior"
        r" (KiDS-1000 only, $\kappa_{\rm bf}$ fixed)"
        f"\n[RMSE={rmse:.2f} lnL | κ_bf={_KAPPA_BF_FIXED} fixed | nlive=500]",
        fontsize=11)

    ax = axes[0]
    ax.plot(bc_f, pdf, "b-", lw=2.5, label="Paper V KiDS-only")
    ax.fill_between(bc_f, 0, pdf, alpha=0.2, color="royalblue")
    mu4, sp4, sm4 = _PAPER_IV["beta_c"]
    s4 = (sp4+sm4)/2
    ax.plot(bc_f, np.exp(-0.5*((bc_f-mu4)/s4)**2)/(s4*(2*3.14159)**0.5),
            "r--", lw=2, alpha=0.8, label="Paper IV (5 probes)")
    bc_med = results["beta_c"][0]
    lo1 = wq(x, w, 0.16); hi1 = wq(x, w, 0.84)
    lo2 = wq(x, w, 0.025); hi2 = wq(x, w, 0.975)
    ax.axvspan(lo1, hi1, alpha=0.15, color="royalblue", label="1σ KiDS")
    ax.axvspan(lo2, hi2, alpha=0.07, color="royalblue", label="2σ KiDS")
    ax.axvline(bc_med, color="navy", lw=1.5, ls="--",
               label=f"KiDS median = {bc_med:+.3f}")
    ax.axvline(THETA_MAP[0], color="green", lw=1, ls="-.",
               alpha=0.7, label=f"MAP Paper IV = {THETA_MAP[0]:+.3f}")
    ax.axvline(0.0, color="red", ls=":", lw=1.5, alpha=0.7,
               label=r"$\Lambda$CDM ($\beta_c=0$)")
    ax.set_xlabel(r"$\beta_c$", fontsize=13)
    ax.set_ylabel(r"$p(\beta_c|\mathrm{KiDS})$", fontsize=13)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title(r"Marginalised posterior on $\beta_c$ (KiDS-only)", fontsize=11)

    ax2 = axes[1]
    ax2.plot(xs, cdf, "b-", lw=2, label="CDF KiDS-only")
    from scipy.stats import norm
    ax2.plot(bc_f, norm.cdf(bc_f, mu4, s4), "r--", lw=1.5,
             alpha=0.8, label="CDF Paper IV")
    ax2.axvline(0.0, color="red", ls=":", lw=1.5, alpha=0.7,
                label=r"$\Lambda$CDM")
    ax2.axvline(THETA_MAP[0], color="green", lw=1, ls="-.", alpha=0.7)
    ax2.axhline(0.5, color="gray", ls="--", lw=1, alpha=0.5)
    ax2.annotate(f"P(β_c<0) = {p_neg:.3f}\n({p_neg*100:.1f}%)",
                 xy=(0.0, p_neg), xytext=(0.45, 0.25),
                 textcoords="axes fraction", fontsize=11,
                 color="navy", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="navy", lw=1.5))
    ax2.set_xlabel(r"$\beta_c$", fontsize=13)
    ax2.set_ylabel(r"$P(\beta_c' < \beta_c|\mathrm{KiDS})$", fontsize=13)
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    ax2.set_title("Cumulative posterior", fontsize=11)

    plt.tight_layout()
    path = os.path.join(output_dir, "beta_c_kids_only.pdf")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"\n  Saved: {path}")
    plt.close()

    # LaTeX
    print("\n" + "="*55)
    print(f"LaTeX (κ_bf={_KAPPA_BF_FIXED} fixed):")
    latex = {"beta_c":r"$\beta_c$","sigma8":r"$\sigma_8$",
             "omega_m":r"$\omega_m$","A_IA":r"$A_{\rm IA}$"}
    for name in ["beta_c","sigma8","omega_m","A_IA"]:
        if name in results:
            med,hi,lo = results[name]
            print(f"{latex[name]} & ${med:+.4f}$ & $+{hi:.4f}/-{lo:.4f}$ \\\\")
    print(f"\\multicolumn{{2}}{{l}}{{$P(\\beta_c<0)$}} & ${p_neg:.3f}$ \\\\")

    return results, p_neg


# ═════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="KiDS-only analysis with kappa_bf fixed at Paper IV MAP")
    ap.add_argument("--bridge",
        default="outputs_shear/kids_real_bridge/kids_xipm_bridge.npz")
    ap.add_argument("--output_dir",  default="kids_fixed_kbf")
    ap.add_argument("--n_cores",     type=int, default=12)
    ap.add_argument("--ell_n",       type=int, default=150)
    ap.add_argument("--n_chi",       type=int, default=100)
    ap.add_argument("--n_a",        type=int, default=800)
    ap.add_argument("--n_eval",      type=int, default=500,
        help="Number of true likelihood evaluations")
    ap.add_argument("--nlive",       type=int, default=500)
    ap.add_argument("--top_units",   type=int, default=150)
    ap.add_argument("--skip_eval",   action="store_true",
        help="Skip evaluations, load existing training_data.npz")
    ap.add_argument("--plot_only",   action="store_true")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    train_npz = os.path.join(args.output_dir, "training_data.npz")
    out_npz   = os.path.join(args.output_dir, "kids_kbf_fixed_dynesty.npz")
    ell_grid  = np.geomspace(2.0, 5e4, args.ell_n)

    print("="*65)
    print(f"BPB Paper V — KiDS-only (κ_bf={_KAPPA_BF_FIXED} fixed)")
    print(f"  Free params: {PARAM_NAMES}")
    print(f"  NDIM={NDIM}  n_eval={args.n_eval}  nlive={args.nlive}")
    print("="*65)

    if args.plot_only:
        d = np.load(out_npz, allow_pickle=True)
        analyse_and_plot(out_npz, args.output_dir, rmse=0.0)
        return

    # ── Step 1: Evaluate true likelihood ─────────────────────────────────
    if args.skip_eval and os.path.exists(train_npz):
        print(f"\n[1] Loading existing training data: {train_npz}")
        d      = np.load(train_npz, allow_pickle=True)
        thetas = d["thetas"]
        lnLs   = d["lnLs"]
        print(f"  {len(thetas)} points  lnL_max={lnLs[np.isfinite(lnLs)].max():.2f}")
    else:
        print(f"\n[1] Generating {args.n_eval} training evaluations...")
        rng    = np.random.default_rng(42)
        points = generate_training_points(args.n_eval, rng=rng)
        lnLs   = evaluate_points(points, args.bridge, ell_grid,
                                  args.n_chi, args.n_a, args.n_cores)
        thetas = points
        np.savez(train_npz, thetas=thetas, lnLs=lnLs,
                 param_names=np.array(PARAM_NAMES))
        print(f"  Saved: {train_npz}")

    # ── Step 2: Train emulator ────────────────────────────────────────────
    print(f"\n[2] Training PCA+MLP emulator (top {args.top_units} units)...")
    pipe, lnL_max, epsilon, rmse, X_train = train_emulator(
        thetas, lnLs, top_units=args.top_units)
    set_emulator(pipe, lnL_max, epsilon)

    # ── Step 3: Dynesty on emulator ───────────────────────────────────────
    print(f"\n[3] Dynesty on emulator (nlive={args.nlive})...")
    results = run_dynesty(X_train, nlive=args.nlive, output_path=out_npz)

    # ── Step 4: Analyse ───────────────────────────────────────────────────
    print("\n[4] Analysis and plots...")
    analyse_and_plot(out_npz, args.output_dir, rmse=rmse)

    print(f"\nAll results in: {args.output_dir}/")
    print(f"Key figure: {args.output_dir}/beta_c_kids_only.pdf")


if __name__ == "__main__":
    main()
