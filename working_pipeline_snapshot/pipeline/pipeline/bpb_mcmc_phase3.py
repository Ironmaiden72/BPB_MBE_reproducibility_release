"""
BPB/MBE Phase 3 standalone MCMC — ACT DR6 + Planck CMB lensing only
====================================================================
Runs MCMC with ONLY the CMB lensing A_lens constraint.
Shows which BPB parameter region is consistent with ACT alone.

Usage:
    python3 bpb_mcmc_phase3.py --bestfit   # MAP only (~2 min)
    python3 bpb_mcmc_phase3.py             # full run (32 walkers × 3000 steps)
    python3 bpb_mcmc_phase3.py --plot chain_p3.h5

Parameters: same 8-vector as Phase 1+2
    [omega_m, omega_b, H0, a_bg, Delta_bg, kappa_bf, beta_c, sigma8_0]
"""

import sys, time
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, '.')
from bpb_phase3_act import lnL_phase3, A_lens_theory, CMB_LENS_MU, CMB_LENS_SIG

try:
    import emcee
    HAS_EMCEE = True
except ImportError:
    HAS_EMCEE = False
    print("WARNING: emcee not found. pip install emcee h5py")

try:
    import corner
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

# ── Parameters & priors ──────────────────────────────────────
PARAM_NAMES = ['omega_m','omega_b','H0','a_bg','Delta_bg',
               'kappa_bf','beta_c','sigma8_0']
LABELS = [r'$\omega_m$', r'$\omega_b$', r'$H_0$',
          r'$a_{\rm bg}$', r'$\Delta_{\rm bg}$',
          r'$\kappa_{\rm bf}$', r'$\beta_c$', r'$\sigma_8^{(0)}$']

PRIOR_BOUNDS = np.array([
    [0.10,  0.20 ],
    [0.019, 0.025],
    [60.0,  80.0 ],
    [0.10,  0.90 ],
    [0.05,  2.0  ],
    [0.0,   5.0  ],
    [-0.5,  0.5  ],
    [0.65,  1.05 ],
])
OB_MU, OB_SIG = 0.02237, 0.00015

# Starting point — near Planck to help convergence
THETA0 = np.array([0.143, 0.02237, 67.4, 0.5, 0.4, 0.3, 0.0, 0.81])

def ln_prior(theta):
    for i,(lo,hi) in enumerate(PRIOR_BOUNDS):
        if not (lo <= theta[i] <= hi): return -np.inf
    if 1.0 + theta[6] <= 0: return -np.inf
    return -0.5*((theta[1]-OB_MU)/OB_SIG)**2

def ln_posterior(theta):
    lp = ln_prior(theta)
    if not np.isfinite(lp): return -np.inf
    try:
        ll = lnL_phase3(theta)
        return lp + ll if np.isfinite(ll) else -np.inf
    except Exception:
        return -np.inf

# ── MAP ──────────────────────────────────────────────────────
def find_map(theta0=THETA0, verbose=True):
    if verbose: print("Finding MAP (ACT only)...")
    t0 = time.time()
    res = minimize(lambda t: -ln_posterior(t), theta0,
                   method='Nelder-Mead',
                   options={'maxiter':5000,'xatol':1e-4,'fatol':1e-4,'adaptive':True})
    if verbose:
        print(f"  MAP lnP = {-res.fun:.3f}  ({time.time()-t0:.1f}s)")
        for n,v in zip(PARAM_NAMES, res.x):
            print(f"  {n:12s} = {v:.5f}")
        A = A_lens_theory(res.x)
        print(f"\n  A_lens at MAP = {A:.4f}  (data = {CMB_LENS_MU:.4f} ± {CMB_LENS_SIG:.4f})")
        print(f"  Pull = {(A-CMB_LENS_MU)/CMB_LENS_SIG:+.2f} sigma")
    return res.x, -res.fun

# ── emcee ────────────────────────────────────────────────────
def run_emcee(theta_map, n_walkers=32, n_steps=3000, n_burn=300,
              backend_file="chain_p3.h5", verbose=True):
    if not HAS_EMCEE: raise ImportError("pip install emcee h5py")
    ndim   = len(theta_map)
    scales = 0.005*(PRIOR_BOUNDS[:,1]-PRIOR_BOUNDS[:,0])
    p0     = theta_map + scales*np.random.randn(n_walkers, ndim)
    for iw in range(n_walkers):
        for ip in range(ndim):
            lo,hi = PRIOR_BOUNDS[ip]
            p0[iw,ip] = np.clip(p0[iw,ip], lo+1e-6, hi-1e-6)
    backend = emcee.backends.HDFBackend(backend_file)
    backend.reset(n_walkers, ndim)
    sampler = emcee.EnsembleSampler(n_walkers, ndim, ln_posterior, backend=backend)
    if verbose:
        print(f"\nRunning emcee (ACT only): {n_walkers}w × {n_steps} steps")
        t0 = time.time()
    state = sampler.run_mcmc(p0, n_burn, progress=verbose, store=False)
    sampler.reset()
    sampler.run_mcmc(state, n_steps, progress=verbose)
    if verbose:
        print(f"\nDone in {(time.time()-t0)/60:.1f} min")
        tau = sampler.get_autocorr_time(quiet=True)
        print(f"Autocorr: {np.round(tau,1)}")
        print(f"N/tau_max = {sampler.iteration/max(tau):.1f}")
        print(f"Acceptance: {np.mean(sampler.acceptance_fraction):.3f}")
    return sampler

# ── Summary ──────────────────────────────────────────────────
def print_summary(chain, lp=None):
    print("\n"+"="*60)
    print("PHASE 3 (ACT only) MCMC SUMMARY")
    print("="*60)
    print(f"{'Parameter':12s}  {'Mean':>10s}  {'Std':>9s}  {'16%':>10s}  {'84%':>10s}")
    print("-"*60)
    for i,name in enumerate(PARAM_NAMES):
        v = chain[:,i]
        print(f"{name:12s}  {np.mean(v):10.4f}  {np.std(v):9.4f}"
              f"  {np.percentile(v,16):10.4f}  {np.percentile(v,84):10.4f}")
    if lp is not None:
        idx = np.argmax(lp)
        print(f"\nBest-fit A_lens = {A_lens_theory(chain[idx]):.4f}")
        print(f"  (data: {CMB_LENS_MU:.4f} ± {CMB_LENS_SIG:.4f})")

def make_plots(chain, lp, out_prefix="p3"):
    if not HAS_PLOT:
        print("Install matplotlib+corner for plots")
        return
    # Corner
    fig = corner.corner(chain, labels=LABELS,
                        quantiles=[0.16,0.5,0.84], show_titles=True)
    fig.savefig(f"{out_prefix}_corner.pdf", dpi=150, bbox_inches='tight')
    print(f"Corner: {out_prefix}_corner.pdf")
    plt.close()
    # A_lens marginal
    A_arr = np.array([A_lens_theory(chain[i]) for i in
                      np.random.choice(len(chain), min(500,len(chain)), replace=False)])
    fig,ax = plt.subplots(figsize=(6,4))
    ax.hist(A_arr, bins=40, density=True, alpha=0.7, color='steelblue', label='BPB posterior')
    x = np.linspace(0.8, 1.4, 200)
    ax.plot(x, np.exp(-0.5*((x-CMB_LENS_MU)/CMB_LENS_SIG)**2)/
            (CMB_LENS_SIG*np.sqrt(2*np.pi)), 'r-', lw=2, label='ACT+Planck data')
    ax.axvline(1.0, color='k', ls='--', lw=1, label='ΛCDM (=1)')
    ax.set_xlabel(r'$A_{\rm lens}$'); ax.set_ylabel('Posterior density')
    ax.legend(fontsize=9); ax.set_title('CMB lensing amplitude — Phase 3')
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_Alens.pdf", dpi=150, bbox_inches='tight')
    print(f"A_lens: {out_prefix}_Alens.pdf")
    plt.close()

# ── Main ─────────────────────────────────────────────────────
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "--full"
    np.random.seed(42)

    if mode == "--bestfit":
        theta_map, _ = find_map()

    elif mode == "--plot":
        chainfile = sys.argv[2] if len(sys.argv) > 2 else "chain_p3.h5"
        reader = emcee.backends.HDFBackend(chainfile, read_only=True)
        n = reader.iteration
        discard = min(300, n//10)
        chain = reader.get_chain(flat=True, discard=discard, thin=3)
        lp    = reader.get_log_prob(flat=True, discard=discard, thin=3)
        print(f"Chain: {n} steps → {chain.shape[0]} samples after thinning")
        print_summary(chain, lp)
        make_plots(chain, lp)

    else:  # --full
        if not HAS_EMCEE:
            print("pip install emcee h5py"); sys.exit(1)
        theta_map, _ = find_map()
        sampler = run_emcee(theta_map, n_walkers=32, n_steps=3000,
                            n_burn=300, backend_file="chain_p3.h5")
        chain = sampler.get_chain(flat=True, discard=100, thin=3)
        lp    = sampler.get_log_prob(flat=True, discard=100, thin=3)
        print_summary(chain, lp)
        make_plots(chain, lp)
        np.savez("p3_results.npz", chain=chain, lnpost=lp,
                 param_names=np.array(PARAM_NAMES))
        print("Saved: p3_results.npz")
