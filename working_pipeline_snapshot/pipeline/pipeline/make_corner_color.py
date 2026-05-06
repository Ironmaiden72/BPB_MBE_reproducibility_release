"""
BPB/MBE — Corner plot en couleur pour Paper IV
Usage: python3 make_corner_color.py bpb_allphases_chain.h5 [output.pdf]
"""
import sys, numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    import corner # type: ignore
    import emcee # type: ignore
except ImportError:
    print("pip install corner emcee h5py"); sys.exit(1)

# ── Style ─────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.family':     'serif',
    'font.size':       11,
    'axes.labelsize':  13,
    'axes.titlesize':  11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'text.usetex':     False,
})

PARAM_NAMES = ['omega_m','omega_b','H0','a_bg','Delta_bg',
               'kappa_bf','beta_c','sigma8_0']

LABELS = [
    r'$\omega_m$',
    r'$\omega_b$',
    r'$H_0$',
    r'$a_{\rm bg}$',
    r'$\Delta_{\rm bg}$',
    r'$\kappa_{\rm bf}$',
    r'$\beta_c$',
    r'$\sigma_8^{(0)}$',
]

# ── Couleurs BPB ──────────────────────────────────────────────
COLOR_1S  = '#1f6fba'   # bleu foncé — 1σ
COLOR_2S  = '#5ba3d9'   # bleu moyen — 2σ
COLOR_3S  = '#aed4f0'   # bleu clair — 3σ
COLOR_HIST = '#1f6fba'
COLOR_TRUTH = '#e63946'  # rouge — valeur zéro (LCDM)

def make_corner_color(chainfile, outfile=None, discard=None, thin=5):
    reader = emcee.backends.HDFBackend(chainfile, read_only=True)
    n = reader.iteration
    print(f"Chain: {n} steps, {reader.shape[0]} walkers")

    if discard is None:
        discard = min(500, n // 10)

    chain = reader.get_chain(flat=True, discard=discard, thin=thin)
    lp    = reader.get_log_prob(flat=True, discard=discard, thin=thin)
    print(f"Samples after burn-in/thin: {chain.shape[0]}")

    # Posteriors
    print("\nPosteriors:")
    for i, name in enumerate(PARAM_NAMES):
        v = chain[:,i]
        med = np.median(v)
        lo  = med - np.percentile(v, 16)
        hi  = np.percentile(v, 84) - med
        print(f"  {name:12s} = {med:.4f} +{hi:.4f} -{lo:.4f}")

    # ── Corner plot ───────────────────────────────────────────
    # Truths: LCDM values (beta_c=0, kappa_bf=0)
    truths = [None, None, None, None, None, 0.0, 0.0, None]

    fig = corner.corner(
        chain,
        labels=LABELS,
        quantiles=[0.16, 0.50, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 10},
        title_fmt='.3f',
        truths=truths,
        truth_color=COLOR_TRUTH,
        color=COLOR_1S,
        hist_kwargs={
            'color': COLOR_HIST,
            'linewidth': 1.5,
        },
        contour_kwargs={
            'linewidths': [1.5, 1.2, 0.8],
        },
        levels=(0.393, 0.865, 0.989),  # 1σ, 2σ, 3σ
        fill_contours=True,
        contourf_kwargs={
            'colors': [COLOR_3S, COLOR_2S, COLOR_1S],
            'alpha': 0.8,
        },
        plot_datapoints=False,
        smooth=1.0,
        smooth1d=1.0,
    )

    # ── Titre et annotation ───────────────────────────────────
    fig.suptitle(
        r'BPB/MBE All-phases posterior (Phases 1+2+3+4+5)',
        fontsize=13, y=1.01, fontweight='bold'
    )

    # Annotation beta_c
    bc = chain[:,6]
    bc_med = np.median(bc)
    bc_lo  = bc_med - np.percentile(bc, 16)
    bc_hi  = np.percentile(bc, 84) - bc_med
    sigma  = bc_med / np.std(bc)

    fig.text(0.72, 0.85,
        f'$\\beta_c = {bc_med:.3f}^{{+{bc_hi:.3f}}}_{{-{bc_lo:.3f}}}$\n'
        f'$\\beta_c < 0$ at ${abs(bc_med/np.mean([bc_lo,bc_hi])):.1f}\\sigma$\n'
        f'5 independent probes',
        fontsize=12,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f4ff',
                  edgecolor='#1f6fba', linewidth=1.5),
        ha='center', va='center'
    )

    fig.text(0.72, 0.70,
        r'$\sigma_8^{(0)} = ' + f'{np.median(chain[:,7]):.3f}'
        r'\pm' + f'{np.std(chain[:,7]):.3f}$',
        fontsize=11, ha='center', color='#333333'
    )

    # Red dashed = LCDM (beta_c=0, kappa_bf=0)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0],[0], color=COLOR_TRUTH, ls='--', lw=1.5, label=r'$\Lambda$CDM ($\beta_c=0$)'),
        mpl.patches.Patch(facecolor=COLOR_1S, label=r'BPB/MBE $1\sigma$'), # type: ignore
        mpl.patches.Patch(facecolor=COLOR_2S, label=r'BPB/MBE $2\sigma$'), # type: ignore
        mpl.patches.Patch(facecolor=COLOR_3S, label=r'BPB/MBE $3\sigma$'), # type: ignore
    ]
    fig.legend(handles=legend_elements, loc='upper right',
               bbox_to_anchor=(0.99, 0.99), fontsize=10,
               framealpha=0.95, edgecolor='#cccccc')

    # ── Save ─────────────────────────────────────────────────
    if outfile is None:
        outfile = chainfile.replace('.h5', '_corner_color.pdf')

    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {outfile}")
    plt.close(fig)
    return outfile

if __name__ == "__main__":
    chainfile = sys.argv[1] if len(sys.argv) > 1 else "bpb_allphases_chain.h5"
    outfile   = sys.argv[2] if len(sys.argv) > 2 else None
    make_corner_color(chainfile, outfile)
