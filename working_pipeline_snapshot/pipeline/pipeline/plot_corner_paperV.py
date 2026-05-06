"""
plot_corner_paperV.py
=====================
Corner plot 3x3 (beta_c, sigma8, omega_m) pour Paper V.

Usage:
    python plot_corner_paperV.py \
        --npz kids_fixed_kbf_10k/kids_kbf_fixed_dynesty.npz \
        --output_dir kids_fixed_kbf_10k
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

# ── Paramètres ────────────────────────────────────────────────────────────
PARAMS   = ["beta_c", "sigma8", "omega_m"]
LABELS   = [r"$\beta_c$", r"$\sigma_8$", r"$\omega_m \equiv \Omega_m h^2$"]
# H0 fixed at Paper IV MAP — needed to convert omega_m → Omega_m for S8
_H0_FIXED = 70.35
_h2       = (_H0_FIXED / 100.0)**2

PAPER_IV = {
    "beta_c" : (-0.146, 0.033, 0.039),
    "sigma8" : ( 0.819, 0.011, 0.011),
    "omega_m": ( 0.1717, 0.012, 0.012),   # ω_m = Ω_m h² ≈ 0.172 at Paper IV MAP
}

BLUE  = "#1f6bb0"
LBLUE = "#a8c8e8"
RED   = "#cc3333"


def weighted_quantile(x, q, w):
    idx = np.argsort(x)
    xs, ws = x[idx], w[idx]
    cdf = np.cumsum(ws); cdf /= cdf[-1]
    return np.interp(q, cdf, xs)


def kde2d(x, y, w, ngrid=80):
    """Weighted 2D KDE on a grid."""
    xg = np.linspace(x.min(), x.max(), ngrid)
    yg = np.linspace(y.min(), y.max(), ngrid)
    xx, yy = np.meshgrid(xg, yg)
    pts = np.vstack([x, y])
    kde = gaussian_kde(pts, weights=w, bw_method="scott")
    zz  = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(ngrid, ngrid)
    return xg, yg, zz


def contour_levels(z, levels=(0.6827, 0.9545)):
    """Compute z-thresholds enclosing given probability mass."""
    z_flat = np.sort(z.ravel())[::-1]
    cumsum  = np.cumsum(z_flat); cumsum /= cumsum[-1]
    thresholds = []
    for lv in levels:
        idx = np.searchsorted(cumsum, lv)
        thresholds.append(z_flat[min(idx, len(z_flat)-1)])
    return thresholds[::-1]


def plot_corner(npz_path, output_dir):
    d       = np.load(npz_path, allow_pickle=True)
    samples = d["samples"]          # (N, 14)
    weights = d["weights"]
    pnames  = list(d["param_names"])

    w = weights / weights.sum()

    # Extraire les 3 colonnes + S8 dérivé
    idx  = [pnames.index(p) for p in PARAMS]
    base = samples[:, idx]          # (N, 3)

    # S8_eff : dérivé avec Omega_m = omega_m / h²  (H0 fixé)
    # NOTE: omega_m est le paramètre physique ω_m = Ω_m h², PAS Ω_m
    # S8_eff n'est PAS directement comparable au S8 KiDS standard sans cette correction
    Om_col  = base[:, 2] / _h2              # Ω_m vrai
    S8_col  = base[:, 1] * np.sqrt(Om_col / 0.3)   # S8_eff correct
    data   = np.column_stack([base, S8_col])   # (N, 4)
    ALL_PARAMS = PARAMS + ["S8_eff"]
    ALL_LABELS = LABELS + [r"$S_8^{\rm eff} \equiv \sigma_8\sqrt{\Omega_m/0.3}$" + "\n" + r"($\Omega_m = \omega_m/h^2$, $H_0$ fixed)"]
    N    = len(ALL_PARAMS)

    fig, axes = plt.subplots(N, N, figsize=(10, 9))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    fig.suptitle(
        r"BPB/MBE Paper V — Corner ($\beta_c$, $\sigma_8$, $\omega_m$, $S_8^{\rm eff}$)"
        "\nKiDS-1000 only, $\\kappa_{\\rm bf}=2.706$ fixed",
        fontsize=11, y=1.01)

    for i in range(N):       # row
        for j in range(N):   # col
            ax = axes[i, j]

            if j > i:        # upper triangle → hide
                ax.set_visible(False)
                continue

            xi = data[:, j]
            yi = data[:, i]

            if i == j:
                # ── Diagonal : 1D posterior ──────────────────────────────
                kde1 = gaussian_kde(xi, weights=w, bw_method="scott")
                xg   = np.linspace(xi.min(), xi.max(), 400)
                _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
                pdf  = kde1(xg); pdf /= _trapz(pdf, xg)

                ax.plot(xg, pdf, color=BLUE, lw=2)
                ax.fill_between(xg, 0, pdf, alpha=0.25, color=BLUE)

                # 1σ shading
                q16 = weighted_quantile(xi, 0.16, w)
                q84 = weighted_quantile(xi, 0.84, w)
                mask = (xg >= q16) & (xg <= q84)
                ax.fill_between(xg[mask], 0, pdf[mask], alpha=0.45, color=BLUE)

                # Median
                med = weighted_quantile(xi, 0.50, w)
                ax.axvline(med, color=BLUE, lw=1.5, ls="--")

                # Paper IV value
                pname = ALL_PARAMS[i]
                if pname in PAPER_IV:
                    piv_med = PAPER_IV[pname][0]
                    ax.axvline(piv_med, color=RED, lw=1.2, ls=":", alpha=0.8)

                # beta_c = 0 line
                if pname == "beta_c":
                    ax.axvline(0, color="gray", lw=1, ls=":", alpha=0.6)
                # S8_eff line (KiDS S8=0.759 mapped to eff space)
                if pname == "S8_eff":
                    ax.axvline(1.52, color="orange", lw=1.2, ls="--", alpha=0.8)

                ax.set_yticks([])
                ymax = pdf.max()*1.15
                ax.set_ylim(0, ymax)

            else:
                # ── Lower triangle : 2D contours ─────────────────────────
                xg, yg, zz = kde2d(xi, yi, w, ngrid=80)
                lvls = contour_levels(zz, levels=(0.6827, 0.9545))

                cmap = LinearSegmentedColormap.from_list(
                    "bpb", ["white", LBLUE, BLUE])

                ax.contourf(xg, yg, zz, levels=[lvls[0], zz.max()],
                            colors=[LBLUE], alpha=0.45)
                ax.contourf(xg, yg, zz, levels=[lvls[1], zz.max()],
                            colors=[BLUE],  alpha=0.30)
                ax.contour(xg, yg, zz,  levels=lvls,
                           colors=[BLUE, BLUE], linewidths=[1.0, 1.5])

                # Paper IV cross-hair
                pi_x = PAPER_IV.get(ALL_PARAMS[j])
                pi_y = PAPER_IV.get(ALL_PARAMS[i])
                if pi_x and pi_y:
                    ax.axvline(pi_x[0], color=RED, lw=0.8, ls=":", alpha=0.7)
                    ax.axhline(pi_y[0], color=RED, lw=0.8, ls=":", alpha=0.7)
                    ax.plot(pi_x[0], pi_y[0], "x",
                            color=RED, ms=6, mew=1.5, alpha=0.8)

                # beta_c = 0 line
                if ALL_PARAMS[j] == "beta_c":
                    ax.axvline(0, color="gray", lw=0.8, ls=":", alpha=0.5)
                if ALL_PARAMS[i] == "beta_c":
                    ax.axhline(0, color="gray", lw=0.8, ls=":", alpha=0.5)
                # S8_eff reference line
                if ALL_PARAMS[j] == "S8_eff":
                    ax.axvline(1.52, color="orange", lw=0.8, ls="--", alpha=0.6)
                if ALL_PARAMS[i] == "S8_eff":
                    ax.axhline(1.52, color="orange", lw=0.8, ls="--", alpha=0.6)

            # ── Axes labels ───────────────────────────────────────────────
            if i == N-1:
                ax.set_xlabel(ALL_LABELS[j], fontsize=11)
            else:
                ax.set_xticklabels([])

            if j == 0 and i != 0:
                ax.set_ylabel(ALL_LABELS[i], fontsize=11)
            else:
                ax.set_yticklabels([])

            ax.tick_params(labelsize=8)

    # ── Légende globale ───────────────────────────────────────────────────
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=BLUE,  alpha=0.5, label="KiDS-only 1σ"),
        Patch(facecolor=LBLUE, alpha=0.5, label="KiDS-only 2σ"),
        Line2D([0],[0], color=RED,    ls=":", lw=1.5, label="Paper IV MAP"),
        Line2D([0],[0], color="gray", ls=":", lw=1,   label=r"$\beta_c=0$ (ΛCDM)"),
        Line2D([0],[0], color="orange", ls="--", lw=1.2,
               label=r"KiDS $S_8$ reference mapped to internal $S_8^{\rm eff}$ (see text)"),
    ]
    fig.legend(handles=legend_elements, loc="upper right",
               bbox_to_anchor=(0.98, 0.98), fontsize=9, framealpha=0.9)

    # ── Stats en titre ────────────────────────────────────────────────────
    bc    = data[:, 0]
    med   = weighted_quantile(bc, 0.50, w)
    q16   = weighted_quantile(bc, 0.16, w)
    q84   = weighted_quantile(bc, 0.84, w)
    p_neg = float(np.interp(0.0,
                            np.sort(bc),
                            np.cumsum(w[np.argsort(bc)])))

    s8_col  = data[:, 3]
    s8_med  = weighted_quantile(s8_col, 0.50, w)
    s8_q16  = weighted_quantile(s8_col, 0.16, w)
    s8_q84  = weighted_quantile(s8_col, 0.84, w)

    stats_txt = (
        rf"$\beta_c = {med:+.3f}^{{+{q84-med:.3f}}}_{{-{med-q16:.3f}}}$"
        f"\n$P(\\beta_c<0) = {p_neg*100:.1f}\\%$"
        f"\n$S_8^{{\\rm eff}} = {s8_med:.3f}^{{+{s8_q84-s8_med:.3f}}}_{{-{s8_med-s8_q16:.3f}}}$"
        "\n(see text for $S_8$ caveat)"
    )
    fig.text(0.72, 0.72, stats_txt, fontsize=10,
             ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.4", fc="white",
                       ec=BLUE, alpha=0.9))

    import os
    out = os.path.join(output_dir, "corner_beta_s8_om_S8.pdf")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    plt.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz",        required=True)
    ap.add_argument("--output_dir", default=".")
    args = ap.parse_args()
    plot_corner(args.npz, args.output_dir)
