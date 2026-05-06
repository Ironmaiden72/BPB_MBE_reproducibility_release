#!/usr/bin/env python3
"""
Paper VI / E_G T7 — Wenzl model-space collinearity obstruction diagnostic.

Purpose
-------
T6c built a model-space inner product over the Wenzl projected 19-bin response
curves. It found near-perfect collinearity among W_mu, W_ACTSigma, and W_EGSigma.

T7 formalizes that as an obstruction:
- The Wenzl 19-bin model-response space is effectively one-dimensional in shape.
- A naive projector using e_m = normalize(W_mu) gives overlap fractions ~1.
- Therefore this space cannot by itself expose the two-mode fraction a_bg.
- The r_Wdelta≈sqrt(a_bg) information must live in the auto/cross amplitude
  diagnostic space or a richer Weyl/proxy space, not in the Wenzl shape curves alone.

Non-destructive:
- Reads existing cosmology outputs.
- Writes only outputs/paperVI_EG_theory/.
- Touches no NS files.
- Runs no likelihood.
"""
from pathlib import Path
import csv
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIMARY = Path("outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.csv")
AUTO_CROSS = Path("outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv")
FREE_R_SCAN = Path("outputs/cosmology/autocross_csigma_free_r_scan_final.csv")
T1_JSON = OUTDIR / "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json"

VECTOR_COLS = {
    "W_LCDM": "EG_lcdm_projected",
    "W_mu": "EG_bpb_mu_projected",
    "W_ACTSigma": "EG_bpb_act_sigma_projected",
    "W_EGSigma": "EG_bpb_eg_sigma_projected",
}


def read_csv_dict(path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        return list(reader), reader.fieldnames or []


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except Exception:
        return default


def dot(a, b, w=None):
    if w is None:
        return sum(x*y for x, y in zip(a, b)) / len(a)
    return sum(wi*x*y for wi, x, y in zip(w, a, b))


def norm(a, w=None):
    v = dot(a, a, w)
    return math.sqrt(v) if v >= 0 else float("nan")


def cosine(a, b, w=None):
    na = norm(a, w)
    nb = norm(b, w)
    if na == 0 or nb == 0:
        return float("nan")
    return dot(a, b, w) / (na * nb)


def project_fraction(v, e, w=None):
    # fraction of ||v||^2 captured by projection onto e.
    # e need not be normalized.
    ee = dot(e, e, w)
    vv = dot(v, v, w)
    ev = dot(e, v, w)
    if ee <= 0 or vv <= 0:
        return float("nan"), float("nan"), float("nan")
    coeff = ev / ee
    frac = (ev * ev) / (ee * vv)
    # residual norm fraction
    resid2 = max(0.0, vv - coeff * coeff * ee)
    resid_frac = resid2 / vv
    return frac, resid_frac, coeff


def mat_eig_4x4_power_approx(gram, n_iter=100):
    # No numpy dependency needed, but repository has numpy. Use pure python for portability.
    # Returns a simple dominant-eigenvalue estimate and trace fraction.
    n = len(gram)
    v = [1.0 / math.sqrt(n)] * n
    lam = 0.0
    for _ in range(n_iter):
        y = [sum(gram[i][j] * v[j] for j in range(n)) for i in range(n)]
        ny = math.sqrt(sum(x*x for x in y))
        if ny == 0:
            break
        v = [x/ny for x in y]
        lam = sum(v[i] * sum(gram[i][j]*v[j] for j in range(n)) for i in range(n))
    tr = sum(gram[i][i] for i in range(n))
    return lam, lam / tr if tr else float("nan")


def write_dict_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def read_t1_values():
    a_bg = 0.403687948
    sqrt_a = math.sqrt(a_bg)
    r_free = 0.630864746
    dchi2 = 0.008545
    if T1_JSON.exists():
        try:
            data = json.loads(T1_JSON.read_text(encoding="utf-8"))
            a_bg = float(data.get("inputs", {}).get("a_bg", a_bg))
            sqrt_a = float(data.get("derived", {}).get("sqrt_a_bg", math.sqrt(a_bg)))
            r_free = float(data.get("inputs", {}).get("r_wdelta_free", r_free))
            dchi2 = float(data.get("derived", {}).get("delta_chi2_fixed_sqrt_abg_minus_free", dchi2))
        except Exception:
            pass
    return a_bg, sqrt_a, r_free, dchi2


def main():
    a_bg, sqrt_a_bg, r_free_t1, dchi2_t1 = read_t1_values()

    if not PRIMARY.exists():
        status = "P6_T7_PRIMARY_WENZL_FILE_MISSING_FAIL"
        result = {
            "status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO"},
            "error": f"Missing primary file: {PRIMARY}",
        }
        (OUTDIR / "paperVI_EG_T7_collinearity_obstruction.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(status)
        return

    rows, cols = read_csv_dict(PRIMARY)
    ell_widths = [max(0.0, fnum(r.get("ell_hi")) - fnum(r.get("ell_lo"))) for r in rows]
    sw = sum(ell_widths)
    w_ell = [x / sw for x in ell_widths] if sw > 0 else None

    vectors = {}
    for name, col in VECTOR_COLS.items():
        vectors[name] = [fnum(r.get(col)) for r in rows]

    schemes = {
        "unweighted_mean": None,
        "ell_width_weighted_mean": w_ell,
    }

    corr_rows = []
    projection_rows = []
    rank_rows = []

    names = list(VECTOR_COLS.keys())

    for scheme, weights in schemes.items():
        gram = []
        for a in names:
            row = []
            for b in names:
                row.append(dot(vectors[a], vectors[b], weights))
                corr_rows.append({
                    "scheme": scheme,
                    "A": a,
                    "B": b,
                    "cosine": cosine(vectors[a], vectors[b], weights),
                })
            gram.append(row)

        lam1, lam1_trace_frac = mat_eig_4x4_power_approx(gram)
        rank_rows.append({
            "scheme": scheme,
            "dominant_eigenvalue_estimate": lam1,
            "trace": sum(gram[i][i] for i in range(len(gram))),
            "dominant_trace_fraction": lam1_trace_frac,
            "interpretation": "near_1D_if_close_to_1",
        })

        e_m = vectors["W_mu"]
        for target in ["W_LCDM", "W_ACTSigma", "W_EGSigma"]:
            frac, resid_frac, coeff = project_fraction(vectors[target], e_m, weights)
            projection_rows.append({
                "scheme": scheme,
                "matter_direction": "W_mu",
                "target": target,
                "parallel_fraction": frac,
                "orthogonal_fraction": resid_frac,
                "projection_coeff_on_W_mu": coeff,
                "delta_parallel_fraction_minus_a_bg": frac - a_bg,
                "obstruction": "YES" if abs(frac - a_bg) > 0.1 and frac > 0.99 else "CHECK",
            })

    # Supporting auto/cross diagnostic audit.
    autocross_rows = []
    r_values = []
    if AUTO_CROSS.exists():
        ac_rows, ac_cols = read_csv_dict(AUTO_CROSS)
        grouped = {}
        for r in ac_rows:
            ds = r.get("dataset", "unknown")
            grouped.setdefault(ds, []).append(r)
            rv = fnum(r.get("r_global"))
            if math.isfinite(rv):
                r_values.append(rv)
        for ds, grows in grouped.items():
            chi2 = 0.0
            n = 0
            rset = []
            for r in grows:
                pull = fnum(r.get("pull_global_auto_cross"))
                rv = fnum(r.get("r_global"))
                if math.isfinite(pull):
                    chi2 += pull * pull
                    n += 1
                if math.isfinite(rv):
                    rset.append(rv)
            rmean = sum(rset)/len(rset) if rset else float("nan")
            autocross_rows.append({
                "dataset": ds,
                "n": n,
                "r_global_mean": rmean,
                "r_global_squared": rmean*rmean if math.isfinite(rmean) else float("nan"),
                "delta_r2_minus_a_bg": (rmean*rmean - a_bg) if math.isfinite(rmean) else float("nan"),
                "chi2_diag_from_pulls": chi2,
                "chi2_per_point": chi2/n if n else float("nan"),
            })

    r_global_mean = sum(r_values)/len(r_values) if r_values else float("nan")
    r_global_sq = r_global_mean*r_global_mean if math.isfinite(r_global_mean) else float("nan")

    # Free-r scan supporting rows.
    free_scan_summary = {}
    if FREE_R_SCAN.exists():
        fs_rows, _ = read_csv_dict(FREE_R_SCAN)
        best = None
        for r in fs_rows:
            # Choose minimum chi2_EG_physical or delta_physical_minus_free if available.
            key = fnum(r.get("chi2_EG_physical"))
            if not math.isfinite(key):
                key = fnum(r.get("delta_physical_minus_free"))
            if math.isfinite(key):
                if best is None or key < best[0]:
                    best = (key, r)
        if best is not None:
            br = best[1]
            free_scan_summary = {
                "c_sigma_auto_best_by_chi2_EG_physical": fnum(br.get("c_sigma_auto")),
                "a_bg": fnum(br.get("a_bg")),
                "sqrt_a_bg": fnum(br.get("sqrt_a_bg")),
                "r_hat_free": fnum(br.get("r_hat_free")),
                "sigma_r": fnum(br.get("sigma_r")),
                "z_sqrt_abg_minus_rhat": fnum(br.get("z_sqrt_abg_minus_rhat")),
                "chi2_EG_free_r": fnum(br.get("chi2_EG_free_r")),
                "chi2_EG_physical": fnum(br.get("chi2_EG_physical")),
                "delta_physical_minus_free": fnum(br.get("delta_physical_minus_free")),
            }

    # Obstruction criteria.
    min_parallel = min(float(r["parallel_fraction"]) for r in projection_rows if r["target"] in {"W_ACTSigma", "W_EGSigma"})
    max_orth = max(float(r["orthogonal_fraction"]) for r in projection_rows if r["target"] in {"W_ACTSigma", "W_EGSigma"})
    one_d = all(float(r["dominant_trace_fraction"]) > 0.99999 for r in rank_rows)
    obstruction = (min_parallel > 0.999999 and max_orth < 1e-5 and one_d)

    if obstruction:
        status = "P6_T7_WENZL_MODEL_SPACE_COLLINEARITY_OBSTRUCTION_LOCKED"
    else:
        status = "P6_T7_WENZL_MODEL_SPACE_COLLINEARITY_CHECK"

    decisions = [
        {
            "item": "wenzl_model_space_projector_trial",
            "decision": "OBSTRUCTION_LOCKED" if obstruction else "CHECK",
            "basis": "W_mu, W_ACTSigma, and W_EGSigma are nearly perfectly collinear; naive projection gives parallel fraction near 1, not a_bg.",
        },
        {
            "item": "A4_test_in_wenzl_shape_space",
            "decision": "DO_NOT_USE_AS_A4_PROOF",
            "basis": "The 19-bin Wenzl projected response curves have essentially one shape degree of freedom; they cannot expose the matter/private two-mode split.",
        },
        {
            "item": "rWdelta_information_location",
            "decision": "AUTO_CROSS_OR_RAW_WEYL_SPACE",
            "basis": "The r_global diagnostic and free-r scan carry amplitude/cross information absent from the collinear Wenzl shape curves.",
        },
        {
            "item": "paper_use",
            "decision": "USE_AS_OBSTRUCTION_DIAGNOSTIC",
            "basis": "This gate prevents a false projector proof and redirects the derivation to the correct vector space.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_OR_RAW_WEYL_SPACE",
            "basis": "Build a vector space from auto/cross amplitude components or raw Weyl proxy, not from Wenzl shape curves alone.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "COLLINEARITY_OBSTRUCTION_DIAGNOSTIC_NOT_A4_PROOF",
        },
        "inputs": {
            "primary_wenzl_file": str(PRIMARY),
            "auto_cross_file": str(AUTO_CROSS),
            "free_r_scan_file": str(FREE_R_SCAN),
            "a_bg": a_bg,
            "sqrt_a_bg": sqrt_a_bg,
            "r_free_T1": r_free_t1,
            "delta_chi2_T1_fixed_sqrt_abg_minus_free": dchi2_t1,
        },
        "rank_diagnostics": rank_rows,
        "projection_diagnostics": projection_rows,
        "autocross_support": autocross_rows,
        "r_global_summary": {
            "r_global_mean": r_global_mean,
            "r_global_squared": r_global_sq,
            "delta_r_global_squared_minus_a_bg": r_global_sq - a_bg if math.isfinite(r_global_sq) else float("nan"),
        },
        "free_r_scan_support": free_scan_summary,
        "decisions": decisions,
    }

    # Write files.
    (OUTDIR / "paperVI_EG_T7_collinearity_obstruction.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    write_dict_csv(
        OUTDIR / "paperVI_EG_T7_cosine_matrix_rows.csv",
        corr_rows,
        ["scheme", "A", "B", "cosine"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T7_projection_fractions.csv",
        projection_rows,
        ["scheme", "matter_direction", "target", "parallel_fraction", "orthogonal_fraction", "projection_coeff_on_W_mu", "delta_parallel_fraction_minus_a_bg", "obstruction"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T7_rank_diagnostics.csv",
        rank_rows,
        ["scheme", "dominant_eigenvalue_estimate", "trace", "dominant_trace_fraction", "interpretation"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T7_autocross_support.csv",
        autocross_rows,
        ["dataset", "n", "r_global_mean", "r_global_squared", "delta_r2_minus_a_bg", "chi2_diag_from_pulls", "chi2_per_point"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T7_decisions.csv",
        decisions,
        ["item", "decision", "basis"],
    )

    # Markdown.
    md = f"""# Paper VI / E_G T7 - Wenzl model-space collinearity obstruction

## Status

**{status}**

Non-destructive diagnostic. No likelihood run. NS untouched.

## Purpose

T6c built a model-space inner product over the Wenzl projected response curves.  
T7 tests whether that space can carry the two-mode matter/private Weyl split required by A4:

    <P_m W,P_m W>_EG / <W,W>_EG = a_bg

It cannot, because the Wenzl 19-bin response curves are essentially collinear.

## Locked inputs

| quantity | value |
|---|---:|
| a_bg | {a_bg:.12f} |
| sqrt(a_bg) | {sqrt_a_bg:.12f} |
| T1 free r_Wdelta | {r_free_t1:.12f} |
| T1 delta chi2 sqrt(a_bg) vs free | {dchi2_t1:.6f} |

## Rank / one-dimensionality diagnostic

| scheme | dominant trace fraction |
|---|---:|
"""

    for r in rank_rows:
        md += f"| {r['scheme']} | {float(r['dominant_trace_fraction']):.12f} |\n"

    md += """
## Naive projection fractions using e_m = normalize(W_mu)

| scheme | target | parallel fraction | orthogonal fraction | delta parallel - a_bg |
|---|---|---:|---:|---:|
"""

    for r in projection_rows:
        if r["target"] in {"W_ACTSigma", "W_EGSigma"}:
            md += (
                f"| {r['scheme']} | {r['target']} | "
                f"{float(r['parallel_fraction']):.12f} | "
                f"{float(r['orthogonal_fraction']):.3e} | "
                f"{float(r['delta_parallel_fraction_minus_a_bg']):.12f} |\n"
            )

    md += """
## Interpretation

A naive Wenzl model-space projector gives parallel fractions essentially equal to 1.  
But the target A4 fraction is:

    a_bg ~= 0.404

Therefore the Wenzl projected model curves cannot by themselves expose the two-mode decomposition.  
They are useful model response curves, but not the right vector space for proving A4.

## Supporting auto/cross diagnostic

The auto/cross file carries the global-r information directly:

    outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv

| dataset | r_global | r_global^2 | r_global^2-a_bg | chi2/pt |
|---|---:|---:|---:|---:|
"""

    for r in autocross_rows:
        md += (
            f"| {r['dataset']} | {float(r['r_global_mean']):.12f} | "
            f"{float(r['r_global_squared']):.12f} | "
            f"{float(r['delta_r2_minus_a_bg']):.12f} | "
            f"{float(r['chi2_per_point']):.12f} |\n"
        )

    if free_scan_summary:
        md += f"""
## Free-r scan support

| quantity | value |
|---|---:|
| c_sigma_auto best | {free_scan_summary.get('c_sigma_auto_best_by_chi2_EG_physical', float('nan')):.12f} |
| sqrt(a_bg) | {free_scan_summary.get('sqrt_a_bg', float('nan')):.12f} |
| r_hat_free | {free_scan_summary.get('r_hat_free', float('nan')):.12f} |
| sigma_r | {free_scan_summary.get('sigma_r', float('nan')):.12f} |
| z_sqrt_abg_minus_rhat | {free_scan_summary.get('z_sqrt_abg_minus_rhat', float('nan')):.12f} |
| chi2_EG_free_r | {free_scan_summary.get('chi2_EG_free_r', float('nan')):.12f} |
| chi2_EG_physical | {free_scan_summary.get('chi2_EG_physical', float('nan')):.12f} |
| delta physical-free | {free_scan_summary.get('delta_physical_minus_free', float('nan')):.12e} |
"""
    md += """
## Decisions

| item | decision | basis |
|---|---|---|
"""
    for d in decisions:
        md += f"| {d['item']} | {d['decision']} | {d['basis']} |\n"

    md += """
## Locked conclusion

The Wenzl 19-bin projected model-response space is effectively one-dimensional in shape.  
It cannot be used to prove the A4 matter/private Weyl overlap fraction, because a naive projection onto W_mu returns overlap fractions near 1 rather than a_bg.

This is an obstruction diagnostic, not a failure of the r_Wdelta=sqrt(a_bg) closure.  
It shows where the relevant information lives: in the auto/cross amplitude diagnostic, the free-r scan, or a richer raw Weyl/proxy space.

## Next gate

P6-T8 should build the candidate vector space from either:

1. auto/cross amplitude components, especially `eg_projected_auto_cross_global_r_seed01.csv`; or
2. raw Weyl/lensing proxy vectors, especially `cmb_lensing_weyl_proxy_map1782_Planck_PR4_raw.csv`.

The goal is to find a space where the matter-correlated and private Weyl components are not collapsed into one collinear model-response shape.
"""

    (OUTDIR / "paperVI_EG_T7_collinearity_obstruction.md").write_text(md, encoding="utf-8")

    # Figure.
    schemes = sorted(set(r["scheme"] for r in projection_rows))
    targets = ["W_ACTSigma", "W_EGSigma"]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    x = []
    heights = []
    labels = []
    for s in schemes:
        for t in targets:
            row = next(r for r in projection_rows if r["scheme"] == s and r["target"] == t)
            x.append(len(x))
            heights.append(float(row["parallel_fraction"]))
            labels.append(f"{s}\n{t}")
    ax.bar(x, heights)
    ax.axhline(a_bg, linestyle="--", linewidth=1, label="a_bg target")
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("parallel fraction onto W_mu")
    ax.set_title("Paper VI T7: Wenzl shape-space projector obstruction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_P6_T7_collinearity_obstruction.png", dpi=180)
    fig.savefig(OUTDIR / "fig_P6_T7_collinearity_obstruction.pdf")
    plt.close(fig)

    print(status)
    print("min_parallel_sigma_on_W_mu:", min_parallel)
    print("max_orthogonal_fraction:", max_orth)
    print("r_global_mean:", r_global_mean)
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
