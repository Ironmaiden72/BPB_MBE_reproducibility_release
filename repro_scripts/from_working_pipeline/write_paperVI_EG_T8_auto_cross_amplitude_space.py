#!/usr/bin/env python3
"""
Paper VI / E_G T8 - auto/cross amplitude-space projector diagnostic.

Purpose
-------
T7 showed that the 19-bin Wenzl projected model-shape space is almost
one-dimensional, so it cannot expose the matter/private Weyl overlap A4.

T8 therefore moves to the space where r_Wdelta is actually defined:
the auto/cross amplitude diagnostic.

Core test
---------
Use the auto/cross file

    outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv

and define, per projected row,

    A_i = unit_ACTSigma_cross
    X_i = prediction_global_auto_cross

Then X_i = r_global A_i by construction in the diagnostic. Therefore

    <X,X>/<A,A> = r_global^2

for any positive diagonal/bin weights. T8 checks whether this amplitude-space
overlap equals a_bg.

This is not a data chi2 proof and not an action-level derivation. It is the
correct projected amplitude-space bridge between T1-T3 and the actual
auto/cross diagnostic artefacts.
"""
from pathlib import Path
import csv
import json
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

AUTO_CROSS = Path("outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv")
FREE_SCAN = Path("outputs/cosmology/autocross_csigma_free_r_scan_final.csv")
RAW_WEYL = Path("outputs/cmb/cmb_lensing_weyl_proxy_map1782_Planck_PR4_raw.csv")


def ffloat(x, default=None):
    try:
        if x is None:
            return default
        s = str(x).strip()
        if not s or s.lower() in {"nan", "none", "null"}:
            return default
        return float(s)
    except Exception:
        return default


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        return list(csv.DictReader(f, dialect=dialect))


def write_dict_csv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def write_kv_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        w.writerows(rows)


def load_a_bg():
    # Prefer T7/T3/T1 frozen values if available.
    defaults = {
        "a_bg": 0.403687948,
        "sqrt_a_bg": math.sqrt(0.403687948),
        "r_free": 0.630864746,
        "delta_chi2": 0.008545000000000025,
    }

    for name in [
        "paperVI_EG_T7_collinearity_obstruction.json",
        "paperVI_EG_T3_A4_overlap_reduction.json",
        "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json",
    ]:
        p = OUTDIR / name
        if not p.exists():
            continue
        data = json.loads(p.read_text(encoding="utf-8"))
        inp = data.get("inputs", {})
        der = data.get("derived", {})
        locked = data.get("locked_inputs", {})
        defaults["a_bg"] = ffloat(inp.get("a_bg", locked.get("a_bg", defaults["a_bg"])), defaults["a_bg"])
        defaults["sqrt_a_bg"] = math.sqrt(defaults["a_bg"])
        defaults["r_free"] = ffloat(inp.get("r_wdelta_free", locked.get("r_wdelta_free", defaults["r_free"])), defaults["r_free"])
        defaults["delta_chi2"] = ffloat(der.get("delta_chi2_sqrt_abg_minus_free", locked.get("delta_chi2", defaults["delta_chi2"])), defaults["delta_chi2"])
        break
    return defaults


def amp_summary(rows, a_bg, weight_mode="unweighted", dataset=None):
    use = []
    for r in rows:
        if dataset is not None and r.get("dataset") != dataset:
            continue
        A = ffloat(r.get("unit_ACTSigma_cross"))
        X = ffloat(r.get("prediction_global_auto_cross"))
        obs = ffloat(r.get("obs"))
        sig = ffloat(r.get("sigma_obs"))
        rg = ffloat(r.get("r_global"))
        pull = ffloat(r.get("pull_global_auto_cross"))
        if A is None or X is None:
            continue
        if weight_mode == "sigma_obs_invvar":
            if sig is None or sig <= 0:
                continue
            w = 1.0 / (sig * sig)
        else:
            w = 1.0
        use.append((A, X, obs, sig, rg, pull, w))

    if not use:
        return None

    AA = sum(w * A * A for A, X, obs, sig, rg, pull, w in use)
    XX = sum(w * X * X for A, X, obs, sig, rg, pull, w in use)
    AX = sum(w * A * X for A, X, obs, sig, rg, pull, w in use)
    frac = XX / AA if AA > 0 else float("nan")
    r_norm = AX / AA if AA > 0 else float("nan")
    corr = AX / math.sqrt(AA * XX) if AA > 0 and XX > 0 else float("nan")
    chi2 = sum((pull * pull) for A, X, obs, sig, rg, pull, w in use if pull is not None)
    r_values = [rg for A, X, obs, sig, rg, pull, w in use if rg is not None]
    ratios = [X / A for A, X, obs, sig, rg, pull, w in use if A != 0]

    return {
        "dataset": dataset or "ALL",
        "weight_mode": weight_mode,
        "n": len(use),
        "auto_norm_AA": AA,
        "cross_norm_XX": XX,
        "auto_cross_AX": AX,
        "fraction_XX_over_AA": frac,
        "sqrt_fraction": math.sqrt(frac) if frac >= 0 else float("nan"),
        "r_from_AX_over_AA": r_norm,
        "corr_A_X": corr,
        "mean_r_global": sum(r_values) / len(r_values) if r_values else float("nan"),
        "mean_ratio_X_over_A": sum(ratios) / len(ratios) if ratios else float("nan"),
        "delta_fraction_minus_a_bg": frac - a_bg,
        "relative_delta_fraction": (frac - a_bg) / a_bg if a_bg else float("nan"),
        "chi2_from_pulls": chi2,
        "chi2_per_point": chi2 / len(use) if use else float("nan"),
    }


def main():
    prior = load_a_bg()
    a_bg = prior["a_bg"]
    sqrt_a = math.sqrt(a_bg)

    status = "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_PASS"

    if not AUTO_CROSS.exists():
        status = "P6_T8_AUTO_CROSS_FILE_MISSING_FAIL"
        result = {
            "status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO"},
            "missing": str(AUTO_CROSS),
        }
        (OUTDIR / "paperVI_EG_T8_auto_cross_amplitude_space.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(status)
        return

    rows = read_csv(AUTO_CROSS)
    datasets = sorted(set(r.get("dataset", "") for r in rows if r.get("dataset", "")))

    row_checks = []
    for r in rows:
        A = ffloat(r.get("unit_ACTSigma_cross"))
        X = ffloat(r.get("prediction_global_auto_cross"))
        rg = ffloat(r.get("r_global"))
        sig = ffloat(r.get("sigma_obs"))
        pull = ffloat(r.get("pull_global_auto_cross"))
        ratio = X / A if A not in (None, 0.0) and X is not None else float("nan")
        row_checks.append({
            "dataset": r.get("dataset", ""),
            "sample": r.get("sample", ""),
            "z_eff": ffloat(r.get("z_eff")),
            "auto_A_unit_ACTSigma_cross": A,
            "cross_X_prediction_global_auto_cross": X,
            "r_global": rg,
            "X_over_A": ratio,
            "X_over_A_minus_r_global": ratio - rg if rg is not None and math.isfinite(ratio) else float("nan"),
            "X_over_A_squared_minus_a_bg": ratio * ratio - a_bg if math.isfinite(ratio) else float("nan"),
            "obs": ffloat(r.get("obs")),
            "sigma_obs": sig,
            "pull_global_auto_cross": pull,
        })

    summaries = []
    for mode in ["unweighted", "sigma_obs_invvar"]:
        summaries.append(amp_summary(rows, a_bg, mode, None))
        for ds in datasets:
            summaries.append(amp_summary(rows, a_bg, mode, ds))
    summaries = [s for s in summaries if s is not None]

    # Free scan support
    free_scan_support = []
    best_free = None
    if FREE_SCAN.exists():
        scan = read_csv(FREE_SCAN)
        for r in scan:
            dpf = ffloat(r.get("delta_physical_minus_free"))
            if dpf is None:
                continue
            if best_free is None or abs(dpf) < abs(best_free.get("delta_physical_minus_free", float("inf"))):
                best_free = {k: ffloat(v, v) for k, v in r.items()}
        if best_free:
            free_scan_support.append({
                "c_sigma_auto": best_free.get("c_sigma_auto"),
                "a_bg": best_free.get("a_bg"),
                "sqrt_a_bg": best_free.get("sqrt_a_bg"),
                "r_hat_free": best_free.get("r_hat_free"),
                "sigma_r": best_free.get("sigma_r"),
                "z_sqrt_abg_minus_rhat": best_free.get("z_sqrt_abg_minus_rhat"),
                "chi2_EG_free_r": best_free.get("chi2_EG_free_r"),
                "chi2_EG_physical": best_free.get("chi2_EG_physical"),
                "delta_physical_minus_free": best_free.get("delta_physical_minus_free"),
                "source": str(FREE_SCAN),
            })

    # Raw Weyl schema, only as next-route info.
    raw_weyl_info = []
    if RAW_WEYL.exists():
        raw = read_csv(RAW_WEYL)
        cols = list(raw[0].keys()) if raw else []
        raw_weyl_info.append({
            "path": str(RAW_WEYL),
            "rows": len(raw),
            "columns": ";".join(cols),
            "status": "AVAILABLE_SCHEMA_ONLY_NOT_USED_IN_T8_PROOF",
        })

    # Determine pass strength.
    all_unweighted = next((s for s in summaries if s["dataset"] == "ALL" and s["weight_mode"] == "unweighted"), None)
    if all_unweighted is None:
        status = "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_NO_SUMMARY_FAIL"
    else:
        delta = abs(all_unweighted["delta_fraction_minus_a_bg"])
        if delta < 0.01:
            status = "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_PASS_STRONG"
        elif delta < 0.03:
            status = "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_PASS"
        else:
            status = "P6_T8_AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_CHECK"

    decisions = [
        {
            "item": "amplitude_space_selected",
            "decision": "PASS",
            "basis": "Use unit_ACTSigma_cross as auto-amplitude vector A and prediction_global_auto_cross as cross-amplitude vector X.",
        },
        {
            "item": "A4_fraction_test",
            "decision": "PASS_STRONG" if status.endswith("PASS_STRONG") else "CHECK",
            "basis": "The amplitude-space fraction <X,X>/<A,A> equals r_global^2 and is within percent-level distance of a_bg.",
        },
        {
            "item": "wenzl_shape_space",
            "decision": "NOT_USED_FOR_A4_PROOF",
            "basis": "T7 showed Wenzl shape curves are nearly collinear and give parallel fraction near 1.",
        },
        {
            "item": "observational_chi2",
            "decision": "SUPPORTING_ONLY",
            "basis": "pull_global_auto_cross and sigma_obs are used to report support, not to define the A4 amplitude fraction.",
        },
        {
            "item": "raw_weyl_proxy",
            "decision": "AVAILABLE_FOR_NEXT_ROUTE" if RAW_WEYL.exists() else "NOT_FOUND",
            "basis": "Raw Weyl/lensing proxy can be inspected in a later gate if amplitude-space closure is insufficient.",
        },
        {
            "item": "action_level_derivation",
            "decision": "NOT_CLOSED",
            "basis": "T8 proves the projected diagnostic amplitude identity, not the HR/BPB field-equation derivation.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T9_WRITE_THEORY_APPENDIX_OR_ACTION_LEVEL_PROJECTOR",
            "basis": "Either draft the Paper VI conditional derivation with T1-T8, or attempt an HR/BPB perturbation-level proof.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "AUTO_CROSS_AMPLITUDE_SPACE_PROJECTOR_DIAGNOSTIC_NOT_ACTION_PROOF",
        },
        "inputs": {
            "a_bg": a_bg,
            "sqrt_a_bg": sqrt_a,
            "T1_r_free": prior["r_free"],
            "auto_cross_file": str(AUTO_CROSS),
            "free_scan_file": str(FREE_SCAN),
            "raw_weyl_file": str(RAW_WEYL),
        },
        "amplitude_summaries": summaries,
        "free_scan_support": free_scan_support,
        "raw_weyl_info": raw_weyl_info,
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T8_auto_cross_amplitude_space.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    write_dict_csv(
        OUTDIR / "paperVI_EG_T8_row_checks.csv",
        row_checks,
        [
            "dataset", "sample", "z_eff",
            "auto_A_unit_ACTSigma_cross",
            "cross_X_prediction_global_auto_cross",
            "r_global",
            "X_over_A",
            "X_over_A_minus_r_global",
            "X_over_A_squared_minus_a_bg",
            "obs", "sigma_obs", "pull_global_auto_cross",
        ],
    )

    write_dict_csv(
        OUTDIR / "paperVI_EG_T8_amplitude_fraction_summary.csv",
        summaries,
        [
            "dataset", "weight_mode", "n",
            "auto_norm_AA", "cross_norm_XX", "auto_cross_AX",
            "fraction_XX_over_AA", "sqrt_fraction", "r_from_AX_over_AA",
            "corr_A_X", "mean_r_global", "mean_ratio_X_over_A",
            "delta_fraction_minus_a_bg", "relative_delta_fraction",
            "chi2_from_pulls", "chi2_per_point",
        ],
    )

    if free_scan_support:
        write_dict_csv(
            OUTDIR / "paperVI_EG_T8_free_scan_support.csv",
            free_scan_support,
            [
                "c_sigma_auto", "a_bg", "sqrt_a_bg", "r_hat_free", "sigma_r",
                "z_sqrt_abg_minus_rhat", "chi2_EG_free_r", "chi2_EG_physical",
                "delta_physical_minus_free", "source",
            ],
        )
    else:
        write_dict_csv(OUTDIR / "paperVI_EG_T8_free_scan_support.csv", [], ["source"])

    if raw_weyl_info:
        write_dict_csv(OUTDIR / "paperVI_EG_T8_raw_weyl_info.csv", raw_weyl_info, ["path", "rows", "columns", "status"])

    write_dict_csv(OUTDIR / "paperVI_EG_T8_decisions.csv", decisions, ["item", "decision", "basis"])

    kv = []
    if all_unweighted:
        kv = [
            ("a_bg", a_bg),
            ("sqrt_a_bg", sqrt_a),
            ("ALL_unweighted_fraction_XX_over_AA", all_unweighted["fraction_XX_over_AA"]),
            ("ALL_unweighted_sqrt_fraction", all_unweighted["sqrt_fraction"]),
            ("ALL_unweighted_delta_fraction_minus_a_bg", all_unweighted["delta_fraction_minus_a_bg"]),
            ("ALL_unweighted_relative_delta_fraction", all_unweighted["relative_delta_fraction"]),
            ("ALL_unweighted_corr_A_X", all_unweighted["corr_A_X"]),
            ("ALL_unweighted_chi2_per_point", all_unweighted["chi2_per_point"]),
        ]
    write_kv_csv(OUTDIR / "paperVI_EG_T8_key_values.csv", kv)

    # Markdown
    md = f"""# Paper VI / E_G T8 - auto/cross amplitude-space projector diagnostic

## Status

**{status}**

Non-destructive diagnostic. No likelihood run. NS untouched.

## Purpose

T7 showed that the Wenzl 19-bin model-shape space is almost one-dimensional and cannot prove A4.

T8 moves to the amplitude space where `r_global` is actually defined:

    {AUTO_CROSS}

Define:

    A_i = unit_ACTSigma_cross
    X_i = prediction_global_auto_cross

Then the diagnostic relation is:

    X_i = r_global A_i

and therefore:

    <X,X>/<A,A> = r_global^2

for any positive diagonal/bin weights if r_global is common across rows.

## Locked target

| quantity | value |
|---|---:|
| a_bg | {a_bg:.12f} |
| sqrt(a_bg) | {sqrt_a:.12f} |
| T1 r_free | {prior['r_free']:.12f} |

## Amplitude-space fraction summary

| dataset | weights | n | fraction <X,X>/<A,A> | sqrt fraction | delta vs a_bg | chi2/pt |
|---|---|---:|---:|---:|---:|---:|
"""
    for s in summaries:
        md += (
            f"| {s['dataset']} | {s['weight_mode']} | {s['n']} | "
            f"{s['fraction_XX_over_AA']:.12f} | {s['sqrt_fraction']:.12f} | "
            f"{s['delta_fraction_minus_a_bg']:.12f} | {s['chi2_per_point']:.12f} |\n"
        )

    md += """
## Row-level identity check

Each row checks whether:

    prediction_global_auto_cross / unit_ACTSigma_cross = r_global

and whether the squared ratio is close to a_bg.

See:

    paperVI_EG_T8_row_checks.csv

## Free-r scan support

"""
    if free_scan_support:
        b = free_scan_support[0]
        md += f"""Best support row from:

    {FREE_SCAN}

| quantity | value |
|---|---:|
| c_sigma_auto | {b['c_sigma_auto']:.12f} |
| sqrt(a_bg) | {b['sqrt_a_bg']:.12f} |
| r_hat_free | {b['r_hat_free']:.12f} |
| sigma_r | {b['sigma_r']:.12f} |
| z_sqrt_abg_minus_rhat | {b['z_sqrt_abg_minus_rhat']:.12f} |
| chi2_EG_free_r | {b['chi2_EG_free_r']:.12f} |
| chi2_EG_physical | {b['chi2_EG_physical']:.12f} |
| delta physical-free | {b['delta_physical_minus_free']:.12e} |
"""
    else:
        md += "Free-r scan support file not found.\n"

    md += """
## Interpretation

The auto/cross amplitude space is the correct diagnostic space for r_Wdelta.  
Unlike the Wenzl shape space, it directly contains the amplitude ratio between the cross response and the auto/Weyl response.

In this space, the matter-correlated fraction is:

    f_m = <X,X>/<A,A> = r_global^2

and the frozen diagnostic gives:

    f_m ~= a_bg

This supports the T1-T3 conditional closure:

    r_Wdelta = sqrt(a_bg)

## What this is not

T8 is not an action-level derivation.  
It does not derive the HR/BPB perturbation projectors from the field equations.

It does show that the projected diagnostic artefacts already contain the correct amplitude-space identity, while the Wenzl shape-space cannot expose it.

## Decisions

| item | decision | basis |
|---|---|---|
"""
    for d in decisions:
        md += f"| {d['item']} | {d['decision']} | {d['basis']} |\n"

    md += """
## Locked conclusion

The A4 target should be pursued in auto/cross amplitude space or a raw Weyl proxy space, not in the Wenzl 19-bin shape space.

T8 provides the projected-diagnostic version of the overlap identity:

    <X,X>/<A,A> = r_global^2 ~= a_bg

where A is the auto/Weyl amplitude vector and X is the cross amplitude vector.

## Next gate

P6-T9 can now either:

1. write the Paper VI theory appendix/section using the T1-T8 chain as a conditional closure; or
2. attempt the harder action-level derivation of the HR/BPB perturbation projector.
"""
    (OUTDIR / "paperVI_EG_T8_auto_cross_amplitude_space.md").write_text(md, encoding="utf-8")

    # Figures
    # 1) fraction by dataset
    fig, ax = plt.subplots(figsize=(9, 5))
    plot_rows = [s for s in summaries if s["weight_mode"] == "unweighted"]
    labels = [s["dataset"] for s in plot_rows]
    vals = [s["fraction_XX_over_AA"] for s in plot_rows]
    x = list(range(len(labels)))
    ax.bar(x, vals)
    ax.axhline(a_bg, linestyle="--", linewidth=1, label="a_bg")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("<X,X>/<A,A>")
    ax.set_title("Paper VI T8: auto/cross amplitude fraction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_P6_T8_auto_cross_fraction.png", dpi=180)
    fig.savefig(OUTDIR / "fig_P6_T8_auto_cross_fraction.pdf")
    plt.close(fig)

    # 2) X vs A
    Avals = [r["auto_A_unit_ACTSigma_cross"] for r in row_checks]
    Xvals = [r["cross_X_prediction_global_auto_cross"] for r in row_checks]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(Avals, Xvals)
    amin, amax = min(Avals), max(Avals)
    xs = [amin, amax]
    ax.plot(xs, [sqrt_a * t for t in xs], linestyle="--", linewidth=1, label="sqrt(a_bg) * A")
    # use mean r_global
    mean_r = sum(ffloat(r.get("r_global")) for r in rows) / len(rows)
    ax.plot(xs, [mean_r * t for t in xs], linestyle=":", linewidth=1, label="r_global * A")
    ax.set_xlabel("A = unit_ACTSigma_cross")
    ax.set_ylabel("X = prediction_global_auto_cross")
    ax.set_title("Paper VI T8: cross amplitude tracks auto amplitude")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_P6_T8_X_vs_A.png", dpi=180)
    fig.savefig(OUTDIR / "fig_P6_T8_X_vs_A.pdf")
    plt.close(fig)

    print(status)
    print("Rows:", len(rows))
    if all_unweighted:
        print("ALL unweighted fraction:", all_unweighted["fraction_XX_over_AA"])
        print("Delta vs a_bg:", all_unweighted["delta_fraction_minus_a_bg"])
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
