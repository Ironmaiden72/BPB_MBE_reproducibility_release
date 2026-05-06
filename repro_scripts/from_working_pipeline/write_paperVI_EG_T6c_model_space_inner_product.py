#!/usr/bin/env python3
"""
Paper VI / E_G T6c — select model-space inner product.

Purpose
-------
T6/T6b established that the primary Wenzl projected file is a theory/model
projection space, not a data+covariance space. T6c therefore defines a
model-space inner product over the projected bins for Weyl-response vectors.

Important guardrail
-------------------
This is NOT an observational chi2 and must not be described as one.
The columns EG_bpb_act_sigma_projected and EG_bpb_eg_sigma_projected are
model responses, not observational error bars.

Writes only:
    outputs/paperVI_EG_theory/

Does not touch NS files and runs no likelihood.
"""
from __future__ import annotations

from pathlib import Path
import csv
import json
import math
from typing import Dict, List, Any, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTDIR = Path("outputs/paperVI_EG_theory")

PRIMARY_PROJECTION = Path("outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.csv")
SUPPORTING_DATA_SPACE = Path("outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv")

MODEL_COLUMNS = {
    "W_LCDM": "EG_lcdm_projected",
    "W_mu": "EG_bpb_mu_projected",
    "W_ACTSigma": "EG_bpb_act_sigma_projected",
    "W_EGSigma": "EG_bpb_eg_sigma_projected",
}

AXIS_COLUMNS = ["sample", "ell_lo", "ell_hi", "ell_eff", "k_eff"]


def ffloat(x: Any) -> float:
    if x is None:
        return float("nan")
    s = str(x).strip()
    if s == "":
        return float("nan")
    try:
        return float(s)
    except Exception:
        return float("nan")


def read_csv_dict(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        return reader.fieldnames or [], list(reader)


def dot(a: List[float], b: List[float], weights: List[float]) -> float:
    return sum(w * x * y for x, y, w in zip(a, b, weights))


def norm(a: List[float], weights: List[float]) -> float:
    v = dot(a, a, weights)
    return math.sqrt(v) if v >= 0 else float("nan")


def corr(a: List[float], b: List[float], weights: List[float]) -> float:
    na = norm(a, weights)
    nb = norm(b, weights)
    if not math.isfinite(na) or not math.isfinite(nb) or na == 0 or nb == 0:
        return float("nan")
    return dot(a, b, weights) / (na * nb)


def sub(a: List[float], b: List[float]) -> List[float]:
    return [x-y for x, y in zip(a, b)]


def unit(a: List[float], weights: List[float]) -> List[float]:
    n = norm(a, weights)
    if not math.isfinite(n) or n == 0:
        return [float("nan") for _ in a]
    return [x/n for x in a]


def projection_coeff(a: List[float], basis: List[float], weights: List[float]) -> float:
    denom = dot(basis, basis, weights)
    if denom == 0:
        return float("nan")
    return dot(a, basis, weights) / denom


def write_dict_csv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def write_matrix_csv(path: Path, names: List[str], mat: Dict[Tuple[str, str], float]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([""] + names)
        for r in names:
            w.writerow([r] + [mat[(r, c)] for c in names])


def build_supporting_data_diagnostic() -> Dict[str, Any]:
    if not SUPPORTING_DATA_SPACE.exists():
        return {"exists": False, "path": str(SUPPORTING_DATA_SPACE)}

    cols, rows = read_csv_dict(SUPPORTING_DATA_SPACE)
    required = ["dataset", "sample", "obs", "sigma_obs", "prediction_global_auto_cross", "r_global"]
    missing = [c for c in required if c not in cols]
    if missing:
        return {"exists": True, "path": str(SUPPORTING_DATA_SPACE), "missing": missing}

    out_rows = []
    by_dataset: Dict[str, Dict[str, float]] = {}
    r_values = []

    for r in rows:
        obs = ffloat(r.get("obs"))
        sig = ffloat(r.get("sigma_obs"))
        pred = ffloat(r.get("prediction_global_auto_cross"))
        rg = ffloat(r.get("r_global"))
        if math.isfinite(rg):
            r_values.append(rg)
        if not (math.isfinite(obs) and math.isfinite(sig) and sig > 0 and math.isfinite(pred)):
            continue
        pull = (obs - pred) / sig
        chi2 = pull * pull
        ds = str(r.get("dataset", "UNKNOWN"))
        by_dataset.setdefault(ds, {"chi2": 0.0, "n": 0})
        by_dataset[ds]["chi2"] += chi2
        by_dataset[ds]["n"] += 1
        out_rows.append({
            "dataset": ds,
            "sample": r.get("sample", ""),
            "obs": obs,
            "sigma_obs": sig,
            "prediction_global_auto_cross": pred,
            "pull": pull,
            "chi2_contribution": chi2,
            "r_global": rg,
        })

    summary_rows = []
    for ds, vals in by_dataset.items():
        summary_rows.append({
            "dataset": ds,
            "n": int(vals["n"]),
            "chi2_diag": vals["chi2"],
            "chi2_per_point": vals["chi2"]/vals["n"] if vals["n"] else float("nan"),
        })

    return {
        "exists": True,
        "path": str(SUPPORTING_DATA_SPACE),
        "columns": cols,
        "rows": out_rows,
        "summary": summary_rows,
        "r_global_min": min(r_values) if r_values else None,
        "r_global_max": max(r_values) if r_values else None,
        "r_global_mean": sum(r_values)/len(r_values) if r_values else None,
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    if not PRIMARY_PROJECTION.exists():
        status = "P6_T6C_PRIMARY_PROJECTION_FILE_MISSING_FAIL"
        result = {
            "status": {
                "gate_status": status,
                "likelihood_run": "NO",
                "ns_touched": "NO",
            },
            "missing_file": str(PRIMARY_PROJECTION),
        }
        (OUTDIR / "paperVI_EG_T6c_model_space_inner_product.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(status)
        print("Missing:", PRIMARY_PROJECTION)
        return

    cols, rows = read_csv_dict(PRIMARY_PROJECTION)
    missing_model_cols = [c for c in MODEL_COLUMNS.values() if c not in cols]
    missing_axis_cols = [c for c in AXIS_COLUMNS if c not in cols]

    if missing_model_cols:
        status = "P6_T6C_MODEL_COLUMNS_MISSING_FAIL"
        result = {
            "status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO"},
            "missing_model_cols": missing_model_cols,
            "columns": cols,
        }
        (OUTDIR / "paperVI_EG_T6c_model_space_inner_product.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(status)
        print("Missing model columns:", missing_model_cols)
        return

    # Build vectors.
    labels = []
    ell_widths = []
    vectors: Dict[str, List[float]] = {name: [] for name in MODEL_COLUMNS}

    for i, r in enumerate(rows):
        sample = r.get("sample", f"row{i}")
        ell_lo = ffloat(r.get("ell_lo"))
        ell_hi = ffloat(r.get("ell_hi"))
        ell_eff = ffloat(r.get("ell_eff"))
        labels.append(f"{sample}:{ell_eff:.1f}" if math.isfinite(ell_eff) else f"{sample}:row{i}")
        width = ell_hi - ell_lo if math.isfinite(ell_hi) and math.isfinite(ell_lo) and ell_hi > ell_lo else 1.0
        ell_widths.append(width)
        for name, col in MODEL_COLUMNS.items():
            vectors[name].append(ffloat(r.get(col)))

    n = len(rows)
    unweighted = [1.0 / n for _ in rows]
    sw = sum(ell_widths)
    ell_weighted = [w / sw for w in ell_widths] if sw > 0 else unweighted

    weight_schemes = {
        "unweighted_mean": unweighted,
        "ell_width_weighted_mean": ell_weighted,
    }

    metric_rows = []
    corr_rows = []
    vector_norm_rows = []
    residual_rows = []

    names = list(MODEL_COLUMNS.keys())

    corr_mats = {}
    dot_mats = {}

    for scheme, weights in weight_schemes.items():
        corr_mat = {}
        dot_mat = {}
        for a in names:
            for b in names:
                dot_mat[(a, b)] = dot(vectors[a], vectors[b], weights)
                corr_mat[(a, b)] = corr(vectors[a], vectors[b], weights)
        corr_mats[scheme] = corr_mat
        dot_mats[scheme] = dot_mat

        for a in names:
            vector_norm_rows.append({
                "scheme": scheme,
                "vector": a,
                "column": MODEL_COLUMNS[a],
                "norm": norm(vectors[a], weights),
                "mean": sum(w*x for w, x in zip(weights, vectors[a])),
                "min": min(vectors[a]),
                "max": max(vectors[a]),
            })

        for a in names:
            for b in names:
                metric_rows.append({
                    "scheme": scheme,
                    "A": a,
                    "B": b,
                    "inner_product": dot_mat[(a, b)],
                    "correlation": corr_mat[(a, b)],
                })

        # Useful differences relative to matter-like vector W_mu.
        for target in ["W_ACTSigma", "W_EGSigma", "W_LCDM"]:
            dvec = sub(vectors[target], vectors["W_mu"])
            residual_rows.append({
                "scheme": scheme,
                "residual": f"{target}-W_mu",
                "norm": norm(dvec, weights),
                "corr_with_W_mu": corr(dvec, vectors["W_mu"], weights),
                "projection_coeff_on_W_mu": projection_coeff(dvec, vectors["W_mu"], weights),
            })

        # Differences between Sigma branches.
        dvec = sub(vectors["W_ACTSigma"], vectors["W_EGSigma"])
        residual_rows.append({
            "scheme": scheme,
            "residual": "W_ACTSigma-W_EGSigma",
            "norm": norm(dvec, weights),
            "corr_with_W_mu": corr(dvec, vectors["W_mu"], weights),
            "projection_coeff_on_W_mu": projection_coeff(dvec, vectors["W_mu"], weights),
        })

    # Canonical decision: use model-space, unweighted first, ell-width as robustness.
    support = build_supporting_data_diagnostic()

    support_summary = support.get("summary", [])
    if support.get("rows"):
        write_dict_csv(
            OUTDIR / "paperVI_EG_T6c_supporting_autocross_rows.csv",
            support["rows"],
            ["dataset", "sample", "obs", "sigma_obs", "prediction_global_auto_cross", "pull", "chi2_contribution", "r_global"],
        )
    if support_summary:
        write_dict_csv(
            OUTDIR / "paperVI_EG_T6c_supporting_autocross_summary.csv",
            support_summary,
            ["dataset", "n", "chi2_diag", "chi2_per_point"],
        )

    decisions = [
        {
            "item": "route",
            "decision": "MODEL_SPACE_SELECTED",
            "basis": "Primary Wenzl file is theory-space; for Weyl projectors we need a response-vector space, not a data chi2.",
        },
        {
            "item": "canonical_inner_product",
            "decision": "UNWEIGHTED_BIN_MEAN_WITH_ELL_WIDTH_ROBUSTNESS",
            "basis": "No covariance/window weights are present in the primary file; unweighted mean is transparent, ell-width weighted is reported as robustness.",
        },
        {
            "item": "not_observational_chi2",
            "decision": "LOCKED_GUARDRAIL",
            "basis": "This norm is a model-space norm over projected bins, not A^T C^-1 B for observed data.",
        },
        {
            "item": "canonical_vectors",
            "decision": "W_MU_AS_MATTER_LIKE_DIRECTION_AND_SIGMA_COLUMNS_AS_WEYL_CANDIDATES",
            "basis": "W_mu is the projected matter/growth-like BPB response; ACTSigma/EGSigma are projected Weyl/Sigma-like responses.",
        },
        {
            "item": "supporting_data_space",
            "decision": "AUTO_CROSS_GLOBAL_R_USED_AS_SUPPORTING_DIAGNOSTIC",
            "basis": "eg_projected_auto_cross_global_r_seed01.csv contains obs, sigma_obs, r_global and prediction_global_auto_cross.",
        },
        {
            "item": "P_m_projector",
            "decision": "NOT_BUILT_YET",
            "basis": "T6c defines the model-space norm and vectors; T7 can construct/project candidate matter/private directions.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T7_MODEL_SPACE_PROJECTOR_TRIAL",
            "basis": "Use W_mu as candidate e_m and decompose Weyl/Sigma responses into parallel/orthogonal components.",
        },
    ]

    status = "P6_T6C_MODEL_SPACE_INNER_PRODUCT_SELECTED_AND_BUILT"

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "MODEL_SPACE_INNER_PRODUCT_BUILT_NOT_DATA_CHI2_NOT_PROJECTOR_PROOF",
        },
        "primary_projection": {
            "path": str(PRIMARY_PROJECTION),
            "n_rows": n,
            "columns": cols,
            "missing_axis_cols": missing_axis_cols,
        },
        "model_columns": MODEL_COLUMNS,
        "weight_schemes": {
            "unweighted_mean": "weights_i = 1/N",
            "ell_width_weighted_mean": "weights_i = (ell_hi-ell_lo)/sum_j(ell_hi-ell_lo)",
        },
        "supporting_data_space": {
            "path": str(SUPPORTING_DATA_SPACE),
            "exists": support.get("exists", False),
            "r_global_min": support.get("r_global_min"),
            "r_global_max": support.get("r_global_max"),
            "r_global_mean": support.get("r_global_mean"),
            "summary": support_summary,
        },
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T6c_model_space_inner_product.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    # Write CSV outputs.
    write_dict_csv(
        OUTDIR / "paperVI_EG_T6c_model_vectors.csv",
        [
            {
                "row": i,
                "label": labels[i],
                "sample": rows[i].get("sample", ""),
                "ell_lo": rows[i].get("ell_lo", ""),
                "ell_hi": rows[i].get("ell_hi", ""),
                "ell_eff": rows[i].get("ell_eff", ""),
                "k_eff": rows[i].get("k_eff", ""),
                **{name: vectors[name][i] for name in names},
                "weight_unweighted": unweighted[i],
                "weight_ell_width": ell_weighted[i],
            }
            for i in range(n)
        ],
        ["row", "label", "sample", "ell_lo", "ell_hi", "ell_eff", "k_eff"] + names + ["weight_unweighted", "weight_ell_width"],
    )

    write_dict_csv(
        OUTDIR / "paperVI_EG_T6c_inner_products.csv",
        metric_rows,
        ["scheme", "A", "B", "inner_product", "correlation"],
    )

    write_dict_csv(
        OUTDIR / "paperVI_EG_T6c_vector_norms.csv",
        vector_norm_rows,
        ["scheme", "vector", "column", "norm", "mean", "min", "max"],
    )

    write_dict_csv(
        OUTDIR / "paperVI_EG_T6c_residual_directions.csv",
        residual_rows,
        ["scheme", "residual", "norm", "corr_with_W_mu", "projection_coeff_on_W_mu"],
    )

    write_dict_csv(
        OUTDIR / "paperVI_EG_T6c_decisions.csv",
        decisions,
        ["item", "decision", "basis"],
    )

    for scheme, mat in corr_mats.items():
        write_matrix_csv(OUTDIR / f"paperVI_EG_T6c_corr_matrix_{scheme}.csv", names, mat)
        write_matrix_csv(OUTDIR / f"paperVI_EG_T6c_inner_matrix_{scheme}.csv", names, dot_mats[scheme])

    # Markdown report.
    un_corr = corr_mats["unweighted_mean"]
    ell_corr = corr_mats["ell_width_weighted_mean"]

    md = f"""# Paper VI / E_G T6c - model-space inner product selection

## Status

**{status}**

Non-destructive gate. No likelihood run. NS untouched.

## Why model-space route is selected

T6/T6b showed that the primary Wenzl projected file is a projected theory/model space, not an observed data/covariance space:

    {PRIMARY_PROJECTION}

Therefore T6c defines a model-space norm over the projected bins. This is used only to compare/decompose Weyl response vectors. It is **not** an observational chi2.

## Canonical model-space vectors

| symbol | column | interpretation |
|---|---|---|
| W_LCDM | EG_lcdm_projected | LCDM projected E_G response |
| W_mu | EG_bpb_mu_projected | BPB matter/growth-like projected response |
| W_ACTSigma | EG_bpb_act_sigma_projected | BPB ACT/Sigma-like Weyl response |
| W_EGSigma | EG_bpb_eg_sigma_projected | BPB E_G/Sigma-like Weyl response |

## Inner product definition

Canonical model-space norm:

    <A,B>_model = (1/N) sum_i A_i B_i

Robustness norm:

    <A,B>_ell = sum_i w_i A_i B_i
    w_i = (ell_hi_i - ell_lo_i) / sum_j (ell_hi_j - ell_lo_j)

Neither is a data chi2. They are projected-bin response-vector norms.

## Unweighted correlations

| A | B | corr |
|---|---|---:|
"""
    for a in names:
        for b in names:
            if a < b:
                md += f"| {a} | {b} | {un_corr[(a,b)]:.12f} |\n"

    md += """
## Ell-width weighted correlations

| A | B | corr |
|---|---|---:|
"""
    for a in names:
        for b in names:
            if a < b:
                md += f"| {a} | {b} | {ell_corr[(a,b)]:.12f} |\n"

    md += """
## Supporting data-space diagnostic

T6c also audits the supporting auto/cross file:

    outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv

This file contains observed values and errors (`obs`, `sigma_obs`) and the fitted/global r diagnostic. It remains supporting evidence for r_Wdelta, while the Wenzl projected file supplies the model-space response vectors.
"""

    if support_summary:
        md += """
| dataset | n | chi2_diag | chi2/pt |
|---|---:|---:|---:|
"""
        for row in support_summary:
            md += f"| {row['dataset']} | {row['n']} | {row['chi2_diag']:.12f} | {row['chi2_per_point']:.12f} |\n"

    md += f"""
## Guardrail

The column `EG_bpb_act_sigma_projected` is not an observational sigma/error column. It is a projected model response. It must not be used as sigma_i in a data-space chi2.

## Decisions

| item | decision | basis |
|---|---|---|
"""
    for row in decisions:
        md += f"| {row['item']} | {row['decision']} | {row['basis']} |\n"

    md += """
## Locked conclusion

T6c selects the model-space route and builds a transparent projected-bin norm for the Wenzl response vectors. This is the correct space for preparing Weyl/matter projectors.

It does not prove A4 and does not build P_m. The next gate can now perform a first projector trial, using W_mu as a candidate matter-like direction and Sigma responses as Weyl candidates.

## Next gate

P6-T7 = model-space projector trial:

    e_m = normalized W_mu
    decompose W_ACTSigma and W_EGSigma into components parallel and orthogonal to e_m
    compare the squared overlap fraction with a_bg
"""

    (OUTDIR / "paperVI_EG_T6c_model_space_inner_product.md").write_text(md, encoding="utf-8")

    # Figure.
    fig, ax = plt.subplots(figsize=(9, 5))
    x = list(range(n))
    for name in names:
        ax.plot(x, vectors[name], marker="o", label=name)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("projected E_G response")
    ax.set_title("Paper VI T6c: model-space projected response vectors")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_P6_T6c_model_space_vectors.png", dpi=180)
    fig.savefig(OUTDIR / "fig_P6_T6c_model_space_vectors.pdf")
    plt.close(fig)

    print(status)
    print("Primary:", PRIMARY_PROJECTION)
    print("Rows:", n)
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
