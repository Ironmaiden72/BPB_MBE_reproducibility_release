#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_SETS = [
    {
        "name": "wenzl_act_dr6_planck_pr4_boss_final",
        "rank": 1,
        "csv": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.csv",
        "json": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.json",
        "role": "PRIMARY_CANONICAL_CANDIDATE",
        "reason": "Highest T5 score; includes ACT DR6 + Planck PR4 + BOSS projected E_G space.",
    },
    {
        "name": "wenzl_act_dr6_planck_pr4_boss_frozen_point",
        "rank": 2,
        "csv": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_frozen_point.csv",
        "json": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_frozen_point.json",
        "role": "SUPPORTING_FROZEN_POINT",
        "reason": "Same projected family, frozen point variant.",
    },
    {
        "name": "wenzl_act_dr6_planck_pr4_boss_seed01",
        "rank": 3,
        "csv": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_seed01.csv",
        "json": "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_seed01.json",
        "role": "SUPPORTING_SEED01",
        "reason": "Seed01 diagnostic variant.",
    },
    {
        "name": "wenzl_planck_pr4_boss_final",
        "rank": 4,
        "csv": "outputs/cosmology/eg_projected_wenzl_planck_pr4_boss_final.csv",
        "json": "outputs/cosmology/eg_projected_wenzl_planck_pr4_boss_final.json",
        "role": "SUPPORTING_PLANCK_PR4_ONLY",
        "reason": "Planck PR4 + BOSS projected E_G space.",
    },
    {
        "name": "auto_cross_global_r_seed01",
        "rank": 5,
        "csv": "outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv",
        "json": "outputs/cosmology/eg_projected_auto_cross_global_r_seed01.json",
        "role": "SUPPORTING_RFIT_DIAGNOSTIC",
        "reason": "Auto/cross global-r diagnostic relevant to r_Wdelta.",
    },
    {
        "name": "autocross_csigma_free_r_scan_final",
        "rank": 6,
        "csv": "outputs/cosmology/autocross_csigma_free_r_scan_final.csv",
        "json": "outputs/cosmology/autocross_csigma_free_r_scan_final.json",
        "role": "SUPPORTING_FREE_R_SCAN",
        "reason": "Free-r scan diagnostic.",
    },
    {
        "name": "cmb_lensing_weyl_proxy_planck_pr4_raw",
        "rank": 7,
        "csv": "outputs/cmb/cmb_lensing_weyl_proxy_map1782_Planck_PR4_raw.csv",
        "json": "",
        "role": "SUPPORTING_WEYL_PROXY",
        "reason": "Raw Weyl/lensing proxy candidate, not primary E_G compressed space.",
    },
]


def parse_float(x):
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null"}:
        return None
    try:
        return float(s)
    except Exception:
        return None


def flatten_keys(obj, prefix="", out=None, max_keys=500):
    if out is None:
        out = []
    if len(out) >= max_keys:
        return out
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.append(key)
            flatten_keys(v, key, out, max_keys)
            if len(out) >= max_keys:
                break
    elif isinstance(obj, list) and obj:
        key = f"{prefix}[0]" if prefix else "[0]"
        out.append(key)
        flatten_keys(obj[0], key, out, max_keys)
    return out


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        rows = list(reader)
        columns = reader.fieldnames or []
    return columns, rows


def numeric_profile(columns, rows):
    prof = []
    n = max(1, len(rows))
    for col in columns:
        vals = [parse_float(r.get(col)) for r in rows]
        nums = [v for v in vals if v is not None and math.isfinite(v)]
        nn = len(nums)
        prof.append({
            "column": col,
            "numeric_count": nn,
            "n_rows": len(rows),
            "numeric_fraction": nn / n,
            "min": min(nums) if nums else None,
            "max": max(nums) if nums else None,
        })
    return prof


def score_column(col, role):
    c = col.lower()
    if role == "bin":
        tokens = ["bin", "dataset", "sample", "z", "zeff", "z_eff", "lowz", "cmass", "label"]
        return sum(t in c for t in tokens)
    if role == "data":
        positives = ["obs", "data", "meas", "eg_obs", "e_g_obs", "eg_data", "value"]
        negatives = ["model", "pred", "theory", "bpb", "lcdm", "chi2", "sigma", "err", "cov"]
        return 3 * sum(t in c for t in positives) - 4 * sum(t in c for t in negatives)
    if role == "model":
        positives = ["model", "pred", "theory", "bpb", "map1782", "strict", "cs_surv", "wenzl"]
        negatives = ["obs", "data", "meas", "sigma", "err", "cov", "chi2"]
        return 3 * sum(t in c for t in positives) - 4 * sum(t in c for t in negatives)
    if role == "error":
        positives = ["sigma", "err", "error", "unc", "std"]
        negatives = ["model", "pred", "theory", "obs", "data", "chi2"]
        return 3 * sum(t in c for t in positives) - 4 * sum(t in c for t in negatives)
    if role == "weyl":
        positives = ["weyl", "sigma", "phi", "lens", "kappa", "wdelta", "auto", "cross"]
        negatives = ["err", "error", "chi2"]
        return 3 * sum(t in c for t in positives) - 3 * sum(t in c for t in negatives)
    return 0


def select_column(columns, rows, role, exclude=None):
    exclude = set(exclude or [])
    prof = numeric_profile(columns, rows)
    candidates = []
    for item in prof:
        col = item["column"]
        if col in exclude:
            continue
        if role in {"data", "model", "error", "weyl"} and item["numeric_fraction"] < 0.7:
            continue
        sc = score_column(col, role)
        if sc <= 0:
            continue
        if role == "error":
            vals = [parse_float(r.get(col)) for r in rows]
            good = [v for v in vals if v is not None and math.isfinite(v)]
            if not good or any(v <= 0 for v in good):
                continue
        candidates.append((sc, item["numeric_fraction"], col, item))
    candidates.sort(key=lambda x: (-x[0], -x[1], x[2]))
    if not candidates:
        return None, []
    return candidates[0][2], candidates


def build_vectors(rows, data_col, model_col, err_col):
    d, m, sig, labels = [], [], [], []
    for i, r in enumerate(rows):
        dv = parse_float(r.get(data_col))
        mv = parse_float(r.get(model_col))
        sv = parse_float(r.get(err_col))
        if dv is None or mv is None or sv is None or sv <= 0:
            continue
        d.append(dv)
        m.append(mv)
        sig.append(sv)
        label = None
        for key in ["bin", "dataset", "sample", "label", "z", "zeff", "z_eff"]:
            if key in r and str(r.get(key)).strip():
                label = str(r.get(key)).strip()
                break
        labels.append(label or f"row_{i}")
    return labels, d, m, sig


def inner_product(a, b, sig):
    return sum((ai * bi) / (si * si) for ai, bi, si in zip(a, b, sig))


def write_dict_csv(path, rows, fields):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    candidate_rows = []
    selected = None
    for c in CANDIDATE_SETS:
        csv_path = Path(c["csv"]) if c["csv"] else None
        json_path = Path(c["json"]) if c["json"] else None
        csv_exists = bool(csv_path and csv_path.exists())
        json_exists = bool(json_path and json_path.exists())
        row = dict(c)
        row["csv_exists"] = csv_exists
        row["json_exists"] = json_exists
        row["selected"] = False
        candidate_rows.append(row)
        if selected is None and csv_exists:
            selected = row
    if selected is None:
        status = "P6_T6_NO_CANONICAL_CSV_FOUND_FAIL"
        result = {"status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO"}, "decisions": [{"item": "canonical_space", "decision": "FAIL", "basis": "No candidate CSV found from the T5-prioritized canonical list."}]}
        (OUTDIR / "paperVI_EG_T6_canonical_space.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(status)
        return
    selected["selected"] = True
    csv_path = Path(selected["csv"])
    json_path = Path(selected["json"]) if selected["json"] else None
    columns, rows = read_csv(csv_path)
    n_rows = len(rows)
    prof = numeric_profile(columns, rows)
    json_keys = []
    if json_path and json_path.exists():
        try:
            data = json.loads(json_path.read_text(encoding="utf-8", errors="replace"))
            json_keys = flatten_keys(data)
        except Exception:
            json_keys = []
    data_col, data_candidates = select_column(columns, rows, "data")
    model_col, model_candidates = select_column(columns, rows, "model", exclude={data_col} if data_col else set())
    err_col, err_candidates = select_column(columns, rows, "error", exclude={data_col, model_col} if model_col else {data_col} if data_col else set())
    weyl_col, weyl_candidates = select_column(columns, rows, "weyl", exclude={data_col, model_col, err_col})
    bin_col, bin_candidates = select_column(columns, rows, "bin")
    can_build = data_col is not None and model_col is not None and err_col is not None
    vector_rows = []
    inner = {}
    if can_build:
        labels, d, m, sig = build_vectors(rows, data_col, model_col, err_col)
        residual = [di - mi for di, mi in zip(d, m)]
        dd = inner_product(d, d, sig)
        mm = inner_product(m, m, sig)
        dm = inner_product(d, m, sig)
        rr = inner_product(residual, residual, sig)
        corr_dm = dm / math.sqrt(dd * mm) if dd > 0 and mm > 0 else None
        inner = {"n_vector": len(d), "data_col": data_col, "model_col": model_col, "error_col": err_col, "weyl_col_candidate": weyl_col, "bin_col_candidate": bin_col, "inner_data_data": dd, "inner_model_model": mm, "inner_data_model": dm, "chi2_model_vs_data_diag": rr, "corr_data_model_diag": corr_dm}
        for lab, di, mi, si, ri in zip(labels, d, m, sig, residual):
            vector_rows.append({"label": lab, "data": di, "model": mi, "sigma": si, "residual_data_minus_model": ri, "pull": ri / si})
        status = "P6_T6_CANONICAL_EG_SPACE_AND_DIAG_INNER_PRODUCT_BUILT" if len(d) >= 2 else "P6_T6_CANONICAL_SPACE_SELECTED_NEEDS_COLUMN_MAP"
    else:
        status = "P6_T6_CANONICAL_SPACE_SELECTED_NEEDS_COLUMN_MAP"
    selected_columns = [
        {"role": "data_vector_d", "column": data_col or "", "status": "FOUND" if data_col else "MISSING"},
        {"role": "model_vector_m", "column": model_col or "", "status": "FOUND" if model_col else "MISSING"},
        {"role": "diagonal_error_sigma", "column": err_col or "", "status": "FOUND" if err_col else "MISSING"},
        {"role": "weyl_response_W_candidate", "column": weyl_col or "", "status": "FOUND" if weyl_col else "MISSING_OR_NOT_IN_PRIMARY_CSV"},
        {"role": "bin_label_candidate", "column": bin_col or "", "status": "FOUND" if bin_col else "MISSING"},
    ]
    decisions = [
        {"item": "canonical_candidate", "decision": "SELECTED", "basis": f"{selected['name']} selected as first available highest-priority T5 candidate."},
        {"item": "projected_inner_product", "decision": "BUILT_DIAGONAL" if can_build else "NEEDS_COLUMN_MAP", "basis": "Built as sum A_i B_i/sigma_i^2 from detected data/model/error columns." if can_build else "Could not unambiguously detect data/model/error columns."},
        {"item": "canonical_W_vector", "decision": "CANDIDATE_ONLY" if weyl_col else "NOT_BUILT", "basis": "T6 may identify a Weyl-like column, but P_m/W decomposition is deferred to T7."},
        {"item": "P_m_projector", "decision": "NOT_ATTEMPTED", "basis": "T6 only builds the projected space; P_m requires matter/private basis construction."},
        {"item": "next_gate", "decision": "P6_T7_CANONICAL_W_VECTOR_AND_PROJECTOR_PREP" if can_build else "P6_T6B_MANUAL_COLUMN_MAP", "basis": "If the inner product is built, next select W and prepare P_m; otherwise map columns manually."},
    ]
    result = {"status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO", "theory_status": "CANONICAL_EG_SPACE_SELECTED_INNER_PRODUCT_BUILT_IF_COLUMNS_FOUND"}, "selected_candidate": selected, "csv": {"path": str(csv_path), "n_rows": n_rows, "columns": columns}, "json": {"path": str(json_path) if json_path else "", "n_keys_sample": len(json_keys), "keys_sample": json_keys[:160]}, "selected_columns": selected_columns, "inner_product": inner, "decisions": decisions}
    (OUTDIR / "paperVI_EG_T6_canonical_space.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    write_dict_csv(OUTDIR / "paperVI_EG_T6_candidate_selection.csv", candidate_rows, ["name", "rank", "csv", "json", "role", "reason", "csv_exists", "json_exists", "selected"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6_numeric_column_profile.csv", prof, ["column", "numeric_count", "n_rows", "numeric_fraction", "min", "max"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6_selected_columns.csv", selected_columns, ["role", "column", "status"])
    if vector_rows:
        write_dict_csv(OUTDIR / "paperVI_EG_T6_projected_vector_rows.csv", vector_rows, ["label", "data", "model", "sigma", "residual_data_minus_model", "pull"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6_decisions.csv", decisions, ["item", "decision", "basis"])
    all_cand_rows = []
    for role, cands in [("data", data_candidates), ("model", model_candidates), ("error", err_candidates), ("weyl", weyl_candidates), ("bin", bin_candidates)]:
        for sc, frac, col, item in cands[:10]:
            all_cand_rows.append({"role": role, "column": col, "score": sc, "numeric_fraction": frac, "min": item["min"], "max": item["max"]})
    write_dict_csv(OUTDIR / "paperVI_EG_T6_column_role_candidates.csv", all_cand_rows, ["role", "column", "score", "numeric_fraction", "min", "max"])
    md = f"""# Paper VI / E_G T6 - canonical projected E_G space and inner product\n\n## Status\n\n**{status}**\n\nNon-destructive gate. No likelihood run. NS untouched.\n\n## Purpose\n\nT5 found many E_G/lensing/Weyl artefacts.  \nT6 selects the first canonical projected E_G candidate and attempts to build the first explicit projected inner product:\n\n    <A,B>_EG = sum_i A_i B_i / sigma_i^2\n\nThis is not yet the matter projector P_m and not an action-level proof.\n\n## Selected canonical candidate\n\n| field | value |\n|---|---|\n| name | {selected['name']} |\n| role | {selected['role']} |\n| CSV | `{selected['csv']}` |\n| JSON | `{selected['json']}` |\n| reason | {selected['reason']} |\n\n## CSV schema\n\n| item | value |\n|---|---:|\n| rows | {n_rows} |\n| columns | {len(columns)} |\n\nColumns:\n\n```text\n{', '.join(columns)}\n```\n\n## Selected columns\n\n| role | column | status |\n|---|---|---|\n"""
    for row in selected_columns:
        md += f"| {row['role']} | `{row['column']}` | {row['status']} |\n"
    if can_build and inner:
        md += f"""\n## Built diagonal projected inner product\n\n| quantity | value |\n|---|---:|\n| vector length | {inner['n_vector']} |\n| <d,d> | {inner['inner_data_data']:.12g} |\n| <m,m> | {inner['inner_model_model']:.12g} |\n| <d,m> | {inner['inner_data_model']:.12g} |\n| chi2 = <d-m,d-m> | {inner['chi2_model_vs_data_diag']:.12g} |\n| corr(d,m) | {inner['corr_data_model_diag']:.12g} |\n\nThe diagonal inner product was built from the detected error column.\nThis is sufficient to define the first explicit projected vector space, but not sufficient to construct P_m.\n"""
    else:
        md += """\n## Inner product status\n\nThe script selected a canonical file, but did not build the inner product because the data/model/error columns were not detected unambiguously.\n\nUse the output:\n\n    paperVI_EG_T6_numeric_column_profile.csv\n    paperVI_EG_T6_column_role_candidates.csv\n\nto provide a manual column map in T6b.\n"""
    md += """\n## Decisions\n\n| item | decision | basis |\n|---|---|---|\n"""
    for row in decisions:
        md += f"| {row['item']} | {row['decision']} | {row['basis']} |\n"
    md += """\n## Locked conclusion\n\nT6 selects the canonical projected E_G candidate and either builds the first diagonal projected inner product or requests a manual column map.\n\nEven if the inner product is built, T6 does not define P_m.  \nThe matter/private projectors are deferred until the canonical Weyl vector W and matter-correlated direction are explicitly selected.\n\n## Next gate\n\nIf status is BUILT:\n\n    P6-T7 = canonical W vector and projector preparation\n\nIf status is NEEDS_COLUMN_MAP:\n\n    P6-T6b = manual column map for d, m, sigma, and W\n"""
    (OUTDIR / "paperVI_EG_T6_canonical_space.md").write_text(md, encoding="utf-8")
    if vector_rows:
        labels = [r["label"] for r in vector_rows]
        d = [r["data"] for r in vector_rows]
        m = [r["model"] for r in vector_rows]
        sig = [r["sigma"] for r in vector_rows]
        x = list(range(len(labels)))
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.errorbar(x, d, yerr=sig, fmt="o", label="data")
        ax.plot(x, m, marker="s", label="model")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("projected E_G value")
        ax.set_title("Paper VI T6: canonical projected E_G vector")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUTDIR / "fig_P6_T6_canonical_projected_vector.png", dpi=180)
        fig.savefig(OUTDIR / "fig_P6_T6_canonical_projected_vector.pdf")
        plt.close(fig)
    print(status)
    print("Selected:", selected["name"])
    print("CSV:", selected["csv"])
    print("Columns:", columns)
    print("Selected columns:", selected_columns)
    if inner:
        print("chi2_diag:", inner["chi2_model_vs_data_diag"])
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
