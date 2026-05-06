#!/usr/bin/env python3
"""
Paper VI / E_G T6b — candidate column-map audit

Purpose
-------
T6 selected a canonical projected E_G candidate, but correctly refused to build
<A,B>_EG because the selected CSV contains projected theory curves, not an
observed data vector and errors.

T6b audits the top T5/T6 candidate files and prepares a manual column map.
It does not build P_m, does not run any likelihood, and touches no NS files.
"""
from __future__ import annotations

from pathlib import Path
import csv
import json
import math
from typing import Any, Dict, List, Tuple

OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIORITY_FILES = [
    ("PRIMARY_PROJECTION_THEORY_SPACE", "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.csv"),
    ("PRIMARY_PROJECTION_THEORY_SPACE_JSON", "outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.json"),
    ("PLANCK_PR4_BOSS_PROJECTION_THEORY_SPACE", "outputs/cosmology/eg_projected_wenzl_planck_pr4_boss_final.csv"),
    ("PLANCK_PR4_BOSS_PROJECTION_THEORY_SPACE_JSON", "outputs/cosmology/eg_projected_wenzl_planck_pr4_boss_final.json"),
    ("EG_DIAGNOSTIC_BEST_PLANCK_PR4_BOSS", "outputs/cosmology/eg_diagnostic_best_planck_pr4_boss_seed01.csv"),
    ("EG_DIAGNOSTIC_BEST_PLANCK_PR4_BOSS_JSON", "outputs/cosmology/eg_diagnostic_best_planck_pr4_boss_seed01.json"),
    ("AUTO_CROSS_GLOBAL_R", "outputs/cosmology/eg_projected_auto_cross_global_r_seed01.csv"),
    ("AUTO_CROSS_GLOBAL_R_JSON", "outputs/cosmology/eg_projected_auto_cross_global_r_seed01.json"),
    ("FREE_R_SCAN", "outputs/cosmology/autocross_csigma_free_r_scan_final.csv"),
    ("FREE_R_SCAN_JSON", "outputs/cosmology/autocross_csigma_free_r_scan_final.json"),
    ("ACT_EG_AUTOCROSS_JOINT", "outputs/cosmology/act_eg_autocross_joint_summary_seed01.csv"),
    ("ACT_EG_AUTOCROSS_JOINT_JSON", "outputs/cosmology/act_eg_autocross_joint_summary_seed01.json"),
    ("WEYL_PROXY_RAW", "outputs/cmb/cmb_lensing_weyl_proxy_map1782_Planck_PR4_raw.csv"),
]

def parse_float(x: Any):
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "none", "null"}:
        return None
    try:
        v = float(s)
    except Exception:
        return None
    if not math.isfinite(v):
        return None
    return v

def read_csv_file(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(8192)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        rows = list(reader)
        cols = reader.fieldnames or []
    return cols, rows

def flatten_json(obj: Any, prefix: str = "", out=None, max_keys: int = 600):
    if out is None:
        out = []
    if len(out) >= max_keys:
        return out
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.append((key, v if not isinstance(v, (dict, list)) else None))
            flatten_json(v, key, out, max_keys)
            if len(out) >= max_keys:
                break
    elif isinstance(obj, list):
        out.append((prefix + ".__len__", len(obj)))
        if obj:
            flatten_json(obj[0], prefix + "[0]", out, max_keys)
    return out

def numeric_profile(cols: List[str], rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    n = max(1, len(rows))
    prof = []
    for col in cols:
        vals = [parse_float(r.get(col)) for r in rows]
        good = [v for v in vals if v is not None]
        prof.append({
            "column": col,
            "numeric_count": len(good),
            "n_rows": len(rows),
            "numeric_fraction": len(good) / n,
            "min": min(good) if good else "",
            "max": max(good) if good else "",
            "mean": sum(good) / len(good) if good else "",
        })
    return prof

def role_scores(col: str) -> Dict[str, int]:
    c = col.lower()
    def pos(*tokens): return sum(t in c for t in tokens)
    def neg(*tokens): return sum(t in c for t in tokens)
    return {
        "bin_or_axis": 2*pos("sample", "bin", "z", "zeff", "ell", "k_eff", "label", "dataset"),
        "observed_data": 4*pos("obs", "observed", "data", "measured", "measurement") - 3*neg("model", "pred", "theory", "lcdm", "bpb", "sigma", "err", "chi2"),
        "model_prediction": 3*pos("model", "pred", "theory", "projected", "lcdm", "bpb", "map1782", "strict") - 3*neg("obs", "observed", "data", "err", "sigma", "chi2"),
        "uncertainty": 4*pos("err", "error", "sigma", "unc", "std", "cov") - 2*neg("projected", "bpb_", "eg_bpb", "model", "pred", "weyl"),
        "weyl_or_sigma_response": 3*pos("weyl", "sigma", "lens", "kappa", "wdelta", "phi", "auto", "cross", "act_sigma", "eg_sigma") - 2*neg("err", "error", "chi2"),
        "chi2_or_fit_summary": 3*pos("chi2", "best", "delta", "scan", "fit", "r_"),
    }

def classify_file(cols: List[str], rows: List[Dict[str, str]]) -> Dict[str, Any]:
    prof = numeric_profile(cols, rows)
    numeric = {p["column"]: p for p in prof if p["numeric_fraction"] >= 0.7}
    role_candidates: Dict[str, List[Tuple[int, str]]] = {
        "bin_or_axis": [], "observed_data": [], "model_prediction": [],
        "uncertainty": [], "weyl_or_sigma_response": [], "chi2_or_fit_summary": []
    }
    for col in cols:
        scores = role_scores(col)
        for role, score in scores.items():
            if role != "bin_or_axis" and col not in numeric:
                continue
            if score > 0:
                role_candidates[role].append((score, col))
    for role in role_candidates:
        role_candidates[role].sort(key=lambda x: (-x[0], x[1]))

    has_obs = bool(role_candidates["observed_data"])
    has_model = bool(role_candidates["model_prediction"])
    has_unc = bool(role_candidates["uncertainty"])
    has_weyl = bool(role_candidates["weyl_or_sigma_response"])

    if has_obs and has_model and has_unc:
        verdict = "CAN_BUILD_DIAGONAL_INNER_PRODUCT_CANDIDATE"
    elif has_model and has_weyl and not has_obs:
        verdict = "THEORY_PROJECTION_SPACE_NO_OBSERVED_DATA"
    elif has_obs and has_unc and not has_model:
        verdict = "DATA_VECTOR_CANDIDATE_NEEDS_MODEL_LINK"
    elif has_weyl:
        verdict = "WEYL_RESPONSE_SUPPORTING_CANDIDATE"
    else:
        verdict = "SCHEMA_REVIEW_REQUIRED"

    return {
        "verdict": verdict,
        "role_candidates": role_candidates,
        "numeric_profile": prof,
    }

def write_dict_csv(path: Path, rows: List[Dict[str, Any]], fields: List[str]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    file_rows = []
    column_rows = []
    role_rows = []
    json_rows = []
    first_rows_md = []

    for tag, raw_path in PRIORITY_FILES:
        p = Path(raw_path)
        exists = p.exists()
        row = {
            "tag": tag,
            "path": raw_path,
            "exists": exists,
            "suffix": p.suffix,
            "n_rows": "",
            "n_columns_or_keys": "",
            "verdict": "MISSING",
        }

        if not exists:
            file_rows.append(row)
            continue

        if p.suffix == ".csv":
            try:
                cols, rows = read_csv_file(p)
                cls = classify_file(cols, rows)
                row.update({
                    "n_rows": len(rows),
                    "n_columns_or_keys": len(cols),
                    "verdict": cls["verdict"],
                })

                for prof in cls["numeric_profile"]:
                    rr = {"tag": tag, "path": raw_path}
                    rr.update(prof)
                    column_rows.append(rr)

                for role, cands in cls["role_candidates"].items():
                    for score, col in cands[:8]:
                        role_rows.append({
                            "tag": tag,
                            "path": raw_path,
                            "role": role,
                            "column": col,
                            "score": score,
                        })

                first_rows_md.append((tag, raw_path, cols, rows[:5]))
            except Exception as e:
                row["verdict"] = "CSV_READ_ERROR:" + repr(e)

        elif p.suffix == ".json":
            try:
                data = json.loads(p.read_text(encoding="utf-8", errors="replace"))
                flat = flatten_json(data)
                row.update({
                    "n_rows": "",
                    "n_columns_or_keys": len(flat),
                    "verdict": "JSON_KEYS_AUDITED",
                })
                for key, val in flat[:300]:
                    sval = ""
                    if val is not None and not isinstance(val, (dict, list)):
                        sval = str(val)
                        if len(sval) > 160:
                            sval = sval[:160] + "..."
                    json_rows.append({
                        "tag": tag,
                        "path": raw_path,
                        "key": key,
                        "value_sample": sval,
                    })
            except Exception as e:
                row["verdict"] = "JSON_READ_ERROR:" + repr(e)

        file_rows.append(row)

    # Proposed decision from audit.
    buildable = [r for r in file_rows if r["verdict"] == "CAN_BUILD_DIAGONAL_INNER_PRODUCT_CANDIDATE"]
    theory_spaces = [r for r in file_rows if r["verdict"] == "THEORY_PROJECTION_SPACE_NO_OBSERVED_DATA"]
    weyl = [r for r in file_rows if "WEYL" in str(r["verdict"]) or any(rr["tag"] == r["tag"] and rr["role"] == "weyl_or_sigma_response" for rr in role_rows)]

    if buildable:
        status = "P6_T6B_COLUMN_MAP_AUDIT_BUILDABLE_CANDIDATE_FOUND"
    elif theory_spaces and weyl:
        status = "P6_T6B_COLUMN_MAP_AUDIT_THEORY_SPACE_AND_WEYL_FOUND_DATA_LINK_MISSING"
    else:
        status = "P6_T6B_COLUMN_MAP_AUDIT_NEEDS_MANUAL_SOURCE_LINK"

    decisions = [
        {
            "item": "primary_projection_csv",
            "decision": "THEORY_SPACE_NOT_DATA_SPACE",
            "basis": "The primary Wenzl projected CSV has columns sample/ell/k and projected theory curves only; it has no observed data or uncertainty column.",
        },
        {
            "item": "inner_product_build",
            "decision": "DEFER_UNTIL_DATA_ERROR_SOURCE_SELECTED" if not buildable else "BUILDABLE_CANDIDATE_FOUND",
            "basis": "A valid <A,B>_EG requires observed vector and covariance/errors, not just model curves.",
        },
        {
            "item": "W_vector",
            "decision": "USE_PRIMARY_PROJECTION_COLUMNS_AS_W_CANDIDATES_ONLY",
            "basis": "EG_bpb_act_sigma_projected and EG_bpb_eg_sigma_projected are Weyl/Sigma-like projected responses, not observational errors.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T6C_SELECT_DATA_ERROR_SOURCE_OR_USE_MODEL_SPACE_INNER_PRODUCT",
            "basis": "Either connect the Wenzl projected theory space to a separate observed/covariance file, or explicitly define a model-space norm rather than data chi2 inner product.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
        },
        "file_rows": file_rows,
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T6b_column_map_audit.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    write_dict_csv(OUTDIR / "paperVI_EG_T6b_file_verdicts.csv", file_rows, ["tag","path","exists","suffix","n_rows","n_columns_or_keys","verdict"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6b_numeric_columns.csv", column_rows, ["tag","path","column","numeric_count","n_rows","numeric_fraction","min","max","mean"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6b_role_candidates.csv", role_rows, ["tag","path","role","column","score"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6b_json_keys.csv", json_rows, ["tag","path","key","value_sample"])
    write_dict_csv(OUTDIR / "paperVI_EG_T6b_decisions.csv", decisions, ["item","decision","basis"])

    md = f"""# Paper VI / E_G T6b - column map audit

## Status

**{status}**

Non-destructive column-map audit. No likelihood run. NS untouched.

## Why T6 stopped

T6 selected the canonical projected theory-space file:

    outputs/cosmology/eg_projected_wenzl_act_dr6_planck_pr4_boss_final.csv

but correctly refused to build <A,B>_EG because the file contains projected model/theory curves, not an observed data vector plus errors.

## File verdicts

| tag | exists | rows/keys | verdict | path |
|---|---:|---:|---|---|
"""
    for r in file_rows:
        md += f"| {r['tag']} | {r['exists']} | {r['n_rows'] or r['n_columns_or_keys']} | {r['verdict']} | `{r['path']}` |\n"

    md += """
## Important correction

The column:

    EG_bpb_act_sigma_projected

is **not** an observational error column. It is a projected BPB ACT/Sigma model response.  
Therefore it must not be used as sigma_i in:

    <A,B> = sum_i A_i B_i / sigma_i^2

unless a later gate explicitly defines a model-space norm rather than a data chi2 inner product.

## First rows of CSV candidates
"""

    for tag, path, cols, rows in first_rows_md:
        md += f"\n### {tag}\n\n`{path}`\n\nColumns:\n\n```text\n{', '.join(cols)}\n```\n\n"
        if rows:
            md += "| " + " | ".join(cols) + " |\n"
            md += "| " + " | ".join(["---"] * len(cols)) + " |\n"
            for r in rows:
                md += "| " + " | ".join(str(r.get(c, ""))[:80] for c in cols) + " |\n"
        else:
            md += "_No rows._\n"

    md += """
## Decisions

| item | decision | basis |
|---|---|---|
"""
    for d in decisions:
        md += f"| {d['item']} | {d['decision']} | {d['basis']} |\n"

    md += """
## Locked conclusion

T6b prevents a false inner-product construction.  
The primary projected Wenzl CSV is a canonical projected **model space**, not yet a data/covariance space.

Next step: either locate the separate observed E_G vector and errors/covariance that correspond to the Wenzl projected bins, or explicitly define a model-space inner product for Weyl responses only.

## Next gate

P6-T6c should choose one of two routes:

1. Data-space route:
   connect the projected Wenzl bins to observed E_G data and errors/covariance.

2. Model-space route:
   define <A,B>_EG as a normalized window/bin sum over projected Weyl responses,
   independent of observational errors, for the sole purpose of building W/P_m.
"""

    (OUTDIR / "paperVI_EG_T6b_column_map_audit.md").write_text(md, encoding="utf-8")

    print(status)
    print("Wrote outputs to:", OUTDIR)

if __name__ == "__main__":
    main()
