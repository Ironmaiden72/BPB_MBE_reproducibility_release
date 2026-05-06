#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper II / D5-A7 — raw SPARC rotation-curve validation of the A6 kernel.

Purpose
-------
A6/A6b validated a compressed, source-only, row-dependent kernel factor using
per-galaxy chi2 anchors.  A7 injects the same kernel into the point-level SPARC
rotation-curve calculation and recomputes chi2 directly from radius/velocity
points.

This is deliberately not a new fit.  It is a raw validation of fixed source-locked
responses and fixed A6 LOO scales.

Outputs
-------
outputs/paperII_D5/paperII_D5A7_raw_rc_validation_of_A6_kernel.md
outputs/paperII_D5/paperII_D5A7_raw_rc_validation_of_A6_kernel.json
outputs/paperII_D5/paperII_D5A7_candidate_summary.csv
outputs/paperII_D5/paperII_D5A7_per_galaxy_raw_rc.csv
outputs/paperII_D5/paperII_D5A7_decisions.csv
outputs/paperII_D5/paperII_D5A7_paperII_wording.tex
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path.cwd()
OUT_DIR = ROOT / "outputs" / "paperII_D5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

A5_PATH = OUT_DIR / "paperII_D5A5_per_galaxy_forensic.csv"
A6_SUMMARY_PATH = OUT_DIR / "paperII_D5A6_model_summary.csv"
A6_PRED_PATH = OUT_DIR / "paperII_D5A6_per_galaxy_predictions.csv"
A6B_DECISIONS_PATH = OUT_DIR / "paperII_D5A6b_decisions.csv"

OUT_MD = OUT_DIR / "paperII_D5A7_raw_rc_validation_of_A6_kernel.md"
OUT_JSON = OUT_DIR / "paperII_D5A7_raw_rc_validation_of_A6_kernel.json"
OUT_SUMMARY = OUT_DIR / "paperII_D5A7_candidate_summary.csv"
OUT_PER_GAL = OUT_DIR / "paperII_D5A7_per_galaxy_raw_rc.csv"
OUT_DECISIONS = OUT_DIR / "paperII_D5A7_decisions.csv"
OUT_TEX = OUT_DIR / "paperII_D5A7_paperII_wording.tex"

PRIMARY_A6_CANDIDATE = "physical_source_fit__chi2_opt_scale__linear__logM200_plus_c200__loo"

FORBIDDEN_A6_PREDICTORS = [
    "active_deltaR1__dR1",
    "amp_proxy_predictions__Areq",
    "canonical__dchi2",
    "chi2_archive",
    "dchi2_archive",
    "dchi2_clean",
    "dchi2_oracle",
    "gap_vs_archive__r_proj_sqrt_a_proj",
    "response_archive",
    "response_oracle",
    "scale_chi2_opt_over_clean",
    "scale_response_archive_over_clean",
]


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv_required(path: Path, **kwargs: Any) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path, **kwargs)


def find_first_existing(candidates: Iterable[Path], label: str) -> Path:
    tried = []
    for p in candidates:
        tried.append(str(p))
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not locate {label}. Tried:\n" + "\n".join(tried))


def locate_raw_inputs() -> Dict[str, Path]:
    chains_candidates = [
        ROOT / "results" / "chains",
        ROOT / "outputs" / "chains",
        ROOT / "outputs" / "paperII_D5",
        ROOT,
    ]
    mrt = find_first_existing([d / "MassModels_Lelli2016c.mrt" for d in chains_candidates], "SPARC MRT")
    ref = find_first_existing([d / "sparc_reference_halo_table_nfw_lcdm_clean.csv" for d in chains_candidates], "NFW reference halo table")
    active = find_first_existing([d / "sparc_active_sample.csv" for d in chains_candidates], "SPARC active sample")
    return {"mrt": mrt, "ref": ref, "active": active}


def clean_name_series(s: pd.Series) -> pd.Series:
    return s.astype("string").fillna("").str.strip()


def infer_prediction_names(pred: pd.DataFrame, ordered_names: List[str]) -> pd.DataFrame:
    """Recover missing A6 prediction names by candidate-local row order.

    A6 originally emitted rows in candidate blocks of 38.  Some intermediate
    versions had blank/NaN names; A6b patched this in its audit output.  A7 does
    the same locally so it remains robust even if the prediction CSV itself was
    not rewritten.
    """
    p = pred.copy()
    if "candidate" not in p.columns:
        raise RuntimeError("A6 prediction table is missing the 'candidate' column")
    if "name" not in p.columns:
        p["name"] = ""
    p["name"] = clean_name_series(p["name"]).astype(object)
    p["_row_in_candidate"] = p.groupby("candidate", sort=False).cumcount()
    by_order = {i: str(n) for i, n in enumerate(ordered_names)}
    blank = clean_name_series(p["name"]).eq("") | clean_name_series(p["name"]).str.lower().isin(["nan", "none", "<na>"])
    if blank.any():
        # Force object dtype before assigning strings into columns that may have
        # been parsed as float64 because all names were empty.
        p["name"] = p["name"].astype(object)
        p.loc[blank, "name"] = p.loc[blank, "_row_in_candidate"].map(by_order).astype(object)
    p["name"] = clean_name_series(p["name"]).astype(str)
    return p


def nfw_v(r: float, v200: float, c200: float, rs: float) -> float:
    def f(u: float) -> float:
        return math.log(1.0 + u) - u / (1.0 + u)

    if not all(math.isfinite(x) and x > 0 for x in [r, v200, c200, rs]):
        return math.nan
    x = r / (c200 * rs)
    u = r / rs
    fc = f(c200)
    fu = f(u)
    if x <= 0 or fc <= 0 or fu < 0:
        return math.nan
    return math.sqrt(max(v200 * v200 * fu / (x * fc), 0.0))


def parse_mrt(mrt_path: Path, active_names: Iterable[str]) -> Dict[str, List[Dict[str, float]]]:
    wanted = set(active_names)
    pts: Dict[str, List[Dict[str, float]]] = {}
    with mrt_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) < 59 or not line[:1].isalpha():
                continue
            try:
                name = line[0:11].strip()
                if name not in wanted:
                    continue
                pts.setdefault(name, []).append(
                    {
                        "R": float(line[19:25]),
                        "Vobs": float(line[26:32]),
                        "eVobs": float(line[33:38]),
                        "Vgas": float(line[39:45]),
                        "Vdisk": float(line[46:52]),
                        "Vbul": float(line[53:59]),
                    }
                )
            except Exception:
                continue
    return pts


def finite_or_default(x: Any, default: float = 0.0) -> float:
    try:
        v = float(x)
    except Exception:
        return default
    return v if math.isfinite(v) else default


def chi2_at(
    name: str,
    alpha: float,
    response: float,
    ref_row: pd.Series,
    pts: Dict[str, List[Dict[str, float]]],
) -> float:
    v200 = finite_or_default(ref_row.get("V200"), math.nan)
    c200 = finite_or_default(ref_row.get("C200"), math.nan)
    rs = finite_or_default(ref_row.get("rs"), math.nan)
    ydisk = finite_or_default(ref_row.get("Ydisk"), 0.0)
    ybul = finite_or_default(ref_row.get("Ybul"), 0.0)

    if not all(math.isfinite(x) and x > 0 for x in [v200, c200, rs]):
        return 1e10

    ratio = 1.0 + float(alpha) * float(response)
    if not math.isfinite(ratio) or ratio <= 0:
        return 1e10

    v200_b = v200 * ratio
    rs_b = rs * ratio  # same old Paper-II machinery: c200 fixed, rs co-scales with V200

    c2 = 0.0
    for p in pts.get(name, []):
        ev = p["eVobs"]
        if not math.isfinite(ev) or ev <= 0:
            continue
        vh = nfw_v(p["R"], v200_b, c200, rs_b)
        if not math.isfinite(vh):
            return 1e10
        vbar2 = p["Vgas"] ** 2 + ydisk * p["Vdisk"] ** 2 + ybul * p["Vbul"] ** 2
        vpred = math.sqrt(max(vbar2 + vh * vh, 0.0))
        c2 += ((vpred - p["Vobs"]) / ev) ** 2
    return float(c2)


def build_a6_candidate_frame(
    pred: pd.DataFrame,
    candidate: str,
    a5_base: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    rows = pred.loc[pred["candidate"].astype(str).eq(candidate)].copy()
    if rows.empty:
        raise RuntimeError(f"Could not find A6 prediction candidate: {candidate}")
    needed = {"name", "scale_eval"}
    missing = sorted(needed - set(rows.columns))
    if missing:
        raise RuntimeError(f"A6 prediction candidate {candidate} missing columns: {missing}")

    keep_cols = ["name", "scale_eval"]
    for extra in ["target_kind", "feature_label", "model_type", "features", "quad_status"]:
        if extra in rows.columns:
            keep_cols.append(extra)
    rows = rows[keep_cols].copy()
    rows["scale_eval"] = pd.to_numeric(rows["scale_eval"], errors="coerce")

    base = a5_base.copy()
    out = base.merge(rows[["name", "scale_eval"]], on="name", how="left", validate="one_to_one")
    if out["scale_eval"].isna().any():
        missing_names = out.loc[out["scale_eval"].isna(), "name"].tolist()
        raise RuntimeError(f"Missing A6 scale for {candidate}: {missing_names}")
    out["candidate"] = label
    out["source_candidate"] = candidate
    out["alpha_eval"] = out["alpha_clean"]
    out["response_eval"] = out["scale_eval"] * out["response_clean_kw2"]
    return out


def find_best_c200_only(summary: pd.DataFrame) -> Optional[str]:
    if summary.empty or "candidate" not in summary.columns:
        return None
    s = summary.copy()
    for col in ["total_dchi2"]:
        if col in s.columns:
            s[col] = pd.to_numeric(s[col], errors="coerce")
    mask = (
        s.get("family", pd.Series("", index=s.index)).astype(str).eq("physical_source_fit")
        & s.get("mode", pd.Series("", index=s.index)).astype(str).eq("loo")
        & s.get("feature_label", pd.Series("", index=s.index)).astype(str).eq("c200_ratio")
        & s.get("target_kind", pd.Series("", index=s.index)).astype(str).eq("chi2_opt_scale")
    )
    cand = s.loc[mask].sort_values("total_dchi2", ascending=True) if "total_dchi2" in s.columns else s.loc[mask]
    if cand.empty:
        return None
    return str(cand.iloc[0]["candidate"])


def make_control_frame(a5_base: pd.DataFrame, candidate: str, response_col: str, alpha_col: str, scale: float = math.nan) -> pd.DataFrame:
    out = a5_base.copy()
    out["candidate"] = candidate
    out["source_candidate"] = candidate
    out["scale_eval"] = scale
    out["alpha_eval"] = pd.to_numeric(out[alpha_col], errors="coerce")
    out["response_eval"] = pd.to_numeric(out[response_col], errors="coerce")
    return out


def evaluate_candidate(
    frame: pd.DataFrame,
    ref: pd.DataFrame,
    pts: Dict[str, List[Dict[str, float]]],
    source_totals: Dict[str, float],
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    rows: List[Dict[str, Any]] = []
    for _, r in frame.iterrows():
        name = str(r["name"])
        if name not in ref.index:
            raise RuntimeError(f"{name} missing from reference halo table")
        if name not in pts or not pts[name]:
            raise RuntimeError(f"{name} missing from MRT point-level table")

        alpha = float(r["alpha_eval"])
        response = float(r["response_eval"])
        chi2_ref_raw = chi2_at(name, 0.0, 0.0, ref.loc[name], pts)
        chi2_eval_raw = chi2_at(name, alpha, response, ref.loc[name], pts)
        dchi2_raw = chi2_eval_raw - chi2_ref_raw

        rows.append(
            {
                "candidate": str(r["candidate"]),
                "source_candidate": str(r.get("source_candidate", r["candidate"])),
                "name": name,
                "alpha_eval": alpha,
                "response_eval": response,
                "scale_eval": float(r.get("scale_eval", math.nan)) if pd.notna(r.get("scale_eval", math.nan)) else math.nan,
                "response_clean_kw2": float(r.get("response_clean_kw2", math.nan)),
                "response_archive": float(r.get("response_archive", math.nan)),
                "response_oracle": float(r.get("response_oracle", math.nan)),
                "alpha_clean": float(r.get("alpha_clean", math.nan)),
                "alpha_oracle": float(r.get("alpha_oracle", math.nan)),
                "chi2_ref_raw": chi2_ref_raw,
                "chi2_eval_raw": chi2_eval_raw,
                "dchi2_raw": dchi2_raw,
                "dchi2_clean_anchor": float(r.get("dchi2_clean", math.nan)),
                "dchi2_archive_anchor": float(r.get("dchi2_archive", math.nan)),
                "dchi2_oracle_anchor": float(r.get("dchi2_oracle", math.nan)),
                "n_points": len(pts[name]),
                "ratio_eval": 1.0 + alpha * response,
                "negative_chi2_flag": int(chi2_eval_raw < 0),
                "invalid_ratio_flag": int((1.0 + alpha * response) <= 0),
            }
        )

    per = pd.DataFrame(rows)
    total_chi2 = float(per["chi2_eval_raw"].sum())
    total_dchi2 = float(per["dchi2_raw"].sum())
    clean_total = source_totals["clean"]
    archive_total = source_totals["archive"]
    oracle_total = source_totals["oracle"]
    gap = clean_total - archive_total
    recovered = (clean_total - total_dchi2) / gap if gap != 0 else math.nan
    summary = {
        "candidate": str(frame["candidate"].iloc[0]),
        "source_candidate": str(frame["source_candidate"].iloc[0]),
        "n": int(len(per)),
        "total_chi2_raw": total_chi2,
        "total_dchi2_raw": total_dchi2,
        "improvement_positive_raw": -total_dchi2,
        "delta_dchi2_vs_clean_raw": total_dchi2 - clean_total,
        "delta_dchi2_vs_archive_raw": total_dchi2 - archive_total,
        "abs_delta_dchi2_vs_archive_raw": abs(total_dchi2 - archive_total),
        "recovered_gap_fraction_raw": recovered,
        "delta_dchi2_vs_oracle_raw": total_dchi2 - oracle_total,
        "negative_chi2_rows": int(per["negative_chi2_flag"].sum()),
        "invalid_ratio_rows": int(per["invalid_ratio_flag"].sum()),
        "scale_mean": float(per["scale_eval"].mean(skipna=True)) if "scale_eval" in per else math.nan,
        "scale_median": float(per["scale_eval"].median(skipna=True)) if "scale_eval" in per else math.nan,
        "scale_min": float(per["scale_eval"].min(skipna=True)) if "scale_eval" in per else math.nan,
        "scale_max": float(per["scale_eval"].max(skipna=True)) if "scale_eval" in per else math.nan,
    }
    return summary, per


def md_table(df: pd.DataFrame, cols: List[str], n: Optional[int] = None) -> str:
    d = df[cols].copy()
    if n is not None:
        d = d.head(n)
    return d.to_markdown(index=False)


def main() -> None:
    for p in [A5_PATH, A6_SUMMARY_PATH, A6_PRED_PATH]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required A7 input: {p}")

    raw_paths = locate_raw_inputs()
    a5 = read_csv_required(A5_PATH)
    a6_summary = read_csv_required(A6_SUMMARY_PATH)
    a6_pred_raw = read_csv_required(A6_PRED_PATH)
    ref = read_csv_required(raw_paths["ref"])

    required_a5 = ["name", "alpha_clean", "response_clean_kw2", "response_archive", "dchi2_clean", "dchi2_archive"]
    missing = [c for c in required_a5 if c not in a5.columns]
    if missing:
        raise RuntimeError(f"Missing required A5 columns: {missing}")
    if "response_oracle" not in a5.columns:
        a5["response_oracle"] = a5["response_archive"]
    if "dchi2_oracle" not in a5.columns:
        raise RuntimeError("A5 must contain dchi2_oracle for oracle comparison")
    if "amp_proxy_predictions__Areq" not in a5.columns:
        raise RuntimeError("A5 must contain amp_proxy_predictions__Areq to build raw oracle alpha")

    a5 = a5.copy()
    a5["name"] = clean_name_series(a5["name"]).astype(str)
    a5 = a5.loc[a5["name"].ne("")].copy()
    ordered_names = a5["name"].tolist()
    if len(ordered_names) != 38 or len(set(ordered_names)) != 38:
        raise RuntimeError(f"A7 expects 38 unique A5 galaxies; got n={len(ordered_names)}, unique={len(set(ordered_names))}")

    for c in [
        "alpha_clean",
        "response_clean_kw2",
        "response_archive",
        "response_oracle",
        "dchi2_clean",
        "dchi2_archive",
        "dchi2_oracle",
        "amp_proxy_predictions__Areq",
    ]:
        a5[c] = pd.to_numeric(a5[c], errors="coerce")
    a5["alpha_oracle"] = -a5["amp_proxy_predictions__Areq"]

    ref["name"] = clean_name_series(ref["name"]).astype(str)
    ref = ref.set_index("name", drop=False)

    a6_pred = infer_prediction_names(a6_pred_raw, ordered_names)

    # Source totals are used only for comparison metrics.  They are the A5/D5
    # source-locked anchors, while raw totals are recomputed independently below.
    source_totals = {
        "clean": float(a5["dchi2_clean"].sum()),
        "archive": float(a5["dchi2_archive"].sum()),
        "oracle": float(a5["dchi2_oracle"].sum()),
    }

    pts = parse_mrt(raw_paths["mrt"], ordered_names)
    missing_pts = sorted(set(ordered_names) - set(pts))
    if missing_pts:
        raise RuntimeError(f"Missing MRT point-level rows for: {missing_pts}")

    # Build candidates.
    frames: List[pd.DataFrame] = []
    frames.append(make_control_frame(a5, "clean_unscaled_kW2_raw", "response_clean_kw2", "alpha_clean", 1.0))
    frames.append(make_control_frame(a5, "archive_response_lock_raw", "response_archive", "alpha_clean", math.nan))
    frames.append(make_control_frame(a5, "oracle_archive_response_raw", "response_oracle", "alpha_oracle", math.nan))
    frames.append(build_a6_candidate_frame(a6_pred, PRIMARY_A6_CANDIDATE, a5, "A6_primary_logM200_plus_c200_kernel_raw"))

    best_c200 = find_best_c200_only(a6_summary)
    if best_c200:
        frames.append(build_a6_candidate_frame(a6_pred, best_c200, a5, "A6_c200_only_kernel_raw"))

    # Optional: if a robust c200 linear candidate exists, include it too when not
    # identical to the chosen best c200-only candidate.
    linear_c200 = "physical_source_fit__chi2_opt_scale__linear__c200_ratio__loo"
    if linear_c200 != best_c200 and linear_c200 in set(a6_pred["candidate"].astype(str)):
        frames.append(build_a6_candidate_frame(a6_pred, linear_c200, a5, "A6_c200_linear_kernel_raw"))

    summaries: List[Dict[str, Any]] = []
    per_all: List[pd.DataFrame] = []
    for fr in frames:
        summary, per = evaluate_candidate(fr, ref, pts, source_totals)
        summaries.append(summary)
        per_all.append(per)

    summary_df = pd.DataFrame(summaries)
    per_df = pd.concat(per_all, ignore_index=True)

    # Raw control errors against source anchors.  These are expected to be tiny if
    # the same raw machinery/inputs are being used; nonzero values reveal an input
    # or rounding mismatch.
    def total_for(cand: str) -> float:
        return float(summary_df.loc[summary_df["candidate"].eq(cand), "total_dchi2_raw"].iloc[0])

    raw_control_errors = {
        "clean_raw_minus_A5_clean": total_for("clean_unscaled_kW2_raw") - source_totals["clean"],
        "archive_raw_minus_A5_archive": total_for("archive_response_lock_raw") - source_totals["archive"],
        "oracle_raw_minus_A5_oracle": total_for("oracle_archive_response_raw") - source_totals["oracle"],
    }

    primary = summary_df.loc[summary_df["candidate"].eq("A6_primary_logM200_plus_c200_kernel_raw")].iloc[0].to_dict()
    c200_summary = None
    if "A6_c200_only_kernel_raw" in set(summary_df["candidate"]):
        c200_summary = summary_df.loc[summary_df["candidate"].eq("A6_c200_only_kernel_raw")].iloc[0].to_dict()

    # Pass criteria are intentionally explicit and conservative.
    primary_rec = float(primary["recovered_gap_fraction_raw"])
    primary_delta_archive = float(primary["delta_dchi2_vs_archive_raw"])
    control_ok = all(abs(v) < 1e-6 for v in raw_control_errors.values())
    no_invalid = int(primary["invalid_ratio_rows"]) == 0 and int(primary["negative_chi2_rows"]) == 0
    strong_raw = primary_rec >= 0.70 and no_invalid
    c200_ok = bool(c200_summary) and float(c200_summary["recovered_gap_fraction_raw"]) >= 0.65

    if strong_raw and control_ok and c200_ok:
        status = "D5A7_RAW_RC_A6_KERNEL_PASS"
        d5_status = "D5_AMPLITUDE_ROWDEPENDENT_KERNEL_RAW_RC_VALIDATED"
        closure = "PASS_RAW_RC"
    elif strong_raw and control_ok:
        status = "D5A7_RAW_RC_A6_PRIMARY_PASS_C200_CONTROL_WEAK"
        d5_status = "D5_A6_KERNEL_RAW_RC_PRIMARY_VALIDATED_C200_REVIEW"
        closure = "PASS_PRIMARY_RAW_RC_CONTROL_REVIEW"
    elif strong_raw:
        status = "D5A7_RAW_RC_A6_PRIMARY_PASS_CONTROL_MISMATCH_REVIEW"
        d5_status = "D5_A6_KERNEL_RAW_RC_PRIMARY_VALIDATED_INPUT_REVIEW"
        closure = "PASS_PRIMARY_RAW_RC_INPUT_REVIEW"
    else:
        status = "D5A7_RAW_RC_A6_KERNEL_FAIL"
        d5_status = "D5_A6_KERNEL_NOT_VALIDATED_BY_RAW_RC"
        closure = "FAIL_RAW_RC"

    # Sort outputs for readability.
    summary_df = summary_df.sort_values("total_dchi2_raw", ascending=True).reset_index(drop=True)
    per_df = per_df.sort_values(["candidate", "dchi2_raw"], ascending=[True, True]).reset_index(drop=True)

    summary_df.to_csv(OUT_SUMMARY, index=False)
    per_df.to_csv(OUT_PER_GAL, index=False)

    decisions = pd.DataFrame(
        [
            {
                "item": "D5A7_status",
                "decision": status,
                "basis": f"primary recovered_gap_fraction_raw={primary_rec:.6f}; delta_vs_archive={primary_delta_archive:.6f}; control_ok={control_ok}; c200_ok={c200_ok}.",
            },
            {
                "item": "D5_status",
                "decision": d5_status,
                "basis": "A7 recomputes raw point-level SPARC rotation-curve chi2 from MRT using fixed A6 source-only kernel scales.",
            },
            {
                "item": "raw_controls",
                "decision": "PASS" if control_ok else "REVIEW",
                "basis": json.dumps(raw_control_errors, sort_keys=True),
            },
            {
                "item": "primary_A6_kernel_raw",
                "decision": "PASS" if strong_raw else "FAIL",
                "basis": json.dumps(primary, sort_keys=True),
            },
            {
                "item": "c200_only_control_raw",
                "decision": "PASS" if c200_ok else "REVIEW_OR_FAIL",
                "basis": json.dumps(c200_summary, sort_keys=True) if c200_summary else "No c200-only LOO candidate found.",
            },
            {
                "item": "closure_decision",
                "decision": closure,
                "basis": "Raw RC validation threshold: primary recovered gap >=0.70 with no invalid/negative chi2 rows; c200-only control target >=0.65.",
            },
            {
                "item": "next_gate",
                "decision": "D5A8_PAPERII_WORDING_OR_FINAL_D5_LOCK" if strong_raw else "D5A7B_DEBUG_RAW_RC_KERNEL_FAILURE",
                "basis": "If A7 passes, lock paper wording. If it fails, inspect raw residuals and clip rows before any amplitude claim.",
            },
        ]
    )
    decisions.to_csv(OUT_DECISIONS, index=False)

    top_primary = (
        per_df.loc[per_df["candidate"].eq("A6_primary_logM200_plus_c200_kernel_raw")]
        .assign(abs_gap_vs_archive=lambda d: (d["dchi2_raw"] - d["dchi2_archive_anchor"]).abs())
        .sort_values("abs_gap_vs_archive", ascending=False)
        .head(15)
    )

    payload = {
        "status": status,
        "d5_status": d5_status,
        "paths": {k: str(v) for k, v in raw_paths.items()},
        "inputs": {
            "A5": str(A5_PATH),
            "A5_sha256": sha256_path(A5_PATH),
            "A6_summary": str(A6_SUMMARY_PATH),
            "A6_summary_sha256": sha256_path(A6_SUMMARY_PATH),
            "A6_predictions": str(A6_PRED_PATH),
            "A6_predictions_sha256": sha256_path(A6_PRED_PATH),
            "mrt_sha256": sha256_path(raw_paths["mrt"]),
            "ref_sha256": sha256_path(raw_paths["ref"]),
            "active_sha256": sha256_path(raw_paths["active"]),
        },
        "source_totals": source_totals,
        "raw_control_errors": raw_control_errors,
        "primary_summary": primary,
        "c200_summary": c200_summary,
        "best_c200_candidate": best_c200,
        "n_galaxies": len(ordered_names),
        "n_mrt_loaded": len(pts),
        "decisions": decisions.to_dict(orient="records"),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    tex = rf"""
% Paper II / D5-A7 wording block — autogenerated
\paragraph{{D5-A7 raw rotation-curve validation.}}
We injected the A6 source-only row-dependent kernel directly into the point-level SPARC rotation-curve reconstruction and recomputed the NFW+barion $\chi^2$ from the Lelli et al. mass-model table. The primary A6 kernel, trained only on source-side halo/kernel variables $\log M_{{200}}$ and $c_{{200}}/c_{{200}}^{{\Lambda\mathrm{{CDM}}}}$, gives
\[
\Delta\chi^2_{{\rm A6,raw}} = {primary['total_dchi2_raw']:.6f},
\qquad
f_{{\rm gap}} = {primary['recovered_gap_fraction_raw']:.6f},
\]
relative to the clean kW$^2$ and archive-response locks. The raw controls reproduce the clean/archive/oracle anchors with errors
$({raw_control_errors['clean_raw_minus_A5_clean']:.3e}, {raw_control_errors['archive_raw_minus_A5_archive']:.3e}, {raw_control_errors['oracle_raw_minus_A5_oracle']:.3e})$.
The D5-A7 status is \texttt{{{status}}}; final interpretation should follow this gate status and not the compressed A6 result alone.
""".strip() + "\n"
    OUT_TEX.write_text(tex, encoding="utf-8")

    md: List[str] = []
    md.append("# Paper II / D5-A7 — raw RC validation of the A6 kernel\n")
    md.append("## Status\n")
    md.append(f"**{status}**\n\n")
    md.append(f"D5 status: **{d5_status}**\n\n")
    md.append("A7 injects the fixed A6 row-dependent kernel into the point-level SPARC rotation-curve calculation and recomputes raw chi² from the MRT radius/velocity table.\n\n")

    md.append("## Inputs\n\n")
    inp_rows = []
    for label, path in [
        ("A5 forensic", A5_PATH),
        ("A6 model summary", A6_SUMMARY_PATH),
        ("A6 per-galaxy predictions", A6_PRED_PATH),
        ("SPARC MRT", raw_paths["mrt"]),
        ("NFW reference", raw_paths["ref"]),
        ("active sample", raw_paths["active"]),
    ]:
        inp_rows.append({"input": label, "path": str(path), "exists": path.exists(), "sha256": sha256_path(path) if path.exists() else ""})
    md.append(pd.DataFrame(inp_rows).to_markdown(index=False))
    md.append("\n\n")

    md.append("## Source totals and raw control errors\n\n")
    controls_df = pd.DataFrame(
        [
            {"quantity": "source clean Δχ²", "value": source_totals["clean"]},
            {"quantity": "source archive Δχ²", "value": source_totals["archive"]},
            {"quantity": "source oracle Δχ²", "value": source_totals["oracle"]},
            {"quantity": "raw clean minus source clean", "value": raw_control_errors["clean_raw_minus_A5_clean"]},
            {"quantity": "raw archive minus source archive", "value": raw_control_errors["archive_raw_minus_A5_archive"]},
            {"quantity": "raw oracle minus source oracle", "value": raw_control_errors["oracle_raw_minus_A5_oracle"]},
        ]
    )
    md.append(controls_df.to_markdown(index=False))
    md.append("\n\n")

    md.append("## Candidate summary\n\n")
    show_cols = [
        "candidate",
        "total_dchi2_raw",
        "delta_dchi2_vs_archive_raw",
        "recovered_gap_fraction_raw",
        "delta_dchi2_vs_oracle_raw",
        "scale_median",
        "invalid_ratio_rows",
        "negative_chi2_rows",
    ]
    md.append(md_table(summary_df, show_cols))
    md.append("\n\n")

    md.append("## Primary A6 raw residual drivers\n\n")
    top_show_cols = [
        "name",
        "scale_eval",
        "dchi2_raw",
        "dchi2_archive_anchor",
        "abs_gap_vs_archive",
        "chi2_eval_raw",
        "chi2_ref_raw",
        "ratio_eval",
    ]
    md.append(md_table(top_primary, top_show_cols, n=15))
    md.append("\n\n")

    md.append("## Decisions\n\n")
    md.append(decisions.to_markdown(index=False))
    md.append("\n\n")
    md.append("## Interpretation lock\n\n")
    if strong_raw:
        md.append(
            "D5-A7 validates the A6 kernel at raw rotation-curve level for the primary source-only model. "
            "This promotes D5 beyond the compressed chi²-anchor test, subject to the c200-only control and wording gate recorded above.\n"
        )
    else:
        md.append(
            "D5-A7 does not validate the A6 kernel at raw rotation-curve level. The amplitude claim must remain compressed-only until the raw residuals are debugged.\n"
        )

    OUT_MD.write_text("".join(md), encoding="utf-8")

    print(f"Saved: {OUT_MD}")
    print(f"Saved: {OUT_JSON}")
    print(f"Saved: {OUT_SUMMARY}")
    print(f"Saved: {OUT_PER_GAL}")
    print(f"Saved: {OUT_DECISIONS}")
    print(f"Saved: {OUT_TEX}")
    print(f"Status: {status}")
    print(f"D5 status: {d5_status}")
    print(f"Primary recovered gap fraction raw: {primary_rec}")
    print(f"Primary total dchi2 raw: {primary['total_dchi2_raw']}")
    print(f"Primary delta vs archive raw: {primary_delta_archive}")
    if c200_summary:
        print(f"c200-only recovered gap fraction raw: {c200_summary['recovered_gap_fraction_raw']}")
    print(f"control_ok: {control_ok}")


if __name__ == "__main__":
    main()
