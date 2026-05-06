#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper II / D6 — kernel falsification and final D5/D6 lock.

D6 is a non-destructive audit gate after D5-A8.  It tests whether the
raw-RC-validated A6 row-dependent kernel can be dismissed as a row-by-row
artifact.  It performs:
  1) a source-feature permutation null for the A6 primary LOO kernel;
  2) a top residual-driver removal audit;
  3) a source-side feature ablation audit;
  4) a final wording/source-lock gate for Paper II.

Outputs
-------
outputs/paperII_D5/paperII_D6_kernel_falsification_lock.md
outputs/paperII_D5/paperII_D6_kernel_falsification_lock.json
outputs/paperII_D5/paperII_D6_permutation_null.csv
outputs/paperII_D5/paperII_D6_top_driver_removal.csv
outputs/paperII_D5/paperII_D6_feature_ablation.csv
outputs/paperII_D5/paperII_D6_decisions.csv
outputs/paperII_D5/paperII_D6_claim_source_lock.csv
outputs/paperII_D5/paperII_D6_paperII_wording.txt
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path.cwd()
OUT_DIR = ROOT / "outputs" / "paperII_D5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

A5_PATH = OUT_DIR / "paperII_D5A5_per_galaxy_forensic.csv"
A6_SUMMARY_PATH = OUT_DIR / "paperII_D5A6_model_summary.csv"
A6B_RESAMPLING_PATH = OUT_DIR / "paperII_D5A6b_resampling_summary.csv"
A6B_COEFF_PATH = OUT_DIR / "paperII_D5A6b_coefficients_bootstrap.csv"
A7_SUMMARY_PATH = OUT_DIR / "paperII_D5A7_candidate_summary.csv"
A7_PER_GAL_PATH = OUT_DIR / "paperII_D5A7_per_galaxy_raw_rc.csv"
A8_LOCK_PATH = OUT_DIR / "paperII_D5A8_final_D5_lock.md"
A8_NUMERIC_LOCKS_PATH = OUT_DIR / "paperII_D5A8_numeric_locks.csv"
A8_CLAIM_LOCK_PATH = OUT_DIR / "paperII_D5A8_claim_source_lock.csv"
A7_SCRIPT_PATH = ROOT / "runs" / "paperII" / "write_D5A7_raw_rc_validation_of_A6_kernel.py"

OUT_MD = OUT_DIR / "paperII_D6_kernel_falsification_lock.md"
OUT_JSON = OUT_DIR / "paperII_D6_kernel_falsification_lock.json"
OUT_PERM = OUT_DIR / "paperII_D6_permutation_null.csv"
OUT_TOP = OUT_DIR / "paperII_D6_top_driver_removal.csv"
OUT_ABL = OUT_DIR / "paperII_D6_feature_ablation.csv"
OUT_DECISIONS = OUT_DIR / "paperII_D6_decisions.csv"
OUT_CLAIMS = OUT_DIR / "paperII_D6_claim_source_lock.csv"
OUT_WORDING = OUT_DIR / "paperII_D6_paperII_wording.txt"

PRIMARY_RAW_CANDIDATE = "A6_primary_logM200_plus_c200_kernel_raw"
C200_RAW_CANDIDATE = "A6_c200_only_kernel_raw"
PRIMARY_SOURCE_CANDIDATE = "physical_source_fit__chi2_opt_scale__linear__logM200_plus_c200__loo"

FEATURE_LOGM = "kernel_variants__logM200"
FEATURE_C200 = "kernel_variants__c200_lcdm_ratio"
FEATURE_DKW2 = "kernel_variants__delta_kW2"
FEATURE_DW2 = "kernel_variants__delta_W2"
TARGET_SCALE = "scale_chi2_opt_over_clean"
SCALE_MIN = 1.0e-6
SCALE_MAX = 3.0

N_PERM = int(os.environ.get("D6_N_PERM", "1000"))
SEED = int(os.environ.get("D6_SEED", "260505"))

FORBIDDEN_HINTS = [
    "archive",
    "oracle",
    "dchi2",
    "chi2",
    "response_archive",
    "response_oracle",
    "scale_response_archive_over_clean",
    "scale_chi2_opt_over_clean",
    "gap_vs_archive",
    "Areq",
    "active_deltaR1",
]


def sha256_path(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv_required(path: Path, **kwargs: Any) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path, **kwargs)


def clean_name_series(s: pd.Series) -> pd.Series:
    return s.astype("string").fillna("").str.strip()


def md_table(df: pd.DataFrame, cols: Optional[List[str]] = None, n: Optional[int] = None) -> str:
    d = df.copy()
    if cols is not None:
        d = d[cols]
    if n is not None:
        d = d.head(n)
    if d.empty:
        return "_empty_"
    return d.to_markdown(index=False)


def load_a7_module() -> Any:
    if not A7_SCRIPT_PATH.exists():
        raise FileNotFoundError(
            f"Missing A7 script needed for raw RC evaluator: {A7_SCRIPT_PATH}\n"
            "Run/keep D5-A7 before D6."
        )
    spec = importlib.util.spec_from_file_location("paperII_D5A7_raw_eval", A7_SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import A7 script: {A7_SCRIPT_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def to_numeric_inplace(df: pd.DataFrame, cols: Iterable[str]) -> None:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")


def prepare_inputs(a7: Any) -> Dict[str, Any]:
    required = [A5_PATH, A6_SUMMARY_PATH, A7_SUMMARY_PATH, A7_PER_GAL_PATH]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("D6 requires D5-A5/A6/A7 outputs. Missing:\n" + "\n".join(missing))

    a5 = read_csv_required(A5_PATH).copy()
    a6_summary = read_csv_required(A6_SUMMARY_PATH).copy()
    a7_summary = read_csv_required(A7_SUMMARY_PATH).copy()
    a7_per = read_csv_required(A7_PER_GAL_PATH).copy()
    a6b_resampling = read_csv_required(A6B_RESAMPLING_PATH).copy() if A6B_RESAMPLING_PATH.exists() else pd.DataFrame()
    a6b_coeff = read_csv_required(A6B_COEFF_PATH).copy() if A6B_COEFF_PATH.exists() else pd.DataFrame()
    a8_numeric = read_csv_required(A8_NUMERIC_LOCKS_PATH).copy() if A8_NUMERIC_LOCKS_PATH.exists() else pd.DataFrame()

    if "name" not in a5.columns:
        raise RuntimeError("A5 table is missing name column")
    a5["name"] = clean_name_series(a5["name"]).astype(str)
    a5 = a5.loc[a5["name"].ne("")].copy()
    ordered_names = a5["name"].tolist()
    if len(ordered_names) != 38 or len(set(ordered_names)) != 38:
        raise RuntimeError(f"D6 expects 38 unique galaxies from A5; got n={len(ordered_names)}, unique={len(set(ordered_names))}")

    required_cols = [
        "alpha_clean",
        "response_clean_kw2",
        "response_archive",
        "dchi2_clean",
        "dchi2_archive",
        "dchi2_oracle",
        TARGET_SCALE,
        FEATURE_LOGM,
        FEATURE_C200,
    ]
    miss = [c for c in required_cols if c not in a5.columns]
    if miss:
        raise RuntimeError(f"A5 is missing required D6 columns: {miss}")
    if "response_oracle" not in a5.columns:
        a5["response_oracle"] = a5["response_archive"]
    if "amp_proxy_predictions__Areq" in a5.columns:
        a5["alpha_oracle"] = -pd.to_numeric(a5["amp_proxy_predictions__Areq"], errors="coerce")
    else:
        a5["alpha_oracle"] = a5["alpha_clean"]

    numeric_cols = list(set(required_cols + ["response_oracle", "alpha_oracle", "amp_proxy_predictions__Areq", FEATURE_DKW2, FEATURE_DW2]))
    to_numeric_inplace(a5, numeric_cols)

    # Raw SPARC inputs are located/reused by the A7 module.
    raw_paths = a7.locate_raw_inputs()
    ref = a7.read_csv_required(raw_paths["ref"])
    ref["name"] = a7.clean_name_series(ref["name"]).astype(str)
    ref = ref.set_index("name", drop=False)
    pts = a7.parse_mrt(raw_paths["mrt"], ordered_names)
    missing_pts = sorted(set(ordered_names) - set(pts))
    if missing_pts:
        raise RuntimeError(f"Missing MRT point-level rows for: {missing_pts}")

    source_totals = {
        "clean": float(a5["dchi2_clean"].sum()),
        "archive": float(a5["dchi2_archive"].sum()),
        "oracle": float(a5["dchi2_oracle"].sum()),
    }

    return {
        "a5": a5,
        "a6_summary": a6_summary,
        "a6b_resampling": a6b_resampling,
        "a6b_coeff": a6b_coeff,
        "a7_summary": a7_summary,
        "a7_per": a7_per,
        "a8_numeric": a8_numeric,
        "ordered_names": ordered_names,
        "raw_paths": raw_paths,
        "ref": ref,
        "pts": pts,
        "source_totals": source_totals,
    }


def fit_predict_loo(
    df: pd.DataFrame,
    feature_cols: Sequence[str],
    target_col: str = TARGET_SCALE,
    model_type: str = "linear",
    X_override: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Return LOO predicted positive scales for every row.

    Training rows are the rows with finite positive target and quad_status == ok.
    Rows not usable as training rows are still predicted from the fold model; if a
    fold is underdetermined, the finite median target is used.
    """
    missing = [c for c in feature_cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing feature columns for LOO: {missing}")
    if target_col not in df.columns:
        raise RuntimeError(f"Missing target column for LOO: {target_col}")

    X = X_override if X_override is not None else df[list(feature_cols)].to_numpy(dtype=float)
    y = pd.to_numeric(df[target_col], errors="coerce").to_numpy(dtype=float)
    n = len(df)
    if X.shape != (n, len(feature_cols)):
        raise RuntimeError(f"Invalid X_override shape: got {X.shape}, expected {(n, len(feature_cols))}")

    finite_y = np.isfinite(y) & (y > 0)
    if "quad_status" in df.columns:
        fit_ok = clean_name_series(df["quad_status"]).str.lower().eq("ok").to_numpy() & finite_y
    else:
        fit_ok = finite_y
    if fit_ok.sum() < len(feature_cols) + 2:
        raise RuntimeError(f"Not enough fit rows for {feature_cols}: {fit_ok.sum()}")

    fallback = float(np.nanmedian(y[fit_ok]))
    preds = np.full(n, fallback, dtype=float)

    for i in range(n):
        train = fit_ok.copy()
        train[i] = False
        if train.sum() < len(feature_cols) + 2:
            preds[i] = fallback
            continue

        Xtr = X[train, :].astype(float)
        ytr = y[train].astype(float)
        xte = X[i : i + 1, :].astype(float)

        # Median-impute feature NaNs using training medians.
        med = np.nanmedian(Xtr, axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        Xtr = np.where(np.isfinite(Xtr), Xtr, med)
        xte = np.where(np.isfinite(xte), xte, med)

        mu = np.mean(Xtr, axis=0)
        sd = np.std(Xtr, axis=0)
        sd = np.where((np.isfinite(sd)) & (sd > 0), sd, 1.0)
        Xtr_s = (Xtr - mu) / sd
        xte_s = (xte - mu) / sd
        A = np.column_stack([np.ones(Xtr_s.shape[0]), Xtr_s])
        b = np.array([1.0] + list(xte_s[0]), dtype=float)

        if model_type == "loglinear":
            ok = ytr > 0
            if ok.sum() < len(feature_cols) + 2:
                preds[i] = fallback
                continue
            A2 = A[ok]
            yy = np.log(ytr[ok])
            beta, *_ = np.linalg.lstsq(A2, yy, rcond=None)
            pred = math.exp(float(b @ beta))
        elif model_type == "linear":
            beta, *_ = np.linalg.lstsq(A, ytr, rcond=None)
            pred = float(b @ beta)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")

        if not math.isfinite(pred):
            pred = fallback
        preds[i] = min(max(pred, SCALE_MIN), SCALE_MAX)
    return preds


def make_frame_from_scales(a5: pd.DataFrame, scales: Sequence[float], candidate: str, source_candidate: str) -> pd.DataFrame:
    if len(scales) != len(a5):
        raise RuntimeError(f"Scale length mismatch: {len(scales)} vs {len(a5)}")
    out = a5.copy()
    out["candidate"] = candidate
    out["source_candidate"] = source_candidate
    out["scale_eval"] = np.asarray(scales, dtype=float)
    out["alpha_eval"] = out["alpha_clean"]
    out["response_eval"] = out["scale_eval"] * out["response_clean_kw2"]
    return out


def evaluate_scales_raw(
    a7: Any,
    a5: pd.DataFrame,
    scales: Sequence[float],
    candidate: str,
    source_candidate: str,
    ref: pd.DataFrame,
    pts: Dict[str, Any],
    source_totals: Dict[str, float],
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    frame = make_frame_from_scales(a5, scales, candidate, source_candidate)
    summary, per = a7.evaluate_candidate(frame, ref, pts, source_totals)
    return summary, per


def recovered_fraction(total_dchi2: float, clean_total: float, archive_total: float) -> float:
    gap = clean_total - archive_total
    return (clean_total - total_dchi2) / gap if gap != 0 else math.nan


def run_permutation_null(
    a7: Any,
    a5: pd.DataFrame,
    ref: pd.DataFrame,
    pts: Dict[str, Any],
    source_totals: Dict[str, float],
    n_perm: int,
    seed: int,
) -> Tuple[pd.DataFrame, Dict[str, Any], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    features = [FEATURE_LOGM, FEATURE_C200]
    X = a5[features].to_numpy(dtype=float)

    real_scales = fit_predict_loo(a5, features, TARGET_SCALE, "linear", X_override=X)
    real_summary, real_per = evaluate_scales_raw(
        a7, a5, real_scales, "D6_recomputed_primary_LOO_raw", PRIMARY_SOURCE_CANDIDATE, ref, pts, source_totals
    )
    real_total = float(real_summary["total_dchi2_raw"])
    real_rec = float(real_summary["recovered_gap_fraction_raw"])

    rows: List[Dict[str, Any]] = []
    for j in range(n_perm):
        perm = rng.permutation(len(a5))
        Xp = X[perm, :]
        scales = fit_predict_loo(a5, features, TARGET_SCALE, "linear", X_override=Xp)
        summary, _ = evaluate_scales_raw(
            a7,
            a5,
            scales,
            f"D6_perm_{j:05d}",
            "permuted_source_features_logM200_plus_c200",
            ref,
            pts,
            source_totals,
        )
        rows.append(
            {
                "perm_id": j,
                "total_dchi2_raw": float(summary["total_dchi2_raw"]),
                "recovered_gap_fraction_raw": float(summary["recovered_gap_fraction_raw"]),
                "delta_dchi2_vs_archive_raw": float(summary["delta_dchi2_vs_archive_raw"]),
                "scale_mean": float(summary.get("scale_mean", math.nan)),
                "scale_median": float(summary.get("scale_median", math.nan)),
                "scale_min": float(summary.get("scale_min", math.nan)),
                "scale_max": float(summary.get("scale_max", math.nan)),
                "negative_chi2_rows": int(summary.get("negative_chi2_rows", 0)),
                "invalid_ratio_rows": int(summary.get("invalid_ratio_rows", 0)),
            }
        )
    null = pd.DataFrame(rows)
    rec = null["recovered_gap_fraction_raw"].to_numpy(dtype=float)
    dchi = null["total_dchi2_raw"].to_numpy(dtype=float)
    # Higher recovered fraction is better; lower total_dchi2 is better.
    p_rec = (1.0 + float(np.sum(rec >= real_rec))) / (len(rec) + 1.0)
    p_dchi = (1.0 + float(np.sum(dchi <= real_total))) / (len(dchi) + 1.0)
    null_mean = float(np.mean(rec))
    null_std = float(np.std(rec, ddof=1)) if len(rec) > 1 else math.nan
    z = (real_rec - null_mean) / null_std if null_std and math.isfinite(null_std) and null_std > 0 else math.nan
    q = np.quantile(rec, [0.50, 0.90, 0.95, 0.99])
    audit = {
        "n_perm": int(n_perm),
        "seed": int(seed),
        "real_total_dchi2_raw_recomputed": real_total,
        "real_recovered_gap_fraction_raw_recomputed": real_rec,
        "null_recovered_mean": null_mean,
        "null_recovered_std": null_std,
        "null_recovered_p50": float(q[0]),
        "null_recovered_p90": float(q[1]),
        "null_recovered_p95": float(q[2]),
        "null_recovered_p99": float(q[3]),
        "p_perm_recovered_ge_real": p_rec,
        "p_perm_total_dchi2_le_real": p_dchi,
        "z_recovered": z,
    }
    return null, audit, real_per


def top_driver_removal(a7_per: pd.DataFrame) -> pd.DataFrame:
    p = a7_per.loc[a7_per["candidate"].astype(str).eq(PRIMARY_RAW_CANDIDATE)].copy()
    if p.empty:
        raise RuntimeError(f"A7 per-galaxy raw table is missing {PRIMARY_RAW_CANDIDATE}")
    for c in ["dchi2_raw", "dchi2_clean_anchor", "dchi2_archive_anchor", "dchi2_oracle_anchor", "chi2_eval_raw"]:
        if c in p.columns:
            p[c] = pd.to_numeric(p[c], errors="coerce")
    p["abs_gap_vs_archive"] = (p["dchi2_raw"] - p["dchi2_archive_anchor"]).abs()
    ordered = p.sort_values("abs_gap_vs_archive", ascending=False).reset_index(drop=True)
    top_names = ordered["name"].astype(str).tolist()

    rows: List[Dict[str, Any]] = []
    for k in [0, 1, 3, 5, 10]:
        removed = set(top_names[:k])
        sub = p.loc[~p["name"].astype(str).isin(removed)].copy()
        clean = float(sub["dchi2_clean_anchor"].sum())
        archive = float(sub["dchi2_archive_anchor"].sum())
        cand = float(sub["dchi2_raw"].sum())
        oracle = float(sub["dchi2_oracle_anchor"].sum()) if "dchi2_oracle_anchor" in sub.columns else math.nan
        rows.append(
            {
                "removed_top_k_abs_gap": k,
                "n_remaining": int(len(sub)),
                "removed_names": ";".join(top_names[:k]),
                "clean_total_dchi2": clean,
                "archive_total_dchi2": archive,
                "primary_total_dchi2": cand,
                "oracle_total_dchi2": oracle,
                "recovered_gap_fraction": recovered_fraction(cand, clean, archive),
                "delta_dchi2_vs_archive": cand - archive,
                "delta_dchi2_vs_clean": cand - clean,
            }
        )
    return pd.DataFrame(rows)


def feature_ablation(
    a7: Any,
    a5: pd.DataFrame,
    a7_summary: pd.DataFrame,
    ref: pd.DataFrame,
    pts: Dict[str, Any],
    source_totals: Dict[str, float],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []

    # Start with canonical A7 raw rows, exactly as locked by A7/A8.
    for cand in ["clean_unscaled_kW2_raw", "archive_response_lock_raw", PRIMARY_RAW_CANDIDATE, C200_RAW_CANDIDATE, "A6_c200_linear_kernel_raw"]:
        hit = a7_summary.loc[a7_summary["candidate"].astype(str).eq(cand)].copy()
        if not hit.empty:
            r = hit.iloc[0].to_dict()
            rows.append(
                {
                    "audit_candidate": cand,
                    "source": "A7_locked_raw",
                    "features": str(r.get("source_candidate", cand)),
                    "model_type": "locked",
                    "total_dchi2_raw": float(r["total_dchi2_raw"]),
                    "recovered_gap_fraction_raw": float(r["recovered_gap_fraction_raw"]),
                    "delta_dchi2_vs_archive_raw": float(r["delta_dchi2_vs_archive_raw"]),
                    "negative_chi2_rows": int(r.get("negative_chi2_rows", 0)),
                    "invalid_ratio_rows": int(r.get("invalid_ratio_rows", 0)),
                    "scale_median": float(r.get("scale_median", math.nan)) if pd.notna(r.get("scale_median", math.nan)) else math.nan,
                }
            )

    # Additional D6-only source-side ablations.  These are diagnostics; they do
    # not replace the locked A7 rows.
    ablations = [
        ("D6_logM200_only_linear_raw", [FEATURE_LOGM], "linear"),
        ("D6_c200_only_linear_recompute_raw", [FEATURE_C200], "linear"),
        ("D6_c200_only_loglinear_recompute_raw", [FEATURE_C200], "loglinear"),
        ("D6_logM200_plus_c200_linear_recompute_raw", [FEATURE_LOGM, FEATURE_C200], "linear"),
    ]
    if FEATURE_DKW2 in a5.columns:
        ablations.append(("D6_delta_kW2_only_linear_raw", [FEATURE_DKW2], "linear"))
    if FEATURE_DW2 in a5.columns:
        ablations.append(("D6_delta_W2_only_linear_raw", [FEATURE_DW2], "linear"))

    for label, feats, model_type in ablations:
        try:
            scales = fit_predict_loo(a5, feats, TARGET_SCALE, model_type)
            summary, _ = evaluate_scales_raw(a7, a5, scales, label, "+".join(feats), ref, pts, source_totals)
            rows.append(
                {
                    "audit_candidate": label,
                    "source": "D6_recomputed_raw",
                    "features": "+".join(feats),
                    "model_type": model_type,
                    "total_dchi2_raw": float(summary["total_dchi2_raw"]),
                    "recovered_gap_fraction_raw": float(summary["recovered_gap_fraction_raw"]),
                    "delta_dchi2_vs_archive_raw": float(summary["delta_dchi2_vs_archive_raw"]),
                    "negative_chi2_rows": int(summary.get("negative_chi2_rows", 0)),
                    "invalid_ratio_rows": int(summary.get("invalid_ratio_rows", 0)),
                    "scale_median": float(summary.get("scale_median", math.nan)),
                }
            )
        except Exception as e:
            rows.append(
                {
                    "audit_candidate": label,
                    "source": "D6_recomputed_raw",
                    "features": "+".join(feats),
                    "model_type": model_type,
                    "total_dchi2_raw": math.nan,
                    "recovered_gap_fraction_raw": math.nan,
                    "delta_dchi2_vs_archive_raw": math.nan,
                    "negative_chi2_rows": math.nan,
                    "invalid_ratio_rows": math.nan,
                    "scale_median": math.nan,
                    "error": str(e),
                }
            )

    out = pd.DataFrame(rows)
    if "recovered_gap_fraction_raw" in out.columns:
        out = out.sort_values("recovered_gap_fraction_raw", ascending=False, na_position="last").reset_index(drop=True)
    return out


def build_input_lock(paths: Dict[str, Path]) -> pd.DataFrame:
    rows = []
    for label, p in paths.items():
        rows.append({"input": label, "path": str(p), "exists": bool(p.exists()), "sha256": sha256_path(p)})
    return pd.DataFrame(rows)


def main() -> None:
    a7 = load_a7_module()
    data = prepare_inputs(a7)
    a5 = data["a5"]
    a7_summary = data["a7_summary"]
    a7_per = data["a7_per"]
    ref = data["ref"]
    pts = data["pts"]
    source_totals = data["source_totals"]

    # Canonical A7 locked numbers.
    primary_row = a7_summary.loc[a7_summary["candidate"].astype(str).eq(PRIMARY_RAW_CANDIDATE)]
    c200_row = a7_summary.loc[a7_summary["candidate"].astype(str).eq(C200_RAW_CANDIDATE)]
    if primary_row.empty:
        raise RuntimeError(f"Missing primary A7 row: {PRIMARY_RAW_CANDIDATE}")
    primary = primary_row.iloc[0].to_dict()
    c200 = c200_row.iloc[0].to_dict() if not c200_row.empty else None
    primary_rec = float(primary["recovered_gap_fraction_raw"])
    primary_total = float(primary["total_dchi2_raw"])
    primary_delta_archive = float(primary["delta_dchi2_vs_archive_raw"])
    c200_rec = float(c200["recovered_gap_fraction_raw"]) if c200 else math.nan

    null_df, perm_audit, recomputed_primary_per = run_permutation_null(a7, a5, ref, pts, source_totals, N_PERM, SEED)
    top_df = top_driver_removal(a7_per)
    ablation_df = feature_ablation(a7, a5, a7_summary, ref, pts, source_totals)

    # Use A7 locked primary for final claim, but require the D6 recomputed primary
    # to be close enough to ensure the permutation implementation is comparable.
    recomputed_primary_rec = float(perm_audit["real_recovered_gap_fraction_raw_recomputed"])
    recompute_delta = abs(recomputed_primary_rec - primary_rec)

    perm_pass = (
        float(perm_audit["p_perm_recovered_ge_real"]) <= 0.01
        and primary_rec > float(perm_audit["null_recovered_p95"])
    )
    perm_strong = (
        float(perm_audit["p_perm_recovered_ge_real"]) <= 0.001
        and primary_rec > float(perm_audit["null_recovered_p99"])
    )
    # Top-driver guard: the effect should remain positive after removing the
    # major residual contributors; top-10 is allowed to be weaker but should not
    # flip sign.
    top_map = {int(r["removed_top_k_abs_gap"]): float(r["recovered_gap_fraction"]) for _, r in top_df.iterrows()}
    top_pass = (top_map.get(1, 0.0) > 0.50) and (top_map.get(3, 0.0) > 0.45) and (top_map.get(5, 0.0) > 0.35) and (top_map.get(10, 0.0) > 0.0)
    ablation_pass = math.isfinite(c200_rec) and c200_rec >= 0.65
    no_bad_rows = int(primary.get("negative_chi2_rows", 0)) == 0 and int(primary.get("invalid_ratio_rows", 0)) == 0
    recompute_ok = recompute_delta < 0.05

    if perm_strong and top_pass and ablation_pass and no_bad_rows and recompute_ok:
        status = "D6_KERNEL_FALSIFICATION_LOCK_STRONG_PASS"
        d5d6_status = "D5_D6_AMPLITUDE_KERNEL_RAW_RC_VALIDATED_AND_NULL_TESTED"
        closure = "PASS_FINAL_D6_LOCK"
    elif perm_pass and top_pass and ablation_pass and no_bad_rows and recompute_ok:
        status = "D6_KERNEL_FALSIFICATION_LOCK_PASS"
        d5d6_status = "D5_D6_AMPLITUDE_KERNEL_RAW_RC_VALIDATED_AND_NULL_TESTED"
        closure = "PASS_FINAL_D6_LOCK"
    else:
        status = "D6_KERNEL_FALSIFICATION_LOCK_REVIEW"
        d5d6_status = "D5_D6_KERNEL_VALIDATED_BUT_D6_REVIEW_REQUIRED"
        closure = "REVIEW_BEFORE_FINAL_PAPERII_INSERTION"

    # Source-leakage guard for D6 ablation/permutation features.
    d6_features = [FEATURE_LOGM, FEATURE_C200, FEATURE_DKW2, FEATURE_DW2]
    used_features = [f for f in d6_features if f in a5.columns]
    forbidden_used = [f for f in used_features if any(h.lower() in f.lower() for h in FORBIDDEN_HINTS)]
    leakage_pass = len(forbidden_used) == 0

    null_df.to_csv(OUT_PERM, index=False)
    top_df.to_csv(OUT_TOP, index=False)
    ablation_df.to_csv(OUT_ABL, index=False)

    decisions = pd.DataFrame(
        [
            {
                "item": "D6_status",
                "decision": status,
                "basis": f"p_perm={perm_audit['p_perm_recovered_ge_real']:.6g}; z={perm_audit['z_recovered']:.3f}; top_pass={top_pass}; c200_rec={c200_rec:.6f}; no_bad_rows={no_bad_rows}.",
            },
            {
                "item": "D5D6_status",
                "decision": d5d6_status,
                "basis": "D5-A7 raw RC validation plus D6 permutation/top-driver/ablation audit.",
            },
            {
                "item": "permutation_null",
                "decision": "STRONG_PASS" if perm_strong else ("PASS" if perm_pass else "REVIEW"),
                "basis": json.dumps(perm_audit, sort_keys=True),
            },
            {
                "item": "top_driver_removal",
                "decision": "PASS" if top_pass else "REVIEW",
                "basis": json.dumps(top_map, sort_keys=True),
            },
            {
                "item": "feature_ablation",
                "decision": "PASS" if ablation_pass else "REVIEW",
                "basis": f"c200-only A7 raw recovered fraction={c200_rec:.6f}; primary raw recovered fraction={primary_rec:.6f}.",
            },
            {
                "item": "raw_bad_rows",
                "decision": "PASS" if no_bad_rows else "FAIL",
                "basis": f"primary negative_chi2_rows={primary.get('negative_chi2_rows', 'NA')}; invalid_ratio_rows={primary.get('invalid_ratio_rows', 'NA')}.",
            },
            {
                "item": "d6_leakage_guard",
                "decision": "PASS" if leakage_pass else "REVIEW",
                "basis": f"D6 used source-side features={used_features}; forbidden_used={forbidden_used}.",
            },
            {
                "item": "recomputed_primary_consistency",
                "decision": "PASS" if recompute_ok else "REVIEW",
                "basis": f"A7 locked rec={primary_rec:.6f}; D6 recomputed rec={recomputed_primary_rec:.6f}; abs_delta={recompute_delta:.6f}.",
            },
            {
                "item": "closure_decision",
                "decision": closure,
                "basis": "D6 closes if raw A7 success is non-random under feature permutation, not top-driver-only, and c200 source-side ablation remains strong.",
            },
            {
                "item": "next_gate",
                "decision": "PAPERII_MANUSCRIPT_INSERTION" if closure == "PASS_FINAL_D6_LOCK" else "D6B_REVIEW_NULL_OR_DRIVER_FAILURE",
                "basis": "After D6 pass, insert D5/D6 amplitude wording into Paper II. If review, inspect permutation/top-removal tables first.",
            },
        ]
    )
    decisions.to_csv(OUT_DECISIONS, index=False)

    claims = pd.DataFrame(
        [
            {
                "claim_id": "P2.D6.01",
                "status": "SOURCE_LOCKED",
                "level": "NULL_TEST",
                "claim": "The A6 raw-RC kernel beats a source-feature permutation null.",
                "locked_value": f"p_perm={perm_audit['p_perm_recovered_ge_real']:.6g}; z={perm_audit['z_recovered']:.3f}; null_p95={perm_audit['null_recovered_p95']:.6f}; real_raw={primary_rec:.6f}",
            },
            {
                "claim_id": "P2.D6.02",
                "status": "SOURCE_LOCKED",
                "level": "ROBUSTNESS",
                "claim": "The D5 raw-RC amplitude result is not only a single top-galaxy effect.",
                "locked_value": "; ".join([f"remove_top{k}:rec={top_map.get(k, math.nan):.6f}" for k in [1, 3, 5, 10]]),
            },
            {
                "claim_id": "P2.D6.03",
                "status": "SOURCE_LOCKED",
                "level": "ABLATION",
                "claim": "The effect is dominated by source-side halo-shape/concentration information.",
                "locked_value": f"c200-only raw recovered={c200_rec:.6f}; primary raw recovered={primary_rec:.6f}",
            },
            {
                "claim_id": "P2.D6.04",
                "status": "SOURCE_LOCKED",
                "level": "GUARDRAIL",
                "claim": "Do not claim exact archive or oracle recovery.",
                "locked_value": f"primary delta_archive={primary_delta_archive:.6f}; primary total_dchi2={primary_total:.6f}",
            },
            {
                "claim_id": "P2.D6.05",
                "status": "SOURCE_LOCKED",
                "level": "FINAL_LOCK",
                "claim": "D5/D6 is locked for Paper II with guardrails.",
                "locked_value": closure,
            },
        ]
    )
    claims.to_csv(OUT_CLAIMS, index=False)

    input_lock = build_input_lock(
        {
            "A5 forensic": A5_PATH,
            "A6 model summary": A6_SUMMARY_PATH,
            "A6b resampling": A6B_RESAMPLING_PATH,
            "A6b coefficients": A6B_COEFF_PATH,
            "A7 candidate summary": A7_SUMMARY_PATH,
            "A7 per-galaxy raw RC": A7_PER_GAL_PATH,
            "A8 final lock": A8_LOCK_PATH,
            "A8 numeric locks": A8_NUMERIC_LOCKS_PATH,
            "A8 claim lock": A8_CLAIM_LOCK_PATH,
            "A7 raw evaluator script": A7_SCRIPT_PATH,
            "SPARC MRT": data["raw_paths"]["mrt"],
            "NFW reference": data["raw_paths"]["ref"],
            "active sample": data["raw_paths"]["active"],
        }
    )

    payload = {
        "status": status,
        "d5d6_status": d5d6_status,
        "closure_decision": closure,
        "n_perm": N_PERM,
        "seed": SEED,
        "primary_a7_locked": primary,
        "c200_a7_locked": c200,
        "permutation_audit": perm_audit,
        "top_driver_removal": top_df.to_dict(orient="records"),
        "feature_ablation": ablation_df.to_dict(orient="records"),
        "claims": claims.to_dict(orient="records"),
        "decisions": decisions.to_dict(orient="records"),
        "input_lock": input_lock.to_dict(orient="records"),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    wording = f"""D6 falsification and final amplitude-kernel lock.

After D5-A7, the raw point-level SPARC validation of the A6 row-dependent kernel is subjected to a final D6 falsification audit. The locked A7 primary raw model gives total dchi2 = {primary_total:.6f}, recovered gap fraction = {primary_rec:.6f}, and delta versus the archive-response lock = {primary_delta_archive:.6f}. The c200-only source-side control gives recovered fraction = {c200_rec:.6f}.

D6 permutes the source-side feature vectors across the 38 active galaxies and repeats the same LOO kernel construction. The real kernel has recovered fraction = {primary_rec:.6f}; the permutation null has median = {perm_audit['null_recovered_p50']:.6f}, p95 = {perm_audit['null_recovered_p95']:.6f}, p99 = {perm_audit['null_recovered_p99']:.6f}, empirical p = {perm_audit['p_perm_recovered_ge_real']:.6g}, and z = {perm_audit['z_recovered']:.3f}.

Top-driver removal gives recovered fractions: top1 removed = {top_map.get(1, math.nan):.6f}, top3 removed = {top_map.get(3, math.nan):.6f}, top5 removed = {top_map.get(5, math.nan):.6f}, top10 removed = {top_map.get(10, math.nan):.6f}. This is a robustness audit, not a replacement for the full 38-galaxy result.

The final allowed Paper-II claim is: a source-only row-dependent amplitude kernel, dominated by halo shape/concentration information, recovers most of the clean-to-archive SPARC amplitude gap in raw point-level rotation curves and survives a source-feature permutation/null audit. The forbidden claim remains exact archive or oracle recovery.

D6 status: {status}.
Closure decision: {closure}.
"""
    OUT_WORDING.write_text(wording, encoding="utf-8")

    md: List[str] = []
    md.append("# Paper II / D6 — kernel falsification and final lock\n\n")
    md.append("## Status\n\n")
    md.append(f"**{status}**\n\n")
    md.append(f"D5/D6 status: **{d5d6_status}**\n\n")
    md.append(f"Closure decision: **{closure}**\n\n")

    md.append("## Input source lock\n\n")
    md.append(md_table(input_lock, ["input", "path", "exists", "sha256"]))
    md.append("\n\n")

    md.append("## Locked A7 raw numbers\n\n")
    locked_df = pd.DataFrame(
        [
            {"quantity": "A7 primary total dchi2 raw", "value": primary_total},
            {"quantity": "A7 primary recovered gap fraction raw", "value": primary_rec},
            {"quantity": "A7 primary delta vs archive raw", "value": primary_delta_archive},
            {"quantity": "A7 c200-only recovered gap fraction raw", "value": c200_rec},
            {"quantity": "D6 recomputed primary recovered gap fraction raw", "value": recomputed_primary_rec},
            {"quantity": "D6 recomputed-vs-A7 recovered abs delta", "value": recompute_delta},
        ]
    )
    md.append(locked_df.to_markdown(index=False))
    md.append("\n\n")

    md.append("## D6-A1 permutation null\n\n")
    perm_df = pd.DataFrame([perm_audit])
    md.append(perm_df.to_markdown(index=False))
    md.append("\n\n")

    md.append("## D6-A2 top-driver removal\n\n")
    md.append(md_table(top_df))
    md.append("\n\n")

    md.append("## D6-A3 feature ablation\n\n")
    show_ab_cols = [
        "audit_candidate",
        "source",
        "features",
        "model_type",
        "total_dchi2_raw",
        "recovered_gap_fraction_raw",
        "delta_dchi2_vs_archive_raw",
        "negative_chi2_rows",
        "invalid_ratio_rows",
        "scale_median",
    ]
    md.append(md_table(ablation_df, [c for c in show_ab_cols if c in ablation_df.columns]))
    md.append("\n\n")

    md.append("## Claim source lock\n\n")
    md.append(md_table(claims))
    md.append("\n\n")

    md.append("## Decisions\n\n")
    md.append(md_table(decisions))
    md.append("\n\n")

    md.append("## Paper-II wording block\n\n")
    md.append(wording)
    md.append("\n")

    md.append("## Interpretation lock\n\n")
    if closure == "PASS_FINAL_D6_LOCK":
        md.append(
            "D6 closes the statistical falsification layer for the D5 amplitude kernel. "
            "The Paper-II amplitude statement can now cite D5-A7 for raw RC validation and D6 for the null/top-driver/ablation guardrails. "
            "The wording remains bounded: most of the clean-to-archive gap is recovered, not the exact archive/oracle response.\n"
        )
    else:
        md.append(
            "D6 requires review before final Paper-II insertion. Keep the D5-A7 raw validation, but do not promote the D6 null/falsification language until the review item is resolved.\n"
        )

    OUT_MD.write_text("".join(md), encoding="utf-8")

    print(f"Saved: {OUT_MD}")
    print(f"Saved: {OUT_JSON}")
    print(f"Saved: {OUT_PERM}")
    print(f"Saved: {OUT_TOP}")
    print(f"Saved: {OUT_ABL}")
    print(f"Saved: {OUT_DECISIONS}")
    print(f"Saved: {OUT_CLAIMS}")
    print(f"Saved: {OUT_WORDING}")
    print(f"Status: {status}")
    print(f"D5/D6 status: {d5d6_status}")
    print(f"Closure decision: {closure}")
    print(f"Primary raw recovered gap fraction: {primary_rec}")
    print(f"Primary raw total dchi2: {primary_total}")
    print(f"Permutation p-value: {perm_audit['p_perm_recovered_ge_real']}")
    print(f"Permutation z: {perm_audit['z_recovered']}")
    print(f"Top removal recovered fractions: {top_map}")
    print(f"c200-only raw recovered gap fraction: {c200_rec}")


if __name__ == "__main__":
    main()
