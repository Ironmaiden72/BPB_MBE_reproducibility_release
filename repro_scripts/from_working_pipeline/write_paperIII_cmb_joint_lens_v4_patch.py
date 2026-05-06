#!/usr/bin/env python3
"""
Paper III CMB P3-A1b — patch joint ACT+Planck lensing v4 extraction.

Purpose
-------
Fixes the P3.JOINTLENS.01 source-match extraction by accepting the actual
Planck lensing column name used by the v4 grid:

    chi2_Planck_lensing

instead of assuming only:

    chi2_Planck

Safety
------
- Reads only outputs/cmb/*.csv.
- Writes only outputs/cmb_paperIII/*.
- Does not touch outputs/ns, src/ns, or NS checkpoints.
- Does not run ACT/Planck likelihoods; this is extraction/source-lock only.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


GRID_DEFAULT = Path("outputs/cmb/cmb15a_v4_joint_act_planck_scalar_amp_map1782.grid.csv")
REFS_DEFAULT = Path("outputs/cmb/cmb15a_v4_joint_act_planck_scalar_amp_map1782.refs.csv")
OUTDIR_DEFAULT = Path("outputs/cmb_paperIII")


EXPECTED = {
    "act_best": {
        "amp_phi": 1.043,
        "Sigma_eff": 1.021274,
        "chi2_ACT": 14.153599,
        "chi2_Planck_lensing": 10.198317,
        "chi2_joint": 24.351916,
    },
    "planck_best": {
        "amp_phi": 1.075,
        "Sigma_eff": 1.036822,
        "chi2_ACT": 15.908702,
        "chi2_Planck_lensing": 8.742460,
        "chi2_joint": 24.651162,
    },
    "joint_best": {
        "amp_phi": 1.057,
        "Sigma_eff": 1.028105,
        "chi2_ACT": 14.478027,
        "chi2_Planck_lensing": 9.243264,
        "chi2_joint": 23.721292,
    },
    "separate_minima_sum": 22.896060,
    "joint_scalar_penalty": 0.825232,
}


@dataclass
class BestRow:
    label: str
    row_index: int
    amp_phi: float | None
    Sigma_eff: float | None
    chi2_ACT: float
    chi2_Planck_lensing: float
    chi2_joint: float


@dataclass
class ExtractionResult:
    status: str
    timestamp_utc: str
    grid_file: str
    refs_file: str | None
    n_rows: int
    columns_detected: dict[str, str | None]
    act_best: BestRow
    planck_best: BestRow
    joint_best: BestRow
    separate_minima_sum: float
    joint_scalar_penalty: float
    expected_value_check: list[dict[str, object]]
    gate_decisions: dict[str, str]


def pick_column(columns: Iterable[str], candidates: Iterable[str], *, required: bool = True) -> str | None:
    cols = list(columns)
    exact = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in exact:
            return exact[cand.lower()]
    # fallback: case-insensitive substring match, useful for minor naming variants
    for cand in candidates:
        token = cand.lower()
        for col in cols:
            if token in col.lower():
                return col
    if required:
        raise KeyError(f"Could not find any of {list(candidates)} in columns: {cols}")
    return None


def finite_numeric(df: pd.DataFrame, col: str) -> pd.Series:
    s = pd.to_numeric(df[col], errors="coerce")
    if not s.notna().any():
        raise ValueError(f"Column {col!r} has no finite numeric values")
    return s


def row_to_best(label: str, df: pd.DataFrame, idx: int, cols: dict[str, str | None]) -> BestRow:
    row = df.loc[idx]

    def val(col: str | None) -> float | None:
        if col is None:
            return None
        x = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
        return None if pd.isna(x) else float(x)

    return BestRow(
        label=label,
        row_index=int(idx),
        amp_phi=val(cols.get("amp_phi")),
        Sigma_eff=val(cols.get("Sigma_eff")),
        chi2_ACT=float(val(cols["chi2_ACT"])),
        chi2_Planck_lensing=float(val(cols["chi2_Planck_lensing"])),
        chi2_joint=float(val(cols["chi2_joint"])),
    )


def check_expected(result: dict[str, object], tol: float = 7e-4) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []

    def add(name: str, got: float | None, expected: float):
        if got is None:
            ok = False
            diff = None
        else:
            diff = float(got) - float(expected)
            ok = abs(diff) <= tol
        checks.append({"quantity": name, "got": got, "expected": expected, "diff": diff, "ok": ok})

    for block in ["act_best", "planck_best", "joint_best"]:
        got_block = result[block]
        if isinstance(got_block, BestRow):
            got_dict = asdict(got_block)
        else:
            got_dict = got_block
        for key, exp_val in EXPECTED[block].items():
            # Some old summaries call the same column chi2_Planck; in this patch output it is normalized.
            add(f"{block}.{key}", got_dict.get(key), exp_val)

    add("separate_minima_sum", result["separate_minima_sum"], EXPECTED["separate_minima_sum"])
    add("joint_scalar_penalty", result["joint_scalar_penalty"], EXPECTED["joint_scalar_penalty"])
    return checks


def format_float(x: float | None, ndigits: int = 9) -> str:
    return "" if x is None else f"{x:.{ndigits}f}"


def write_outputs(res: ExtractionResult, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for best in [res.act_best, res.planck_best, res.joint_best]:
        summary_rows.append(asdict(best))
    summary_rows.append(
        {
            "label": "separate_minima_sum",
            "row_index": "",
            "amp_phi": "",
            "Sigma_eff": "",
            "chi2_ACT": res.act_best.chi2_ACT,
            "chi2_Planck_lensing": res.planck_best.chi2_Planck_lensing,
            "chi2_joint": res.separate_minima_sum,
        }
    )
    summary_rows.append(
        {
            "label": "joint_scalar_penalty",
            "row_index": "",
            "amp_phi": "",
            "Sigma_eff": "",
            "chi2_ACT": "",
            "chi2_Planck_lensing": "",
            "chi2_joint": res.joint_scalar_penalty,
        }
    )
    summary_csv = outdir / "paperIII_cmb_joint_lens_v4_patch.csv"
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    checks_csv = outdir / "paperIII_cmb_joint_lens_v4_expected_checks.csv"
    pd.DataFrame(res.expected_value_check).to_csv(checks_csv, index=False)

    json_path = outdir / "paperIII_cmb_joint_lens_v4_patch.json"
    json_path.write_text(json.dumps(asdict(res), indent=2, ensure_ascii=False), encoding="utf-8")

    lock_row = pd.DataFrame(
        [
            {
                "claim id": "P3.JOINTLENS.01",
                "status": "SOURCE_MATCHED_PATCHED_P3A1B",
                "level": "JOINT_ACT_PLANCK_SCALAR_LENSING_DIAGNOSTIC",
                "claim": "A single scalar Weyl-auto amplitude gives a coherent joint ACT DR6 + Planck CMB-marginalized lensing closure; the joint scalar best fit is close to the separate ACT-only and Planck-only minima, with a small penalty relative to independent amplitudes.",
                "source files": f"{res.grid_file}; {res.refs_file or ''}",
                "paper wording": "Joint lensing is a scalar-amplitude diagnostic closure, not a full combined official ACT+Planck cosmological likelihood. Canonical values: ACT-only amp_phi=1.043, Planck-only amp_phi=1.075, joint scalar amp_phi=1.057, joint chi2=23.721292, separate-minima sum=22.896060, scalar penalty=0.825232.",
            }
        ]
    )
    lock_patch_csv = outdir / "paperIII_cmb_claim_source_lock_p3a1b_jointlens_patch.csv"
    lock_row.to_csv(lock_patch_csv, index=False)

    md_path = outdir / "paperIII_cmb_joint_lens_v4_patch.md"
    checks_ok = all(bool(c["ok"]) for c in res.expected_value_check)
    md = f"""# Paper III CMB P3-A1b — joint lensing v4 extraction patch

## Status

**{res.status}**

This gate fixes the source-match extraction for `P3.JOINTLENS.01` by using `{res.columns_detected['chi2_Planck_lensing']}` as the Planck lensing column.

Safety: read-only on `outputs/cmb`; writes only to `outputs/cmb_paperIII`; no ACT/Planck likelihood run; no NS touch.

## Source files

- Grid: `{res.grid_file}`
- Refs: `{res.refs_file or 'not found / not used'}`

Rows scanned: **{res.n_rows}**

## Canonical extracted values

| selector | amp_phi | Sigma_eff | chi2_ACT | chi2_Planck_lensing | chi2_joint |
|---|---:|---:|---:|---:|---:|
| ACT-only best | {format_float(res.act_best.amp_phi, 6)} | {format_float(res.act_best.Sigma_eff, 9)} | {format_float(res.act_best.chi2_ACT, 6)} | {format_float(res.act_best.chi2_Planck_lensing, 6)} | {format_float(res.act_best.chi2_joint, 6)} |
| Planck-only best | {format_float(res.planck_best.amp_phi, 6)} | {format_float(res.planck_best.Sigma_eff, 9)} | {format_float(res.planck_best.chi2_ACT, 6)} | {format_float(res.planck_best.chi2_Planck_lensing, 6)} | {format_float(res.planck_best.chi2_joint, 6)} |
| Joint scalar best | {format_float(res.joint_best.amp_phi, 6)} | {format_float(res.joint_best.Sigma_eff, 9)} | {format_float(res.joint_best.chi2_ACT, 6)} | {format_float(res.joint_best.chi2_Planck_lensing, 6)} | {format_float(res.joint_best.chi2_joint, 6)} |

Separate minima sum: **{res.separate_minima_sum:.6f}**  
Joint scalar penalty: **{res.joint_scalar_penalty:.6f}**

Expected-value check: **{'PASS' if checks_ok else 'CHECK'}**

## Gate decisions

| decision | value |
|---|---|
"""
    for k, v in res.gate_decisions.items():
        md += f"| {k} | {v} |\n"
    md += f"""

## Outputs

- `{summary_csv}`
- `{checks_csv}`
- `{json_path}`
- `{lock_patch_csv}`
- `{md_path}`
"""
    md_path.write_text(md, encoding="utf-8")

    print(f"Saved: {summary_csv}")
    print(f"Saved: {checks_csv}")
    print(f"Saved: {json_path}")
    print(f"Saved: {lock_patch_csv}")
    print(f"Saved: {md_path}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Patch Paper III P3.JOINTLENS.01 extraction for v4 grid columns.")
    ap.add_argument("--grid", type=Path, default=GRID_DEFAULT)
    ap.add_argument("--refs", type=Path, default=REFS_DEFAULT)
    ap.add_argument("--outdir", type=Path, default=OUTDIR_DEFAULT)
    args = ap.parse_args()

    if not args.grid.exists():
        raise FileNotFoundError(f"Missing grid CSV: {args.grid}")

    df = pd.read_csv(args.grid)
    if df.empty:
        raise ValueError(f"Grid CSV is empty: {args.grid}")

    cols = {
        "amp_phi": pick_column(df.columns, ["amp_phi", "amp_PhiPhi", "amp_phi_scalar", "amp"], required=False),
        "Sigma_eff": pick_column(df.columns, ["Sigma_eff", "sigma_eff", "Sigma", "sigma"], required=False),
        "chi2_ACT": pick_column(df.columns, ["chi2_ACT", "chi2_act", "act_chi2"], required=True),
        "chi2_Planck_lensing": pick_column(
            df.columns,
            ["chi2_Planck_lensing", "chi2_planck_lensing", "chi2_Planck", "chi2_planck", "planck_chi2"],
            required=True,
        ),
        "chi2_joint": pick_column(df.columns, ["chi2_joint", "joint_chi2", "chi2_total", "total_chi2"], required=False),
    }

    act = finite_numeric(df, cols["chi2_ACT"])
    planck = finite_numeric(df, cols["chi2_Planck_lensing"])

    if cols["chi2_joint"] is None:
        df = df.copy()
        df["chi2_joint__computed_ACT_plus_Planck_lensing"] = act + planck
        cols["chi2_joint"] = "chi2_joint__computed_ACT_plus_Planck_lensing"
    joint = finite_numeric(df, cols["chi2_joint"])

    act_idx = int(act.idxmin())
    planck_idx = int(planck.idxmin())
    joint_idx = int(joint.idxmin())

    result_dict = {
        "act_best": row_to_best("act_best", df, act_idx, cols),
        "planck_best": row_to_best("planck_best", df, planck_idx, cols),
        "joint_best": row_to_best("joint_best", df, joint_idx, cols),
    }
    separate_minima_sum = result_dict["act_best"].chi2_ACT + result_dict["planck_best"].chi2_Planck_lensing
    joint_scalar_penalty = result_dict["joint_best"].chi2_joint - separate_minima_sum
    result_dict["separate_minima_sum"] = float(separate_minima_sum)
    result_dict["joint_scalar_penalty"] = float(joint_scalar_penalty)

    checks = check_expected(result_dict)
    checks_pass = all(bool(c["ok"]) for c in checks)

    res = ExtractionResult(
        status="P3_A1B_JOINT_LENS_V4_EXTRACTION_PATCH_PASS" if checks_pass else "P3_A1B_JOINT_LENS_V4_EXTRACTION_PATCH_CHECK_VALUES",
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        grid_file=str(args.grid),
        refs_file=str(args.refs) if args.refs.exists() else None,
        n_rows=int(len(df)),
        columns_detected=cols,
        act_best=result_dict["act_best"],
        planck_best=result_dict["planck_best"],
        joint_best=result_dict["joint_best"],
        separate_minima_sum=float(separate_minima_sum),
        joint_scalar_penalty=float(joint_scalar_penalty),
        expected_value_check=checks,
        gate_decisions={
            "P3.JOINTLENS.01 extraction": "PASS" if checks_pass else "CHECK_NUMBERS",
            "chi2_Planck_lensing column detected": "PASS",
            "planck_best found": "PASS",
            "joint_best found": "PASS",
            "source-match corrected": "PASS" if checks_pass else "PARTIAL",
            "likelihood run": "NO",
            "NS touched": "NO",
            "write scope": "outputs/cmb_paperIII only",
            "next gate": "P3-A2 canonical table choice / figure inventory",
        },
    )

    write_outputs(res, args.outdir)
    return 0 if checks_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
