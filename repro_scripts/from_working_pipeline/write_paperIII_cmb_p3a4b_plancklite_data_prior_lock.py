#!/usr/bin/env python3
"""
Paper III CMB P3-A4b — Planck-lite data/prior decomposition source lock.

Purpose
-------
Lock/source-match P3.CLAIM.PLANCKLITE.02:
  the residual penalty of the strict face-minus Planck-lite fine-ns scan is
  calibration-normalization driven, not acoustic-shape driven.

This script is non-destructive for NS work:
  - reads only outputs/cmb/* and/or outputs/cmb_paperIII/*
  - writes only outputs/cmb_paperIII/*
  - does not run ACT/Planck likelihoods
  - does not touch outputs/ns, src/ns, or checkpoints

It prefers raw fine-scan JSON files in outputs/cmb/cmb11d_ns_fine_scan/.
If they are not present, it falls back to the flattened audit file
outputs/cmb_paperIII/paperIII_cmb_json_key_values.csv.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

A_PLANCK_SIGMA = 0.0025
EXPECTED = {
    "planck_ns": 0.970,
    "planck_chi2": 791.3148594982862,
    "strict_ns": 0.976,
    "strict_chi2": 817.7914377494669,
    "delta_total": 26.47657825118074,
    "delta_data": -5.054032720707028,
    "delta_prior": 31.530610971887803,
}
TOL = {
    "ns": 5e-7,
    "chi2": 5e-6,
    "delta": 5e-6,
}


@dataclass
class FineNsRow:
    source_file: str
    ns: float
    planck_A: float
    planck_chi2: float
    strict_A: float
    strict_chi2: float
    vanilla_A: Optional[float] = None
    vanilla_chi2: Optional[float] = None
    zbest_A: Optional[float] = None
    zbest_chi2: Optional[float] = None
    verdict: str = ""

    @property
    def planck_prior(self) -> float:
        return ((self.planck_A - 1.0) / A_PLANCK_SIGMA) ** 2

    @property
    def planck_data(self) -> float:
        return self.planck_chi2 - self.planck_prior

    @property
    def strict_prior(self) -> float:
        return ((self.strict_A - 1.0) / A_PLANCK_SIGMA) ** 2

    @property
    def strict_data(self) -> float:
        return self.strict_chi2 - self.strict_prior

    def to_record(self) -> Dict[str, Any]:
        d = asdict(self)
        d.update(
            {
                "planck_prior_A_chi2": self.planck_prior,
                "planck_data_TTTEEE_chi2": self.planck_data,
                "strict_prior_A_chi2": self.strict_prior,
                "strict_data_TTTEEE_chi2": self.strict_data,
                "strict_minus_planck_total_chi2": self.strict_chi2 - self.planck_chi2,
                "strict_minus_planck_data_chi2": self.strict_data - self.planck_data,
                "strict_minus_planck_prior_chi2": self.strict_prior - self.planck_prior,
            }
        )
        return d


def _parse_ns_from_path(path: str) -> Optional[float]:
    m = re.search(r"ns_(\d+)p(\d+)", path)
    if not m:
        return None
    return float(f"{m.group(1)}.{m.group(2)}")


def _num(x: Any) -> float:
    if isinstance(x, (int, float)):
        return float(x)
    return float(str(x))


def load_from_raw_json(repo: Path) -> List[FineNsRow]:
    rows: List[FineNsRow] = []
    scan_dir = repo / "outputs" / "cmb" / "cmb11d_ns_fine_scan"
    if not scan_dir.exists():
        return rows
    for p in sorted(scan_dir.glob("planck_lite_ns_*.json")):
        ns = _parse_ns_from_path(str(p))
        if ns is None:
            continue
        data = json.loads(p.read_text(encoding="utf-8"))
        cases = data.get("cases", {})
        try:
            planck = cases["Planck_PR4"]
            strict = cases["C_s_surv0"]
        except KeyError:
            continue
        rows.append(
            FineNsRow(
                source_file=str(p.relative_to(repo)).replace("\\", "/"),
                ns=ns,
                planck_A=_num(planck["A_planck_best"]),
                planck_chi2=_num(planck["chi2_planck"]),
                strict_A=_num(strict["A_planck_best"]),
                strict_chi2=_num(strict["chi2_planck"]),
                vanilla_A=_num(cases.get("A_vanilla", {}).get("A_planck_best", "nan"))
                if "A_vanilla" in cases
                else None,
                vanilla_chi2=_num(cases.get("A_vanilla", {}).get("chi2_planck", "nan"))
                if "A_vanilla" in cases
                else None,
                zbest_A=_num(cases.get("C_Zbest", {}).get("A_planck_best", "nan"))
                if "C_Zbest" in cases
                else None,
                zbest_chi2=_num(cases.get("C_Zbest", {}).get("chi2_planck", "nan"))
                if "C_Zbest" in cases
                else None,
                verdict=str(data.get("verdict", "")),
            )
        )
    return rows


def load_from_flattened_audit(repo: Path) -> List[FineNsRow]:
    p = repo / "outputs" / "cmb_paperIII" / "paperIII_cmb_json_key_values.csv"
    if not p.exists():
        # Also support running directly on an unpacked cmb_paperIII directory.
        p = repo / "cmb_paperIII" / "paperIII_cmb_json_key_values.csv"
    if not p.exists():
        return []

    grouped: Dict[str, Dict[str, Any]] = {}
    with p.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            file = r.get("file", "")
            if "cmb11d_ns_fine_scan" not in file or not file.endswith(".json"):
                continue
            grouped.setdefault(file, {})[r["key"]] = r["value"]

    rows: List[FineNsRow] = []
    for file, d in sorted(grouped.items()):
        ns = _parse_ns_from_path(file)
        if ns is None:
            continue
        needed = [
            "cases.Planck_PR4.A_planck_best",
            "cases.Planck_PR4.chi2_planck",
            "cases.C_s_surv0.A_planck_best",
            "cases.C_s_surv0.chi2_planck",
        ]
        if any(k not in d for k in needed):
            continue
        rows.append(
            FineNsRow(
                source_file=file,
                ns=ns,
                planck_A=_num(d["cases.Planck_PR4.A_planck_best"]),
                planck_chi2=_num(d["cases.Planck_PR4.chi2_planck"]),
                strict_A=_num(d["cases.C_s_surv0.A_planck_best"]),
                strict_chi2=_num(d["cases.C_s_surv0.chi2_planck"]),
                vanilla_A=_num(d["cases.A_vanilla.A_planck_best"])
                if "cases.A_vanilla.A_planck_best" in d
                else None,
                vanilla_chi2=_num(d["cases.A_vanilla.chi2_planck"])
                if "cases.A_vanilla.chi2_planck" in d
                else None,
                zbest_A=_num(d["cases.C_Zbest.A_planck_best"])
                if "cases.C_Zbest.A_planck_best" in d
                else None,
                zbest_chi2=_num(d["cases.C_Zbest.chi2_planck"])
                if "cases.C_Zbest.chi2_planck" in d
                else None,
                verdict=str(d.get("verdict", "")),
            )
        )
    return rows


def load_rows(repo: Path) -> tuple[List[FineNsRow], str]:
    rows = load_from_raw_json(repo)
    if rows:
        return rows, "raw_json"
    rows = load_from_flattened_audit(repo)
    if rows:
        return rows, "flattened_audit"
    raise FileNotFoundError(
        "Could not find fine ns scan rows in outputs/cmb/cmb11d_ns_fine_scan or "
        "outputs/cmb_paperIII/paperIII_cmb_json_key_values.csv"
    )


def close(a: float, b: float, tol: float) -> bool:
    return abs(float(a) - float(b)) <= tol


def main() -> None:
    repo = Path.cwd()
    out = repo / "outputs" / "cmb_paperIII"
    if not out.exists() and (repo / "cmb_paperIII").exists():
        out = repo / "cmb_paperIII"
    out.mkdir(parents=True, exist_ok=True)

    rows, load_mode = load_rows(repo)
    rows = sorted(rows, key=lambda r: r.ns)
    planck_best = min(rows, key=lambda r: r.planck_chi2)
    strict_best = min(rows, key=lambda r: r.strict_chi2)

    delta_total = strict_best.strict_chi2 - planck_best.planck_chi2
    delta_data = strict_best.strict_data - planck_best.planck_data
    delta_prior = strict_best.strict_prior - planck_best.planck_prior

    checks = [
        ("planck_best_ns", planck_best.ns, EXPECTED["planck_ns"], TOL["ns"]),
        ("planck_best_chi2", planck_best.planck_chi2, EXPECTED["planck_chi2"], TOL["chi2"]),
        ("strict_best_ns", strict_best.ns, EXPECTED["strict_ns"], TOL["ns"]),
        ("strict_best_chi2", strict_best.strict_chi2, EXPECTED["strict_chi2"], TOL["chi2"]),
        ("delta_total", delta_total, EXPECTED["delta_total"], TOL["delta"]),
        ("delta_data", delta_data, EXPECTED["delta_data"], TOL["delta"]),
        ("delta_prior", delta_prior, EXPECTED["delta_prior"], TOL["delta"]),
    ]
    check_records = [
        {"check": name, "value": value, "expected": expected, "tol": tol, "pass": close(value, expected, tol)}
        for name, value, expected, tol in checks
    ]
    all_pass = all(r["pass"] for r in check_records)
    status = "P3_A4B_PLANCKLITE_DATA_PRIOR_LOCK_PASS" if all_pass else "P3_A4B_PLANCKLITE_DATA_PRIOR_LOCK_FAIL"

    # Scan rows CSV
    row_records = [r.to_record() for r in rows]
    rows_csv = out / "paperIII_cmb_p3a4b_plancklite_fine_ns_rows.csv"
    with rows_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row_records[0].keys()))
        writer.writeheader()
        writer.writerows(row_records)

    # Checks CSV
    checks_csv = out / "paperIII_cmb_p3a4b_plancklite_data_prior_checks.csv"
    with checks_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(check_records[0].keys()))
        writer.writeheader()
        writer.writerows(check_records)

    summary = {
        "status": status,
        "load_mode": load_mode,
        "n_fine_ns_rows": len(rows),
        "claim_id": "P3.CLAIM.PLANCKLITE.02",
        "source_files": sorted({planck_best.source_file, strict_best.source_file}),
        "A_planck_prior_sigma": A_PLANCK_SIGMA,
        "planck_best": {
            "source_file": planck_best.source_file,
            "ns": planck_best.ns,
            "A_planck": planck_best.planck_A,
            "chi2_total": planck_best.planck_chi2,
            "chi2_data_TTTEEE": planck_best.planck_data,
            "chi2_A_planck_prior": planck_best.planck_prior,
        },
        "strict_best": {
            "source_file": strict_best.source_file,
            "ns": strict_best.ns,
            "A_planck": strict_best.strict_A,
            "chi2_total": strict_best.strict_chi2,
            "chi2_data_TTTEEE": strict_best.strict_data,
            "chi2_A_planck_prior": strict_best.strict_prior,
        },
        "deltas_strict_minus_planck": {
            "delta_total_chi2": delta_total,
            "delta_data_TTTEEE_chi2": delta_data,
            "delta_A_planck_prior_chi2": delta_prior,
        },
        "interpretation": (
            "The strict face-minus BPB dictionary has a better Planck-lite TTTEEE spectral-shape/data term "
            "than the Planck sanity best along this fine-ns diagnostic, but is penalized by the absolute "
            "A_Planck calibration-normalization prior."
        ),
        "caveat": (
            "Planck-lite TTTEEE diagnostic only; not a full official Planck nuisance/foreground MCMC. "
            "This locks the audit wording, not a final Planck model-selection claim."
        ),
        "ns_touched": False,
        "likelihood_run": False,
        "checks": check_records,
    }

    json_out = out / "paperIII_cmb_p3a4b_plancklite_data_prior_lock.json"
    json_out.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    claim_record = {
        "claim_id": "P3.CLAIM.PLANCKLITE.02",
        "claim_level": "PLANCK_LITE_FINE_NS_DATA_PRIOR_DECOMPOSITION",
        "status": "SOURCE_MATCHED" if all_pass else "CHECK_FAILED",
        "claim": "Planck-lite residual penalty of strict face-minus is calibration-normalization driven, not acoustic-shape driven.",
        "source_files": json.dumps(summary["source_files"]),
        "exact_numbers": json.dumps(
            {
                "Planck_PR4_best_ns": planck_best.ns,
                "Planck_PR4_best_chi2_total": planck_best.planck_chi2,
                "Planck_PR4_best_chi2_data_TTTEEE": planck_best.planck_data,
                "Planck_PR4_best_chi2_A_prior": planck_best.planck_prior,
                "C_s_surv0_best_ns": strict_best.ns,
                "C_s_surv0_best_chi2_total": strict_best.strict_chi2,
                "C_s_surv0_best_chi2_data_TTTEEE": strict_best.strict_data,
                "C_s_surv0_best_chi2_A_prior": strict_best.strict_prior,
                "delta_total": delta_total,
                "delta_data_TTTEEE": delta_data,
                "delta_A_planck_prior": delta_prior,
            },
            sort_keys=True,
        ),
        "paper_wording": (
            "In the Planck-lite fine-ns diagnostic, the residual strict-vs-Planck penalty is not an acoustic-shape failure: "
            f"Δχ²_data(TTTEEE)={delta_data:.3f}, while Δχ²_APlanck={delta_prior:.3f}, giving Δχ²_total={delta_total:.3f}. "
            "The strict branch improves the TT/TE/EE shape term but is penalized by absolute A_Planck normalization."
        ),
    }
    claim_csv = out / "paperIII_cmb_claim_source_lock_p3a4b_plancklite02.csv"
    with claim_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(claim_record.keys()))
        writer.writeheader()
        writer.writerow(claim_record)

    final_audit = {
        "status": "P3_CMB_AUDIT_CLOSED_PRE_LATEX" if all_pass else "P3_CMB_AUDIT_NOT_CLOSED_CHECK_FAILED",
        "gates_closed": [
            "P3-A0 non-destructive audit",
            "P3-A1 source-match claims",
            "P3-A1b joint lensing v4 patch",
            "P3-A2 canonical tables and figure inventory",
            "P3-A3 wording lock",
            "P3-A4 canonical figure generation",
            "P3-A4b Planck-lite data/prior decomposition source lock",
        ],
        "remaining_cmb_items": [
            "Insert locked tables/figures/wording into LaTeX",
            "Full official Planck MCMC remains a later robustness gate, not a missing CMB-theory derivation",
        ],
        "ns_touched": False,
        "likelihood_run": False,
    }
    final_json = out / "paperIII_cmb_final_audit_lock_pre_latex.json"
    final_json.write_text(json.dumps(final_audit, indent=2, sort_keys=True), encoding="utf-8")

    md = out / "paperIII_cmb_p3a4b_plancklite_data_prior_lock.md"
    md.write_text(
        "# Paper III CMB P3-A4b — Planck-lite data/prior decomposition lock\n\n"
        f"## Status\n\n**{status}**\n\n"
        "## Claim\n\n"
        "`P3.CLAIM.PLANCKLITE.02`: the Planck-lite residual penalty of the strict face-minus branch is calibration-normalization driven, not acoustic-shape driven.\n\n"
        "## Source files\n\n"
        + "\n".join(f"- `{s}`" for s in summary["source_files"])
        + "\n\n## Best rows\n\n"
        f"- Planck sanity best: ns={planck_best.ns:.3f}, A_Planck={planck_best.planck_A:.12f}, "
        f"chi2_total={planck_best.planck_chi2:.12f}, chi2_data={planck_best.planck_data:.12f}, "
        f"chi2_prior={planck_best.planck_prior:.12f}\n"
        f"- Strict C_s_surv0 best: ns={strict_best.ns:.3f}, A_Planck={strict_best.strict_A:.12f}, "
        f"chi2_total={strict_best.strict_chi2:.12f}, chi2_data={strict_best.strict_data:.12f}, "
        f"chi2_prior={strict_best.strict_prior:.12f}\n\n"
        "## Decomposition\n\n"
        f"- Δχ²_total = {delta_total:.12f}\n"
        f"- Δχ²_data(TTTEEE) = {delta_data:.12f}\n"
        f"- Δχ²_APlanck_prior = {delta_prior:.12f}\n\n"
        "## Locked interpretation\n\n"
        "The strict branch improves the Planck-lite TT/TE/EE spectral-shape/data term, "
        "but is penalized by the absolute A_Planck calibration-normalization prior. "
        "This strengthens the Paper III wording from 'survives' to 'better shape, calibration-normalization penalty'.\n\n"
        "## Caveat\n\n"
        "Planck-lite TTTEEE diagnostic only; not a full official Planck nuisance/foreground MCMC.\n\n"
        "## Isolation\n\n"
        "NS touched = NO. Likelihood run = NO.\n",
        encoding="utf-8",
    )

    final_md = out / "paperIII_cmb_final_audit_lock_pre_latex.md"
    final_md.write_text(
        "# Paper III CMB — final audit lock before LaTeX\n\n"
        f"## Status\n\n**{final_audit['status']}**\n\n"
        "## Closed gates\n\n"
        + "\n".join(f"- {g}" for g in final_audit["gates_closed"])
        + "\n\n## Remaining CMB items\n\n"
        + "\n".join(f"- {x}" for x in final_audit["remaining_cmb_items"])
        + "\n\n## Isolation\n\nNS touched = NO. Likelihood run = NO.\n",
        encoding="utf-8",
    )

    print(status)
    print(f"Saved: {md}")
    print(f"Saved: {json_out}")
    print(f"Saved: {rows_csv}")
    print(f"Saved: {checks_csv}")
    print(f"Saved: {claim_csv}")
    print(f"Saved: {final_md}")
    print(f"Saved: {final_json}")


if __name__ == "__main__":
    main()
