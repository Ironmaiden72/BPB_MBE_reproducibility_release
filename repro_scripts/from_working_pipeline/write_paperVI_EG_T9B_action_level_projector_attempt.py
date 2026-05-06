#!/usr/bin/env python3
"""
Paper VI / E_G T9B — action-level projector derivation attempt

Non-destructive gate.
- No likelihood run.
- No NS access by design.
- Writes only outputs/paperVI_EG_theory/.

Purpose
-------
Try to promote the T1--T9A projected diagnostic closure

    r_Wdelta = sqrt(a_bg)

from a conditional projected theorem to an action-level HR/BPB perturbation
projector derivation. The strict target is

    <P_m W, P_m W>_EG / <W,W>_EG = a_bg.

This script does NOT declare success unless an explicit HR/BPB perturbation
projector/mixing object is found or supplied. In its default mode, it derives
symbolic closure requirements and inventories the repo for required theoretical
objects, then returns one of:

    P6_T9B_ACTION_LEVEL_PROJECTOR_CLOSED
    P6_T9B_ACTION_LEVEL_NOT_CLOSED_REQUIRE_EXPLICIT_HR_PROJECTOR

The expected honest outcome is usually the second status unless the repo already
contains a usable HR/BPB perturbation matrix/eigenbasis.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Any, Tuple

OUTDIR = Path("outputs/paperVI_EG_theory")

# Avoid heavy or irrelevant tree scans.
SEARCH_ROOTS = [
    Path("papers"),
    Path("paper0"),
    Path("paperI"),
    Path("paperVI"),
    Path("theory"),
    Path("src"),
    Path("runs"),
    Path("pipeline"),
    Path("outputs"),
]

SKIP_PATTERNS = [
    "/outputs/ns/", "/src/ns/", "/eos_", "/.git/", "/.venv/", "/venv/", "/__pycache__/",
]

TEXT_EXTS = {".tex", ".md", ".py", ".json", ".yaml", ".yml", ".txt", ".csv"}

TOKEN_GROUPS = {
    "hr_bigravity_action": ["Hassan", "Rosen", "bigravity", "HR", "beta_n", "bimetric", "g_{\\mu", "f_{\\mu"],
    "scalar_perturbations": ["perturb", "Phi_g", "Psi_g", "Phi_f", "Psi_f", "Newtonian gauge", "Bardeen", "scalar mode"],
    "weyl_response": ["Weyl", "Phi+Psi", "Phi + Psi", "Sigma", "lensing", "E_G", "EG", "r_Wdelta", "Wdelta"],
    "projector_language": ["projector", "projection", "eigenmode", "eigenvector", "adiabatic", "relative mode", "orthogonal", "inner product"],
    "abg_dictionary": ["a_bg", "a_{\\rm bg}", "DeltaOmega", "Delta\\Omega", "s_surv", "s_{\\rm surv}", "y_t"],
}


def safe_read(path: Path, max_bytes: int = 250_000) -> str:
    try:
        data = path.read_bytes()[:max_bytes]
        return data.decode("utf-8", errors="replace")
    except Exception:
        return ""


def should_skip(path: Path) -> bool:
    s = "/" + str(path).replace("\\", "/")
    return any(pat in s for pat in SKIP_PATTERNS)


def load_prior_values() -> Dict[str, float]:
    vals = {
        "a_bg": 0.403687948,
        "sqrt_a_bg": math.sqrt(0.403687948),
        "r_free": 0.630864746,
        "r2": 0.397990327746,
        "delta_chi2": 0.008545,
        "theta_bg_deg": math.degrees(math.acos(math.sqrt(0.403687948))),
    }
    for fn in [
        OUTDIR / "paperVI_EG_T9A_theory_appendix.json",
        OUTDIR / "paperVI_EG_T8_auto_cross_amplitude_space.json",
        OUTDIR / "paperVI_EG_T3_A4_overlap_reduction.json",
        OUTDIR / "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json",
    ]:
        if not fn.exists():
            continue
        try:
            data = json.loads(fn.read_text(encoding="utf-8"))
        except Exception:
            continue
        text = json.dumps(data)
        def find_number(keys):
            for k in keys:
                m = re.search(r'"' + re.escape(k) + r'"\s*:\s*([-+0-9.eE]+)', text)
                if m:
                    try:
                        return float(m.group(1))
                    except Exception:
                        pass
            return None
        v = find_number(["a_bg"])
        if v is not None:
            vals["a_bg"] = v
            vals["sqrt_a_bg"] = math.sqrt(v)
            vals["theta_bg_deg"] = math.degrees(math.acos(math.sqrt(v)))
        v = find_number(["r_wdelta_free", "r_free", "T1_r_free"])
        if v is not None:
            vals["r_free"] = v
            vals["r2"] = v * v
        v = find_number(["delta_chi2_fixed_sqrt_abg_minus_free", "delta_chi2"])
        if v is not None:
            vals["delta_chi2"] = v
    return vals


def symbolic_requirements(a_bg: float) -> Dict[str, float]:
    sqrt_a = math.sqrt(a_bg)
    sqrt_oma = math.sqrt(1.0 - a_bg)
    theta = math.acos(sqrt_a)
    tan_theta = math.tan(theta)
    ratio_private_to_matter_amp = sqrt_oma / sqrt_a
    ratio_private_to_matter_power = (1.0 - a_bg) / a_bg
    # 2x2 symmetric toy matrix where eigenvector angle theta obeys tan(2theta)=2b/(a-d).
    tan_2theta = math.tan(2.0 * theta)
    # If choose (a-d)=1 normalization, b required is 0.5*tan(2theta).
    b_over_delta = 0.5 * tan_2theta
    return {
        "sqrt_a_bg": sqrt_a,
        "sqrt_one_minus_a_bg": sqrt_oma,
        "theta_bg_rad": theta,
        "theta_bg_deg": math.degrees(theta),
        "tan_theta": tan_theta,
        "private_to_matter_amplitude_ratio": ratio_private_to_matter_amp,
        "private_to_matter_power_ratio": ratio_private_to_matter_power,
        "tan_2theta": tan_2theta,
        "offdiag_over_diag_split_for_2x2_symmetric_matrix": b_over_delta,
    }


def inventory_files() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    seen = set()
    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix not in TEXT_EXTS or should_skip(path):
                continue
            rp = str(path)
            if rp in seen:
                continue
            seen.add(rp)
            txt = safe_read(path)
            if not txt:
                continue
            group_hits = {}
            total_hits = 0
            for group, toks in TOKEN_GROUPS.items():
                hits = [tok for tok in toks if tok.lower() in txt.lower()]
                if hits:
                    group_hits[group] = hits
                    total_hits += len(hits)
            if not group_hits:
                continue
            # Score files that contain multiple needed groups higher.
            score = total_hits + 3 * len(group_hits)
            rows.append({
                "path": rp,
                "suffix": path.suffix,
                "score": score,
                "groups_found": ";".join(sorted(group_hits.keys())),
                "tokens_found": ";".join([f"{g}:{'|'.join(v)}" for g, v in sorted(group_hits.items())]),
                "contains_projector_and_abg": bool("projector_language" in group_hits and "abg_dictionary" in group_hits),
                "contains_perturbation_and_weyl": bool("scalar_perturbations" in group_hits and "weyl_response" in group_hits),
            })
    rows.sort(key=lambda r: (-r["score"], r["path"]))
    return rows


def evaluate_action_closure(inventory: List[Dict[str, Any]]) -> Tuple[str, List[Dict[str, str]]]:
    # Conservative closure logic: only close if one file has projector+abg+weyl+perturbation groups.
    strong = []
    for r in inventory:
        groups = set(r["groups_found"].split(";"))
        if {"projector_language", "abg_dictionary", "weyl_response", "scalar_perturbations"}.issubset(groups):
            strong.append(r)

    decisions = []
    if strong:
        # Even if strong files exist, the script cannot parse and validate a full field-equation proof.
        # So require manual validation. Do not close automatically.
        status = "P6_T9B_ACTION_LEVEL_CANDIDATE_FOUND_MANUAL_VALIDATION_REQUIRED"
        decisions.append({
            "item": "candidate_action_level_sources",
            "decision": "FOUND_BUT_NOT_VALIDATED",
            "basis": f"{len(strong)} files contain projector/a_bg/Weyl/perturbation tokens, but this gate does not parse a full proof.",
        })
    else:
        status = "P6_T9B_ACTION_LEVEL_NOT_CLOSED_REQUIRE_EXPLICIT_HR_PROJECTOR"
        decisions.append({
            "item": "candidate_action_level_sources",
            "decision": "NOT_FOUND",
            "basis": "No scanned file simultaneously exposed a usable HR/BPB perturbation projector, Weyl response, and a_bg overlap identity.",
        })

    decisions.extend([
        {
            "item": "algebraic_requirement",
            "decision": "LOCKED",
            "basis": "Full action-level closure must derive cos^2(theta_W)=a_bg or equivalently <P_m W,P_m W>/<W,W>=a_bg without inserting it by hand.",
        },
        {
            "item": "non_circularity_test",
            "decision": "REQUIRED",
            "basis": "A proof is circular if it starts by defining P_m so that the overlap is a_bg.",
        },
        {
            "item": "current_paper_use",
            "decision": "KEEP_T9A_CONDITIONAL_APPENDIX",
            "basis": "Until the HR/BPB projector is derived, Paper VI should use the T9A projected-diagnostic closure wording.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T9C_EXPLICIT_HR_PERTURBATION_MATRIX_OR_MANUAL_SOURCE_REVIEW",
            "basis": "Provide or derive the actual scalar perturbation matrix/eigenbasis; then test the mixing angle against theta_bg.",
        },
    ])
    return status, decisions


def write_dict_csv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    vals = load_prior_values()
    req = symbolic_requirements(vals["a_bg"])
    inv = inventory_files()
    status, decisions = evaluate_action_closure(inv)

    derivation_steps = [
        {
            "step": "T9B.1",
            "statement": "Start from the HR/BPB scalar perturbation sector and identify the Weyl response vector W in the E_G projected space.",
            "needed_to_close": "Explicit perturbation variables and Weyl transfer response.",
            "status": "NOT_CLOSED_BY_THIS_GATE",
        },
        {
            "step": "T9B.2",
            "statement": "Diagonalize or project the perturbation sector into matter-correlated mode e_m and private/relative mode e_x.",
            "needed_to_close": "A non-circular matter-mode projector P_m.",
            "status": "NOT_CLOSED_BY_THIS_GATE",
        },
        {
            "step": "T9B.3",
            "statement": "Show that the Weyl eigenvector/mixing angle obeys cos^2(theta_W)=a_bg.",
            "needed_to_close": "Derive theta_W=acos(sqrt(a_bg)) from field-equation coefficients.",
            "status": "TARGET_IDENTITY",
        },
        {
            "step": "T9B.4",
            "statement": "Then r_Wdelta=sqrt(a_bg) follows from T2/T3 without fitting r in E_G.",
            "needed_to_close": "T2/T3 already provide this algebraic implication.",
            "status": "CONDITIONAL_ALREADY_PROVEN",
        },
    ]

    blockers = [
        {
            "blocker": "missing_explicit_projector",
            "description": "No validated expression for P_m in terms of HR/BPB perturbation variables is available to this gate.",
            "consequence": "Cannot claim action-level derivation of <P_m W,P_m W>/<W,W>=a_bg.",
        },
        {
            "blocker": "mixing_angle_not_from_field_coefficients",
            "description": "theta_bg=acos(sqrt(a_bg)) is numerically known, but not derived here from a perturbation matrix/eigenvector.",
            "consequence": "a_bg remains a projected-diagnostic overlap target, not a proven HR perturbation eigen-overlap.",
        },
        {
            "blocker": "projection_inner_product_dependence",
            "description": "The E_G projected inner product and the raw HR perturbation variables live in different spaces unless a window projection is specified.",
            "consequence": "An action proof must commute or explicitly map HR perturbations through the E_G projection.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "ACTION_LEVEL_ATTEMPT_DONE_NOT_CLOSED_AUTOMATICALLY",
        },
        "prior_values": vals,
        "symbolic_requirements": req,
        "derivation_steps": derivation_steps,
        "blockers": blockers,
        "inventory_counts": {
            "candidate_files": len(inv),
            "top_file": inv[0]["path"] if inv else "",
        },
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T9B_action_level_projector_attempt.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    write_dict_csv(
        OUTDIR / "paperVI_EG_T9B_symbolic_requirements.csv",
        [{"quantity": k, "value": v} for k, v in req.items()],
        ["quantity", "value"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T9B_theory_inventory_candidates.csv",
        inv[:200],
        ["path", "suffix", "score", "groups_found", "tokens_found", "contains_projector_and_abg", "contains_perturbation_and_weyl"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T9B_derivation_steps.csv",
        derivation_steps,
        ["step", "statement", "needed_to_close", "status"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T9B_blockers.csv",
        blockers,
        ["blocker", "description", "consequence"],
    )
    write_dict_csv(
        OUTDIR / "paperVI_EG_T9B_decisions.csv",
        decisions,
        ["item", "decision", "basis"],
    )

    top = inv[:20]
    md = f"""# Paper VI / E_G T9B - action-level projector derivation attempt

## Status

**{status}**

Non-destructive theory gate. No likelihood run. NS untouched.

## Purpose

T9A produced a paper-ready projected-diagnostic closure:

```text
r_Wdelta = sqrt(a_bg)
```

T9B attempts to promote this to an action-level HR/BPB perturbation derivation. The strict target is:

```text
<P_m W,P_m W>_EG / <W,W>_EG = a_bg
```

This gate does not close the action-level proof unless the HR/BPB projector/eigenbasis is explicitly available.

## Locked numerical target

| quantity | value |
|---|---:|
| a_bg | {vals['a_bg']:.12f} |
| sqrt(a_bg) | {vals['sqrt_a_bg']:.12f} |
| theta_bg deg | {vals['theta_bg_deg']:.12f} |
| free r_Wdelta | {vals['r_free']:.12f} |
| delta chi2 projected closure | {vals['delta_chi2']:.12g} |

## Symbolic requirements for a full proof

The action-level derivation must produce the Weyl/matter mixing angle, not assume it.

| requirement | value |
|---|---:|
| sqrt(1-a_bg) | {req['sqrt_one_minus_a_bg']:.12f} |
| private/matter amplitude ratio | {req['private_to_matter_amplitude_ratio']:.12f} |
| private/matter power ratio | {req['private_to_matter_power_ratio']:.12f} |
| tan(theta_bg) | {req['tan_theta']:.12f} |
| tan(2 theta_bg) | {req['tan_2theta']:.12f} |
| 2x2 offdiag/diag-split requirement | {req['offdiag_over_diag_split_for_2x2_symmetric_matrix']:.12f} |

A minimal two-mode perturbation matrix would need an eigenvector whose squared matter overlap is exactly `a_bg`. Equivalently, if `W_hat=cos(theta_W)e_m+sin(theta_W)e_x`, the field equations must imply:

```text
cos^2(theta_W)=a_bg.
```

## Derivation status

| step | status | needed to close |
|---|---|---|
"""
    for row in derivation_steps:
        md += f"| {row['step']} | {row['status']} | {row['needed_to_close']} |\n"

    md += """
## Blockers

| blocker | consequence |
|---|---|
"""
    for row in blockers:
        md += f"| {row['blocker']} | {row['consequence']} |\n"

    md += f"""
## Repository inventory

Candidate theory files scanned: **{len(inv)}**.

Top candidates:

| score | path | groups |
|---:|---|---|
"""
    for row in top:
        md += f"| {row['score']} | `{row['path']}` | {row['groups_found']} |\n"

    md += """
## Decisions

| item | decision | basis |
|---|---|---|
"""
    for row in decisions:
        md += f"| {row['item']} | {row['decision']} | {row['basis']} |\n"

    md += """
## Locked conclusion

T9B does not invalidate the projected closure. It clarifies its exact action-level target.

The T1--T8 chain proves and verifies a projected-diagnostic closure. A complete HR/BPB derivation still requires an explicit perturbation-sector projector or mixing matrix that yields:

```text
cos^2(theta_W)=a_bg
```

without defining the projector to make it true by construction.

Until that object is supplied or derived, Paper VI should keep the T9A wording: strong projected-diagnostic closure, action-level projector derivation open.

## Next gate

P6-T9C should provide one of the following:

1. an explicit HR/BPB scalar perturbation matrix/eigenbasis;
2. an explicit non-circular expression for P_m in terms of perturbation variables;
3. a manual source review of the top inventory candidates if they contain the missing derivation.
"""
    (OUTDIR / "paperVI_EG_T9B_action_level_projector_attempt.md").write_text(md, encoding="utf-8")

    print(status)
    print("Candidate theory files:", len(inv))
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
