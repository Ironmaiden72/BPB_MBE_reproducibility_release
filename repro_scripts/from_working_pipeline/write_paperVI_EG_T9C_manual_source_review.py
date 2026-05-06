#!/usr/bin/env python3
"""
Paper VI / E_G T9C -- manual source review for an action-level HR/BPB projector.

Non-destructive gate.
- Reads candidate files from the repo.
- Excludes newly generated T-series artefacts so they do not self-confirm.
- Does not touch NS paths.
- Runs no likelihood.
- Writes outputs/paperVI_EG_theory only.

Purpose
-------
T9B found candidate files with projector/a_bg/Weyl/perturbation tokens but did not
validate a proof. T9C performs a stricter source review for older/non-generated
candidate files and asks whether any file contains a non-circular derivation of

    <P_m W, P_m W>_EG / <W,W>_EG = a_bg

or equivalently

    cos^2(theta_W) = a_bg.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any

OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Strongly exclude generated theory-chain files to avoid circular source review.
EXCLUDE_SUBSTRINGS = [
    "/outputs/paperVI_EG_theory/",
    "outputs/paperVI_EG_theory/",
    "/runs/paperVI/",
    "runs/paperVI/",
    "/outputs/ns/",
    "outputs/ns/",
    "/src/ns/",
    "src/ns/",
    "/runs/ns/",
    "runs/ns/",
    "checkpoints_strict_global_gate",
    "__pycache__",
    ".git/",
    ".venv/",
    "/venv/",
]

# Hand-prioritized paths from T9B and the running discussion.
PRIORITY_PATTERNS = [
    "outputs/cosmology/autocross_delta_phi_toy_projector",
    "outputs/cosmology/autocross_delta_phi_toy_local_sensitivity",
    "runs/cosmology/run_cmb_lensing_kernel_phase_fraction.py",
    "outputs/cosmology/sigma_window_closure_diagnostic",
    "outputs/cosmology/act_eg_autocross",
    "outputs/cosmology/eg_projected_auto_cross",
    "outputs/cosmology/physical_autocross",
    "outputs/cosmology/autocross_csigma",
    "outputs/cmb/cmb_lensing_weyl_proxy",
    "runs/cosmology/",
]

EXTS = {".py", ".md", ".json", ".csv", ".txt", ".tex"}

TOKEN_GROUPS = {
    "target_identity": [
        "<P_m W", "P_m W", "Pm W", "P_m", "overlap", "fraction", "cos^2", "theta_W", "sqrt(a_bg)",
        "r_Wdelta", "Wdelta", "r_global", "XX/AA", "X/X", "X_i", "A_i",
    ],
    "abg_dictionary": ["a_bg", "abg", "DeltaOmega", "Delta_Omega", "y_t", "s_surv"],
    "projector_language": ["projector", "projection", "eigenmode", "eigenvector", "basis", "orthogonal", "inner product", "overlap"],
    "perturbation_language": ["perturb", "delta_phi", "Phi_g", "Psi_g", "Phi_f", "Psi_f", "Bardeen", "Newtonian", "scalar"],
    "weyl_language": ["Weyl", "Phi+Psi", "Phi + Psi", "Sigma", "lensing", "kappa", "E_G", "EG", "auto", "cross"],
    "action_language": ["Hassan", "Rosen", "bigravity", "bimetric", "HR", "action", "beta_n", "potential"],
    "matrix_language": ["matrix", "Jacobian", "Hessian", "eigen", "diagonal", "mixing", "angle", "2x2", "offdiag"],
}

CLOSURE_TARGETS = {
    "a_bg": 0.403687948,
    "sqrt_a_bg": math.sqrt(0.403687948),
    "theta_bg_deg": math.degrees(math.acos(math.sqrt(0.403687948))),
    "tan_theta": math.sqrt(1.0 - 0.403687948) / math.sqrt(0.403687948),
}


def norm_path(p: Path) -> str:
    return str(p).replace("\\", "/")


def should_exclude(p: Path) -> bool:
    s = norm_path(p)
    for bad in EXCLUDE_SUBSTRINGS:
        if bad in s:
            return True
    return False


def safe_read_text(path: Path, max_chars: int = 2_000_000) -> str:
    try:
        txt = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        try:
            txt = path.read_text(encoding="latin-1", errors="replace")
        except Exception as e:
            return f"__READ_ERROR__ {repr(e)}"
    if len(txt) > max_chars:
        return txt[:max_chars] + "\n__TRUNCATED_FOR_REVIEW__"
    return txt


def score_text(text: str) -> Tuple[int, Dict[str, int], Dict[str, List[str]]]:
    lower = text.lower()
    group_counts: Dict[str, int] = {}
    group_hits: Dict[str, List[str]] = {}
    score = 0
    for group, toks in TOKEN_GROUPS.items():
        hits = []
        for tok in toks:
            if tok.lower() in lower:
                hits.append(tok)
        group_hits[group] = hits
        group_counts[group] = len(hits)
        # Different groups carry different proof value.
        weight = {
            "target_identity": 4,
            "matrix_language": 3,
            "projector_language": 3,
            "perturbation_language": 3,
            "weyl_language": 2,
            "abg_dictionary": 3,
            "action_language": 2,
        }.get(group, 1)
        score += weight * len(hits)
    return score, group_counts, group_hits


def classify_candidate(path: Path, text: str, group_counts: Dict[str, int]) -> Tuple[str, str]:
    s = norm_path(path).lower()
    lower = text.lower()

    has_target = group_counts.get("target_identity", 0) >= 2
    has_abg = group_counts.get("abg_dictionary", 0) >= 1
    has_proj = group_counts.get("projector_language", 0) >= 2
    has_pert = group_counts.get("perturbation_language", 0) >= 2
    has_weyl = group_counts.get("weyl_language", 0) >= 2
    has_matrix = group_counts.get("matrix_language", 0) >= 2
    has_action = group_counts.get("action_language", 0) >= 1

    # Direct explicit proof language is very rare; do not overclaim.
    explicit_identity = (
        "<p_m w" in lower or "pm w" in lower or "cos^2(theta_w)" in lower or
        "cos^2(theta" in lower or "theta_w" in lower
    )
    non_circular_hint = (
        "derive" in lower or "from field" in lower or "field equation" in lower or
        "eigenvector" in lower or "eigenbasis" in lower or "diagonal" in lower
    )

    if explicit_identity and has_abg and has_matrix and has_pert and has_weyl:
        return "HIGH_VALUE_MANUAL_REVIEW", "Contains target/matrix/perturbation/Weyl tokens; inspect manually for non-circular derivation."

    if "autocross_delta_phi_toy_projector" in s:
        return "TOY_PROJECTOR_CANDIDATE", "Older toy projector candidate; likely diagnostic/scaffold, must check whether P_m is non-circular."

    if "delta_phi" in s and (has_proj or has_matrix) and has_weyl:
        return "DELTA_PHI_PROJECTOR_REVIEW", "delta_phi/Weyl/projector language present; candidate for manual review."

    if "run_cmb_lensing_kernel_phase_fraction.py" in s or "kernel_phase_fraction" in s:
        return "KERNEL_PHASE_FRACTION_REVIEW", "May contain phase/window fraction diagnostic; inspect if it defines a true projector or only a kernel split."

    if has_target and has_abg and has_proj and (has_pert or has_weyl):
        return "MEDIUM_VALUE_REVIEW", "Contains several target tokens but not enough evidence for explicit action-level proof."

    if has_weyl and (has_proj or has_matrix):
        return "WEYL_DIAGNOSTIC_REVIEW", "Weyl/projector-like diagnostic, likely not action-level."

    return "LOW_VALUE_OR_CONTEXT", "Tokens present but no clear action-level projector structure."


def extract_snippets(path: Path, text: str, max_snippets: int = 8, radius: int = 260) -> List[Dict[str, Any]]:
    terms = [
        "a_bg", "r_Wdelta", "Wdelta", "projector", "projection", "delta_phi",
        "Weyl", "Sigma", "eigen", "mixing", "theta", "overlap", "Phi_g", "Psi_g", "Phi_f", "Psi_f",
    ]
    snippets = []
    used_spans = []
    lower = text.lower()
    for term in terms:
        idx = lower.find(term.lower())
        while idx != -1 and len(snippets) < max_snippets:
            start = max(0, idx - radius)
            end = min(len(text), idx + radius)
            if any(abs(start - u[0]) < radius for u in used_spans):
                idx = lower.find(term.lower(), idx + len(term))
                continue
            snip = text[start:end].replace("\r", " ").replace("\n", " ")
            snip = re.sub(r"\s+", " ", snip).strip()
            snippets.append({
                "path": norm_path(path),
                "term": term,
                "start": start,
                "snippet": snip,
            })
            used_spans.append((start, end))
            idx = lower.find(term.lower(), idx + len(term))
        if len(snippets) >= max_snippets:
            break
    return snippets


def discover_files() -> List[Path]:
    roots = [Path("outputs/cosmology"), Path("outputs/cmb"), Path("runs/cosmology"), Path("pipeline"), Path("src")]
    files: List[Path] = []
    seen = set()
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if should_exclude(p):
                continue
            if p.suffix not in EXTS:
                continue
            s = norm_path(p)
            if s in seen:
                continue
            seen.add(s)
            # Keep files with some obvious relevance in path first.
            path_lower = s.lower()
            if any(tok in path_lower for tok in ["autocross", "delta_phi", "projector", "weyl", "sigma", "lensing", "eg", "cmb", "kernel", "phase", "bigravity", "hr"]):
                files.append(p)
    return files


def write_csv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    files = discover_files()

    reviewed = []
    snippets_all = []

    for p in files:
        text = safe_read_text(p)
        if text.startswith("__READ_ERROR__"):
            continue
        score, group_counts, group_hits = score_text(text + "\n" + norm_path(p))
        if score <= 0:
            continue
        classification, reason = classify_candidate(p, text, group_counts)
        priority = any(pat.lower() in norm_path(p).lower() for pat in PRIORITY_PATTERNS)
        reviewed.append({
            "path": norm_path(p),
            "suffix": p.suffix,
            "size_bytes": p.stat().st_size if p.exists() else -1,
            "score": score,
            "priority_path": priority,
            "classification": classification,
            "reason": reason,
            "target_identity_hits": group_counts.get("target_identity", 0),
            "abg_hits": group_counts.get("abg_dictionary", 0),
            "projector_hits": group_counts.get("projector_language", 0),
            "perturbation_hits": group_counts.get("perturbation_language", 0),
            "weyl_hits": group_counts.get("weyl_language", 0),
            "action_hits": group_counts.get("action_language", 0),
            "matrix_hits": group_counts.get("matrix_language", 0),
            "tokens_found": ";".join(f"{g}:{'|'.join(h)}" for g, h in group_hits.items() if h),
        })
        if classification not in {"LOW_VALUE_OR_CONTEXT"} or priority:
            snippets_all.extend(extract_snippets(p, text))

    reviewed.sort(key=lambda r: (not r["priority_path"], -r["score"], r["path"]))

    high = [r for r in reviewed if r["classification"] in {"HIGH_VALUE_MANUAL_REVIEW", "TOY_PROJECTOR_CANDIDATE", "DELTA_PHI_PROJECTOR_REVIEW", "KERNEL_PHASE_FRACTION_REVIEW"}]
    medium = [r for r in reviewed if r["classification"] == "MEDIUM_VALUE_REVIEW"]

    # Conservative final status: this gate does not manually validate math, it prepares review.
    if high:
        status = "P6_T9C_MANUAL_SOURCE_REVIEW_CANDIDATES_PRIORITIZED_ACTION_NOT_CLOSED"
    elif medium:
        status = "P6_T9C_MANUAL_SOURCE_REVIEW_MEDIUM_CANDIDATES_ACTION_NOT_CLOSED"
    else:
        status = "P6_T9C_NO_ACTION_LEVEL_PROJECTOR_SOURCE_FOUND"

    decisions = [
        {
            "item": "self_generated_files",
            "decision": "EXCLUDED",
            "basis": "outputs/paperVI_EG_theory and runs/paperVI are excluded to prevent circular source validation.",
        },
        {
            "item": "candidate_sources",
            "decision": "PRIORITIZED" if high or medium else "NOT_FOUND",
            "basis": f"Reviewed {len(reviewed)} relevant non-NS files; high-priority candidates={len(high)}, medium={len(medium)}.",
        },
        {
            "item": "action_level_derivation",
            "decision": "NOT_CLOSED_BY_T9C",
            "basis": "T9C is a source review/inventory gate; it does not validate a non-circular HR/BPB perturbation projector proof.",
        },
        {
            "item": "paper_use",
            "decision": "KEEP_T9A_CONDITIONAL_APPENDIX",
            "basis": "Until a candidate is manually validated, Paper VI should keep projected-diagnostic closure wording.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T9D_MANUAL_VALIDATE_TOP_TOY_PROJECTOR_OR_STOP",
            "basis": "Manually inspect top older toy_projector/delta_phi/kernel candidates; if no explicit matrix/projector exists, stop action-level branch.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "MANUAL_SOURCE_REVIEW_DONE_ACTION_LEVEL_NOT_CLOSED",
        },
        "targets": CLOSURE_TARGETS,
        "counts": {
            "files_discovered": len(files),
            "files_reviewed_with_tokens": len(reviewed),
            "high_priority_candidates": len(high),
            "medium_candidates": len(medium),
            "snippets": len(snippets_all),
        },
        "decisions": decisions,
        "top_candidates": reviewed[:40],
    }

    (OUTDIR / "paperVI_EG_T9C_manual_source_review.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    write_csv(
        OUTDIR / "paperVI_EG_T9C_reviewed_candidates.csv",
        reviewed,
        [
            "path", "suffix", "size_bytes", "score", "priority_path", "classification", "reason",
            "target_identity_hits", "abg_hits", "projector_hits", "perturbation_hits", "weyl_hits", "action_hits", "matrix_hits", "tokens_found",
        ],
    )
    write_csv(
        OUTDIR / "paperVI_EG_T9C_snippets.csv",
        snippets_all,
        ["path", "term", "start", "snippet"],
    )
    write_csv(
        OUTDIR / "paperVI_EG_T9C_decisions.csv",
        decisions,
        ["item", "decision", "basis"],
    )

    # Short prioritized manual-review list.
    manual_rows = high + medium[:20]
    write_csv(
        OUTDIR / "paperVI_EG_T9C_manual_review_shortlist.csv",
        manual_rows,
        [
            "path", "suffix", "size_bytes", "score", "priority_path", "classification", "reason",
            "target_identity_hits", "abg_hits", "projector_hits", "perturbation_hits", "weyl_hits", "action_hits", "matrix_hits", "tokens_found",
        ],
    )

    # Markdown report.
    md = f"""# Paper VI / E_G T9C - manual source review for action-level projector

## Status

**{status}**

Non-destructive source-review gate. No likelihood run. NS untouched.

## Purpose

T9B found candidate files but did not validate a full action-level derivation.  
T9C excludes the newly generated T-series files and reviews older/source candidates for a non-circular HR/BPB projector or perturbation mixing matrix that could prove:

```text
<P_m W,P_m W>_EG / <W,W>_EG = a_bg
```

or equivalently:

```text
cos^2(theta_W)=a_bg.
```

## Locked target values

| quantity | value |
|---|---:|
| a_bg | {CLOSURE_TARGETS['a_bg']:.12f} |
| sqrt(a_bg) | {CLOSURE_TARGETS['sqrt_a_bg']:.12f} |
| theta_bg deg | {CLOSURE_TARGETS['theta_bg_deg']:.12f} |
| tan(theta_bg) | {CLOSURE_TARGETS['tan_theta']:.12f} |

## Counts

| item | count |
|---|---:|
| files discovered | {len(files)} |
| reviewed files with relevant tokens | {len(reviewed)} |
| high-priority candidates | {len(high)} |
| medium candidates | {len(medium)} |
| snippets extracted | {len(snippets_all)} |

## Top candidates

| score | classification | path | reason |
|---:|---|---|---|
"""
    for row in reviewed[:30]:
        md += f"| {row['score']} | {row['classification']} | `{row['path']}` | {row['reason']} |\n"

    md += """
## Manual-review shortlist

The following file classes are the most relevant older candidates:

1. `outputs/cosmology/autocross_delta_phi_toy_projector_*`
2. `outputs/cosmology/autocross_delta_phi_toy_local_sensitivity_*`
3. `runs/cosmology/run_cmb_lensing_kernel_phase_fraction.py`
4. auto/cross amplitude diagnostics, only as projected-diagnostic support

## Non-circularity criterion

A candidate is not sufficient if it merely defines a projector so that the answer is `a_bg`.  
To close the action-level proof, it must derive the mixing angle or squared overlap from field-equation coefficients, a perturbation matrix, or a non-circular HR/BPB eigenbasis:

```text
theta_W = acos(sqrt(a_bg))
```

or:

```text
<P_m W,P_m W>_EG / <W,W>_EG = a_bg.
```

## Decisions

| item | decision | basis |
|---|---|---|
"""
    for row in decisions:
        md += f"| {row['item']} | {row['decision']} | {row['basis']} |\n"

    md += """
## Locked conclusion

T9C does not close the action-level derivation. It produces a filtered manual-review shortlist and keeps the non-circular proof criterion strict.

Until a candidate file is manually validated as containing an explicit HR/BPB perturbation projector or mixing matrix, Paper VI should keep the T9A wording: strong projected-diagnostic closure, action-level projector derivation open.

## Next gate

P6-T9D should manually validate the top older candidates, especially the `autocross_delta_phi_toy_projector` files. If those are only toy/projected ansatz diagnostics, the action-level branch should stop here for Paper VI.
"""

    (OUTDIR / "paperVI_EG_T9C_manual_source_review.md").write_text(md, encoding="utf-8")

    print(status)
    print("Reviewed files with tokens:", len(reviewed))
    print("High-priority candidates:", len(high))
    print("Medium candidates:", len(medium))
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
