#!/usr/bin/env python3
"""
Paper VI / E_G T9E -- source-match density-weighted mixing relation.

Non-destructive gate.
- No likelihood run.
- NS paths are ignored.
- Excludes self-generated Paper VI E_G theory files to avoid circular validation.

Purpose
-------
T9D upgraded the projected closure to a conditional density-weighted HR/BPB
mixing candidate. The remaining D2 source-match target is:

    tan^2(theta_W) = rho_f_eff / rho_g_eff

or equivalently, in a two-mode perturbation/eigenbasis, the Weyl scalar mixing
angle must be fixed by sector weights rather than inserted by hand.

This script searches the repo for non-circular internal support for D2 and
classifies candidates. It does not mark the full action-level proof as closed
unless a candidate contains enough evidence of an explicit HR/BPB scalar
perturbation mixing relation, matrix/eigenvector, and density ratio.
"""
from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

SEARCH_ROOTS = [Path("outputs"), Path("runs"), Path("src"), Path("pipeline"), Path("papers"), Path("docs")]
EXTS = {".py", ".md", ".tex", ".txt", ".json", ".csv", ".yaml", ".yml"}

# Avoid circular validation from the newly generated Paper VI gates and avoid NS.
EXCLUDE_SUBSTRINGS = [
    "outputs/paperVI_EG_theory",
    "runs/paperVI",
    "/ns/",
    "outputs/ns",
    "src/ns",
    "checkpoints_strict_global_gate",
    ".git/",
    ".venv/",
    "venv/",
    "__pycache__",
]

# Strong D2 tokens.
TOKEN_GROUPS: Dict[str, List[str]] = {
    "density_ratio": [
        "rho_f/rho_g", "rho_f_eff/rho_g_eff", "rho_f_eff", "rho_g_eff",
        "rho_f", "rho_g", "density-weighted", "density weighted",
        "sector weight", "sector weights", "effective density", "densities",
    ],
    "mixing_angle": [
        "tan^2", "tan2", "tan(theta", "theta_W", "mixing angle",
        "mixing_angle", "cos^2", "cos2", "acos", "arccos",
    ],
    "two_mode": [
        "two-mode", "two mode", "coupled modes", "eigenmode", "eigenvector",
        "eigenbasis", "mass eigen", "adiabatic", "relative mode", "private mode",
    ],
    "matrix_or_coefficients": [
        "matrix", "Hessian", "mass matrix", "mixing matrix", "offdiag",
        "off-diagonal", "diagonalize", "diagonalise", "eigenvalue", "coefficients",
    ],
    "hr_bpb": [
        "Hassan", "Rosen", "bigravity", "bimetric", "HR", "BPB", "beta_n",
        "y_t", "a_bg", "s_surv", "biface", "branch",
    ],
    "weyl_scalar": [
        "Weyl", "Phi+Psi", "Phi + Psi", "Sigma", "lensing", "scalar perturb",
        "Phi_g", "Psi_g", "Phi_f", "Psi_f", "Bardeen", "Newtonian gauge",
    ],
}

PRIORITY_HINTS = [
    "run_cmb_lensing_kernel_phase_fraction.py",
    "autocross_delta_phi_toy_projector",
    "autocross_delta_phi_toy_local_sensitivity",
    "sigma_window_closure_diagnostic",
    "closure",
    "projector",
    "phase_fraction",
]


@dataclass
class Candidate:
    path: str
    suffix: str
    size_bytes: int
    score: int
    classification: str
    groups_found: str
    tokens_found: str
    priority_path: bool
    has_density_ratio: bool
    has_mixing_angle: bool
    has_matrix_or_coefficients: bool
    has_hr_bpb: bool
    has_weyl_scalar: bool
    has_non_circular_potential: bool


def should_skip(path: Path) -> bool:
    s = str(path).replace("\\", "/")
    sl = s.lower()
    for pat in EXCLUDE_SUBSTRINGS:
        if pat.lower() in sl:
            return True
    return False


def read_text_safe(path: Path, max_bytes: int = 1_000_000) -> str:
    try:
        data = path.read_bytes()[:max_bytes]
        return data.decode("utf-8", errors="replace")
    except Exception:
        return ""


def find_token_hits(text: str) -> Tuple[Dict[str, List[str]], int]:
    hits: Dict[str, List[str]] = {}
    score = 0
    lower = text.lower()
    for group, toks in TOKEN_GROUPS.items():
        ghits = []
        for tok in toks:
            if tok.lower() in lower:
                ghits.append(tok)
        if ghits:
            hits[group] = ghits
            # Weight exact critical groups higher.
            if group in {"density_ratio", "mixing_angle"}:
                score += 5 * len(ghits)
            elif group in {"matrix_or_coefficients", "two_mode"}:
                score += 4 * len(ghits)
            else:
                score += 3 * len(ghits)
    return hits, score


def classify(path: Path, hits: Dict[str, List[str]], text: str) -> str:
    groups = set(hits)
    sl = str(path).lower()
    # Explicit pass candidate requires density+mixing+matrix/eigen + HR/BPB + Weyl/scalar.
    if {"density_ratio", "mixing_angle", "matrix_or_coefficients", "hr_bpb", "weyl_scalar"}.issubset(groups):
        return "STRONG_D2_SOURCE_CANDIDATE_MANUAL_VALIDATE"
    if {"density_ratio", "mixing_angle", "two_mode", "hr_bpb"}.issubset(groups):
        return "DENSITY_MIXING_CANDIDATE_MANUAL_VALIDATE"
    if "run_cmb_lensing_kernel_phase_fraction.py" in sl:
        return "KERNEL_PHASE_FRACTION_REVIEW"
    if "autocross_delta_phi_toy_projector" in sl:
        return "TOY_PROJECTOR_REVIEW"
    if "autocross_delta_phi_toy_local_sensitivity" in sl:
        return "TOY_LOCAL_SENSITIVITY_REVIEW"
    if "sigma_window_closure_diagnostic" in sl:
        return "SIGMA_WINDOW_CONTEXT_REVIEW"
    if "density_ratio" in groups or "mixing_angle" in groups:
        return "PARTIAL_D2_CONTEXT"
    return "LOW_VALUE_OR_CONTEXT"


def iter_files() -> Iterable[Path]:
    seen = set()
    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file() or p.suffix not in EXTS:
                continue
            if should_skip(p):
                continue
            rp = str(p)
            if rp in seen:
                continue
            seen.add(rp)
            yield p


def collect_candidates() -> List[Candidate]:
    out: List[Candidate] = []
    for p in iter_files():
        text = read_text_safe(p)
        if not text.strip():
            continue
        hits, score = find_token_hits(text)
        priority = any(h.lower() in str(p).lower() for h in PRIORITY_HINTS)
        if score <= 0 and not priority:
            continue
        if priority:
            score += 10
        classification = classify(p, hits, text)
        groups = set(hits)
        tokens_found = []
        for g, toks in hits.items():
            tokens_found.append(f"{g}:{'|'.join(toks)}")
        out.append(Candidate(
            path=str(p),
            suffix=p.suffix,
            size_bytes=p.stat().st_size,
            score=score,
            classification=classification,
            groups_found=";".join(sorted(groups)),
            tokens_found=";".join(tokens_found),
            priority_path=priority,
            has_density_ratio="density_ratio" in groups,
            has_mixing_angle="mixing_angle" in groups,
            has_matrix_or_coefficients="matrix_or_coefficients" in groups,
            has_hr_bpb="hr_bpb" in groups,
            has_weyl_scalar="weyl_scalar" in groups,
            has_non_circular_potential=classification in {
                "STRONG_D2_SOURCE_CANDIDATE_MANUAL_VALIDATE",
                "DENSITY_MIXING_CANDIDATE_MANUAL_VALIDATE",
            },
        ))
    out.sort(key=lambda c: (-c.score, c.path))
    return out


def extract_snippets(path: Path, patterns: List[str], max_snippets: int = 8, radius: int = 260) -> List[Dict[str, str]]:
    text = read_text_safe(path, max_bytes=2_000_000)
    snippets = []
    low = text.lower()
    seen = set()
    for pat in patterns:
        start = 0
        pl = pat.lower()
        while len(snippets) < max_snippets:
            idx = low.find(pl, start)
            if idx < 0:
                break
            a = max(0, idx - radius)
            b = min(len(text), idx + len(pat) + radius)
            key = (a, b)
            if key not in seen:
                snip = text[a:b].replace("\r", " ").replace("\n", " ")
                snippets.append({"path": str(path), "pattern": pat, "snippet": snip})
                seen.add(key)
            start = idx + len(pat)
    return snippets


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def read_t9d_values() -> Dict[str, float]:
    defaults = {
        "a_bg": 0.403687948,
        "sqrt_a_bg": math.sqrt(0.403687948),
        "rho_f_over_rho_g": (1 - 0.403687948) / 0.403687948,
        "theta_deg": math.degrees(math.atan(math.sqrt((1 - 0.403687948) / 0.403687948))),
        "delta_chi2": 0.008545,
    }
    path = OUTDIR / "paperVI_EG_T9D_key_values.csv"
    if not path.exists():
        return defaults
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                q = row.get("quantity") or row.get("check") or ""
                v = row.get("value")
                try:
                    x = float(v)
                except Exception:
                    continue
                if q in {"a_bg", "sqrt_a_bg", "rho_f_over_rho_g", "theta_deg", "delta_chi2_fixed_minus_free"}:
                    key = "delta_chi2" if q == "delta_chi2_fixed_minus_free" else q
                    defaults[key] = x
    except Exception:
        pass
    return defaults


def main() -> None:
    values = read_t9d_values()
    candidates = collect_candidates()
    strong = [c for c in candidates if c.classification == "STRONG_D2_SOURCE_CANDIDATE_MANUAL_VALIDATE"]
    density_mix = [c for c in candidates if c.classification == "DENSITY_MIXING_CANDIDATE_MANUAL_VALIDATE"]
    priority = [c for c in candidates if c.priority_path][:80]

    # Extract snippets from top manually relevant candidates.
    snip_patterns = [
        "tan^2", "tan(theta", "theta_W", "mixing angle", "rho_f", "rho_g",
        "density", "projector", "eigen", "matrix", "Weyl", "Sigma", "a_bg",
    ]
    snippet_paths = []
    for c in candidates:
        if c.classification in {
            "STRONG_D2_SOURCE_CANDIDATE_MANUAL_VALIDATE",
            "DENSITY_MIXING_CANDIDATE_MANUAL_VALIDATE",
            "KERNEL_PHASE_FRACTION_REVIEW",
            "TOY_PROJECTOR_REVIEW",
            "TOY_LOCAL_SENSITIVITY_REVIEW",
            "SIGMA_WINDOW_CONTEXT_REVIEW",
        }:
            snippet_paths.append(Path(c.path))
        if len(snippet_paths) >= 12:
            break
    snippets: List[Dict[str, str]] = []
    for p in snippet_paths:
        snippets.extend(extract_snippets(p, snip_patterns, max_snippets=8))

    if strong:
        status = "P6_T9E_DENSITY_WEIGHTED_MIXING_STRONG_SOURCE_CANDIDATE_FOUND_NOT_CLOSED"
    elif density_mix:
        status = "P6_T9E_DENSITY_WEIGHTED_MIXING_PARTIAL_SOURCE_CANDIDATE_FOUND_NOT_CLOSED"
    else:
        status = "P6_T9E_DENSITY_WEIGHTED_MIXING_NOT_SOURCE_MATCHED_INTERNAL"

    decisions = [
        {
            "item": "self_generated_files",
            "decision": "EXCLUDED",
            "basis": "outputs/paperVI_EG_theory and runs/paperVI excluded to avoid circular validation.",
        },
        {
            "item": "D2_source_match",
            "decision": "STRONG_CANDIDATE_FOUND" if strong else "PARTIAL_CANDIDATE_FOUND" if density_mix else "NOT_SOURCE_MATCHED_INTERNAL",
            "basis": f"strong={len(strong)}, density_mix={len(density_mix)}, total_candidates={len(candidates)}.",
        },
        {
            "item": "full_action_level_proof",
            "decision": "NOT_CLOSED_BY_T9E",
            "basis": "T9E is an internal source-match gate; a strong candidate still requires manual validation of non-circular HR/BPB perturbation equations.",
        },
        {
            "item": "paper_use",
            "decision": "USE_T9D_AS_CONDITIONAL_DENSITY_WEIGHTED_DERIVATION",
            "basis": "T9D algebra and non-circular diagnostic status are valid; D2 remains a labelled source-match requirement.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T9F_MANUAL_VALIDATE_D2_CANDIDATE_OR_INSERT_V02",
            "basis": "If a strong/partial candidate exists, inspect snippets manually; otherwise insert T9D v0.2 wording with D2 guardrail.",
        },
    ]

    result = {
        "status": {
            "gate_status": status,
            "likelihood_run": "NO",
            "ns_touched": "NO",
            "theory_status": "D2_SOURCE_MATCH_SEARCH_DONE_NOT_FULL_ACTION_PROOF",
        },
        "t9d_values": values,
        "counts": {
            "total_candidates": len(candidates),
            "strong_d2_candidates": len(strong),
            "density_mixing_candidates": len(density_mix),
            "priority_candidates": len([c for c in candidates if c.priority_path]),
            "snippets_extracted": len(snippets),
        },
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T9E_source_match_density_weighted_mixing.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    cand_rows = [asdict(c) for c in candidates]
    write_csv(
        OUTDIR / "paperVI_EG_T9E_source_match_candidates.csv",
        cand_rows,
        [
            "path", "suffix", "size_bytes", "score", "classification", "groups_found", "tokens_found",
            "priority_path", "has_density_ratio", "has_mixing_angle", "has_matrix_or_coefficients",
            "has_hr_bpb", "has_weyl_scalar", "has_non_circular_potential",
        ],
    )
    write_csv(
        OUTDIR / "paperVI_EG_T9E_priority_candidates.csv",
        [asdict(c) for c in priority],
        [
            "path", "suffix", "size_bytes", "score", "classification", "groups_found", "tokens_found",
            "priority_path", "has_density_ratio", "has_mixing_angle", "has_matrix_or_coefficients",
            "has_hr_bpb", "has_weyl_scalar", "has_non_circular_potential",
        ],
    )
    write_csv(
        OUTDIR / "paperVI_EG_T9E_snippets.csv",
        snippets,
        ["path", "pattern", "snippet"],
    )
    write_csv(
        OUTDIR / "paperVI_EG_T9E_decisions.csv",
        decisions,
        ["item", "decision", "basis"],
    )

    md = f"""# Paper VI / E_G T9E - source-match density-weighted mixing relation

## Status

**{status}**

Non-destructive source-match gate. No likelihood run. NS untouched.

## Purpose

T9D established the conditional density-weighted bridge:

```text
tan^2(theta_W)=rho_f_eff/rho_g_eff
rho_f_eff/rho_g_eff=(1-a_bg)/a_bg
=> cos^2(theta_W)=a_bg
=> r_Wdelta=sqrt(a_bg)
```

T9E searches the non-generated repository sources for internal support for the remaining D2 relation:

```text
tan^2(theta_W)=rho_f_eff/rho_g_eff
```

The gate excludes the self-generated Paper VI theory outputs and `runs/paperVI` to avoid circular validation.

## Locked target values from T9D

| quantity | value |
|---|---:|
| a_bg | {values['a_bg']:.12f} |
| sqrt(a_bg) | {values['sqrt_a_bg']:.12f} |
| rho_f/rho_g | {values['rho_f_over_rho_g']:.12f} |
| theta_W deg | {values['theta_deg']:.12f} |
| delta chi2 fixed-free | {values['delta_chi2']:.12f} |

## Candidate counts

| item | count |
|---|---:|
| total candidates | {len(candidates)} |
| strong D2 candidates | {len(strong)} |
| density-mixing candidates | {len(density_mix)} |
| priority-path candidates | {len([c for c in candidates if c.priority_path])} |
| snippets extracted | {len(snippets)} |

## Top candidates

| score | classification | path | groups |
|---:|---|---|---|
"""
    for c in candidates[:30]:
        md += f"| {c.score} | {c.classification} | `{c.path}` | {c.groups_found} |\n"

    md += """
## Interpretation

T9E is a source-match gate, not a proof gate. A file can be a strong candidate only if it contains the ingredients needed to validate D2 without defining the projector to force the answer:

1. sector density/weight ratio,
2. scalar or Weyl two-mode mixing angle,
3. matrix/eigenmode or coefficient-level origin,
4. HR/BPB or bigravity context,
5. Weyl/lensing scalar response.

## Decisions

| item | decision | basis |
|---|---|---|
"""
    for d in decisions:
        md += f"| {d['item']} | {d['decision']} | {d['basis']} |\n"

    md += """
## Locked conclusion

T9E does not by itself close the full action-level HR/BPB proof. It either identifies internal candidates for manual D2 validation or confirms that the T9D relation should remain a conditional density-weighted mixing assumption.

For Paper VI v0.2, the safe wording remains: density-weighted HR/BPB scalar Weyl mixing gives the closure under the explicit D2 condition, and the E_G auto/cross diagnostic verifies the resulting parameter-free prediction.

## Next gate

P6-T9F should either manually validate the best D2 source candidate(s), or stop the action-level branch and insert the T9D v0.2 conditional wording.
"""
    (OUTDIR / "paperVI_EG_T9E_source_match_density_weighted_mixing.md").write_text(md, encoding="utf-8")

    print(status)
    print("Candidates:", len(candidates))
    print("Strong D2 candidates:", len(strong))
    print("Density-mixing candidates:", len(density_mix))
    print("Snippets:", len(snippets))
    print("Wrote outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
