#!/usr/bin/env python3
"""
Paper VI / E_G theory gate T1 — r_Wdelta ≈ sqrt(a_bg)

Non-destructive audit/closure scaffold.
- Reads no NS files.
- Runs no likelihood.
- Writes only outputs/paperVI_EG_theory/ by default.

Purpose
-------
Lock the numerical status and minimal conditional theory statement for the
candidate closure

    r_Wδ^2 = a_bg,   equivalently   r_Wδ = sqrt(a_bg),

where r_Wδ is the effective projected Weyl--matter correlation coefficient
from the E_G auto/cross diagnostic, and a_bg is the BPB transition parameter
from the low-z/Paper-0 dictionary.

The script uses the already audited locked diagnostic numbers unless an
optional JSON source is supplied later.
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Any, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class LockedInputs:
    # Exact/provisional values from the BPB/MBE diagnostic lock.
    # r_free and a_bg correspond to the provisional seed01/MAP1782 E_G projected diagnostic.
    r_wdelta_free: float = 0.630864746
    a_bg: float = 0.403687948
    chi2_free_r: float = 1.358943
    chi2_fixed_sqrt_abg: float = 1.367488
    chi2_r_equals_abg: float = 23.141918
    chi2_r_equals_sqrt_1_minus_abg: float = 9.791644
    chi2_r_equals_1_minus_abg: float = 1.862853
    n_projected_points: int = 6
    note: str = (
        "Locked diagnostic values: provisional seed01/MAP1782 projected E_G auto/cross global-r audit. "
        "Treat as theory-closure diagnostic, not a final likelihood or full covariance result."
    )


def fmt(x: float, n: int = 9) -> str:
    return f"{x:.{n}f}"


def compute(inp: LockedInputs) -> Dict[str, Any]:
    r = inp.r_wdelta_free
    a = inp.a_bg
    sqrt_a = math.sqrt(a)
    one_minus_a = 1.0 - a
    sqrt_oma = math.sqrt(one_minus_a)
    r2 = r * r

    out = {
        "inputs": asdict(inp),
        "derived": {
            "sqrt_a_bg": sqrt_a,
            "one_minus_a_bg": one_minus_a,
            "sqrt_one_minus_a_bg": sqrt_oma,
            "r_wdelta_free_squared": r2,
            "delta_r_free_minus_sqrt_a_bg": r - sqrt_a,
            "abs_delta_r": abs(r - sqrt_a),
            "relative_delta_r_vs_sqrt_a_bg": (r - sqrt_a) / sqrt_a,
            "delta_r2_minus_a_bg": r2 - a,
            "abs_delta_r2": abs(r2 - a),
            "relative_delta_r2_vs_a_bg": (r2 - a) / a,
            "delta_chi2_fixed_sqrt_abg_minus_free": inp.chi2_fixed_sqrt_abg - inp.chi2_free_r,
            "reduced_chi2_free_r": inp.chi2_free_r / inp.n_projected_points,
            "reduced_chi2_fixed_sqrt_abg": inp.chi2_fixed_sqrt_abg / inp.n_projected_points,
        },
        "status": {},
        "decisions": [],
    }

    dchi = out["derived"]["delta_chi2_fixed_sqrt_abg_minus_free"]
    rel = abs(out["derived"]["relative_delta_r_vs_sqrt_a_bg"])
    if dchi < 0.01 and rel < 0.01:
        status = "P6_T1_RWDELTA_SQRT_ABG_NUMERIC_PASS_STRONG"
    elif dchi < 0.1 and rel < 0.03:
        status = "P6_T1_RWDELTA_SQRT_ABG_NUMERIC_PASS"
    else:
        status = "P6_T1_RWDELTA_SQRT_ABG_NUMERIC_WEAK_OR_FAIL"

    out["status"] = {
        "gate_status": status,
        "numeric_relation": "r_Wdelta ≈ sqrt(a_bg)",
        "theory_status": "CONDITIONAL_CLOSURE_CANDIDATE_NOT_FULL_DERIVATION",
        "likelihood_run": "NO",
        "ns_touched": "NO",
    }
    out["decisions"] = [
        {
            "item": "numeric_closeness",
            "decision": "PASS_STRONG" if status.endswith("PASS_STRONG") else "CHECK",
            "basis": "|r-sqrt(a_bg)| < 0.01 and Δχ² < 0.01 versus free global r",
        },
        {
            "item": "paper_use",
            "decision": "USE_AS_THEORY_CLOSURE_CANDIDATE",
            "basis": "present as conditional theorem/scaffold, not yet as derived action-level equality",
        },
        {
            "item": "claim_strength",
            "decision": "UPGRADE_FROM_EFFECTIVE_R_TO_DICTIONARY_TIED_CANDIDATE",
            "basis": "r^2 is within ~1.4% of a_bg and fixed sqrt(a_bg) is statistically indistinguishable from free r in the compressed diagnostic",
        },
        {
            "item": "next_gate",
            "decision": "DERIVE_OR_DISPROVE_FROM_HR_BPB_WEYL_DECOMPOSITION",
            "basis": "need action-level or perturbation-level identification of a_bg as correlated Weyl-power fraction",
        },
    ]
    return out


def write_csvs(outdir: Path, result: Dict[str, Any]) -> None:
    rows = []
    inputs = result["inputs"]
    der = result["derived"]
    rows.extend([
        ("r_Wdelta_free", inputs["r_wdelta_free"]),
        ("a_bg", inputs["a_bg"]),
        ("sqrt_a_bg", der["sqrt_a_bg"]),
        ("r_Wdelta_free_squared", der["r_wdelta_free_squared"]),
        ("delta_r_free_minus_sqrt_a_bg", der["delta_r_free_minus_sqrt_a_bg"]),
        ("relative_delta_r_vs_sqrt_a_bg", der["relative_delta_r_vs_sqrt_a_bg"]),
        ("delta_r2_minus_a_bg", der["delta_r2_minus_a_bg"]),
        ("relative_delta_r2_vs_a_bg", der["relative_delta_r2_vs_a_bg"]),
        ("chi2_free_r", inputs["chi2_free_r"]),
        ("chi2_fixed_sqrt_abg", inputs["chi2_fixed_sqrt_abg"]),
        ("delta_chi2_fixed_sqrt_abg_minus_free", der["delta_chi2_fixed_sqrt_abg_minus_free"]),
        ("chi2_r_equals_abg", inputs["chi2_r_equals_abg"]),
        ("chi2_r_equals_sqrt_1_minus_abg", inputs["chi2_r_equals_sqrt_1_minus_abg"]),
        ("chi2_r_equals_1_minus_abg", inputs["chi2_r_equals_1_minus_abg"]),
    ])
    with (outdir / "paperVI_EG_T1_rWdelta_sqrt_abg_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        w.writerows(rows)

    comp_rows = [
        {"hypothesis": "free global r", "r_value": inputs["r_wdelta_free"], "chi2": inputs["chi2_free_r"], "delta_chi2_vs_free": 0.0, "verdict": "reference"},
        {"hypothesis": "r = sqrt(a_bg)", "r_value": der["sqrt_a_bg"], "chi2": inputs["chi2_fixed_sqrt_abg"], "delta_chi2_vs_free": der["delta_chi2_fixed_sqrt_abg_minus_free"], "verdict": "statistically indistinguishable in compressed diagnostic"},
        {"hypothesis": "r = 1 - a_bg", "r_value": der["one_minus_a_bg"], "chi2": inputs["chi2_r_equals_1_minus_abg"], "delta_chi2_vs_free": inputs["chi2_r_equals_1_minus_abg"] - inputs["chi2_free_r"], "verdict": "worse but not catastrophic"},
        {"hypothesis": "r = sqrt(1 - a_bg)", "r_value": der["sqrt_one_minus_a_bg"], "chi2": inputs["chi2_r_equals_sqrt_1_minus_abg"], "delta_chi2_vs_free": inputs["chi2_r_equals_sqrt_1_minus_abg"] - inputs["chi2_free_r"], "verdict": "rejected by diagnostic"},
        {"hypothesis": "r = a_bg", "r_value": inputs["a_bg"], "chi2": inputs["chi2_r_equals_abg"], "delta_chi2_vs_free": inputs["chi2_r_equals_abg"] - inputs["chi2_free_r"], "verdict": "strongly rejected by diagnostic"},
    ]
    with (outdir / "paperVI_EG_T1_r_hypothesis_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["hypothesis", "r_value", "chi2", "delta_chi2_vs_free", "verdict"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(comp_rows)

    with (outdir / "paperVI_EG_T1_decisions.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["item", "decision", "basis"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(result["decisions"])


def write_md(outdir: Path, result: Dict[str, Any]) -> None:
    inp = result["inputs"]
    der = result["derived"]
    status = result["status"]["gate_status"]
    md = f"""# Paper VI / E_G T1 — r_Wδ ≈ sqrt(a_bg) theory-closure gate

## Status

**{status}**

This is a non-destructive theory/audit gate. It runs no ACT/Planck/E_G likelihood and touches no NS files.

## Locked numerical relation

| quantity | value |
|---|---:|
| free projected r_Wδ | {inp['r_wdelta_free']:.12f} |
| a_bg | {inp['a_bg']:.12f} |
| sqrt(a_bg) | {der['sqrt_a_bg']:.12f} |
| r_Wδ^2 | {der['r_wdelta_free_squared']:.12f} |
| r_Wδ - sqrt(a_bg) | {der['delta_r_free_minus_sqrt_a_bg']:.12f} |
| relative Δr/sqrt(a_bg) | {100*der['relative_delta_r_vs_sqrt_a_bg']:.6f}% |
| r_Wδ^2 - a_bg | {der['delta_r2_minus_a_bg']:.12f} |
| relative Δr²/a_bg | {100*der['relative_delta_r2_vs_a_bg']:.6f}% |

## χ² comparison

| hypothesis | χ² | Δχ² vs free r |
|---|---:|---:|
| free global r | {inp['chi2_free_r']:.6f} | 0.000000 |
| r = sqrt(a_bg) | {inp['chi2_fixed_sqrt_abg']:.6f} | {der['delta_chi2_fixed_sqrt_abg_minus_free']:.6f} |
| r = 1 - a_bg | {inp['chi2_r_equals_1_minus_abg']:.6f} | {inp['chi2_r_equals_1_minus_abg']-inp['chi2_free_r']:.6f} |
| r = sqrt(1-a_bg) | {inp['chi2_r_equals_sqrt_1_minus_abg']:.6f} | {inp['chi2_r_equals_sqrt_1_minus_abg']-inp['chi2_free_r']:.6f} |
| r = a_bg | {inp['chi2_r_equals_abg']:.6f} | {inp['chi2_r_equals_abg']-inp['chi2_free_r']:.6f} |

The fixed relation `r_Wδ = sqrt(a_bg)` is numerically indistinguishable from the free-r solution in this compressed diagnostic: Δχ² = {der['delta_chi2_fixed_sqrt_abg_minus_free']:.6f}.

## Minimal conditional theory statement

Assume the projected Weyl response can be decomposed into a matter-correlated component and a BPB-private/decorrelated component,

```text
W = sqrt(a_bg) δ_m + sqrt(1-a_bg) ξ_W,
```

with `corr(δ_m, ξ_W)=0` and equal normalized variances in the projected diagnostic window. Then

```text
P_Wδ / sqrt(P_WW P_δδ) = sqrt(a_bg),
```

so the effective cross-correlation coefficient measured by E_G is

```text
r_Wδ = sqrt(a_bg).
```

Under this closure, `a_bg` is not fitted to E_G. It is inherited from the BPB/HR background dictionary and interpreted as the correlated Weyl-power fraction. The residual `1-a_bg` is the decorrelated Weyl auto-power fraction.

## Claim wording candidate

The projected E_G diagnostic does not require a freely fitted Weyl--matter correlation coefficient. The best-fit value `r_Wδ = {inp['r_wdelta_free']:.6f}` is reproduced, to sub-percent accuracy and with only Δχ² = {der['delta_chi2_fixed_sqrt_abg_minus_free']:.6f}, by the dictionary-tied closure `r_Wδ = sqrt(a_bg) = {der['sqrt_a_bg']:.6f}`. This suggests that the BPB transition parameter controls the matter-correlated fraction of the Weyl response, while `1-a_bg` corresponds to decorrelated Weyl auto-power.

## Caveat

This is not yet an action-level derivation. The next required step is to derive the Weyl decomposition from the BPB/HR perturbation sector, or to show why the equality is only an empirical coincidence of the compressed E_G windows.

## Files written

- `paperVI_EG_T1_rWdelta_sqrt_abg_summary.csv`
- `paperVI_EG_T1_r_hypothesis_comparison.csv`
- `paperVI_EG_T1_decisions.csv`
- `paperVI_EG_T1_rWdelta_sqrt_abg_gate.json`
- `paperVI_EG_T1_rWdelta_sqrt_abg_gate.md`
- `fig_P6_T1_rWdelta_sqrt_abg.png/pdf`

"""
    (outdir / "paperVI_EG_T1_rWdelta_sqrt_abg_gate.md").write_text(md, encoding="utf-8")


def make_fig(outdir: Path, result: Dict[str, Any]) -> None:
    inp = result["inputs"]
    der = result["derived"]
    labels = ["free r", "sqrt(a_bg)", "1-a_bg", "sqrt(1-a_bg)", "a_bg"]
    rvals = [inp["r_wdelta_free"], der["sqrt_a_bg"], der["one_minus_a_bg"], der["sqrt_one_minus_a_bg"], inp["a_bg"]]
    chi2 = [inp["chi2_free_r"], inp["chi2_fixed_sqrt_abg"], inp["chi2_r_equals_1_minus_abg"], inp["chi2_r_equals_sqrt_1_minus_abg"], inp["chi2_r_equals_abg"]]

    fig, ax1 = plt.subplots(figsize=(9, 5.2))
    x = list(range(len(labels)))
    bars = ax1.bar(x, chi2, alpha=0.75)
    ax1.set_ylabel("Projected E_G diagnostic χ²")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=20, ha="right")
    ax1.set_title("Paper VI T1: r_Wδ = sqrt(a_bg) is degenerate with free r")
    ax1.axhline(inp["chi2_free_r"], linestyle="--", linewidth=1)

    for i, (b, rv, c2) in enumerate(zip(bars, rvals, chi2)):
        ax1.text(b.get_x()+b.get_width()/2, c2 + 0.35, f"r={rv:.3f}\nχ²={c2:.3f}", ha="center", va="bottom", fontsize=8)

    txt = (
        f"a_bg={inp['a_bg']:.6f}\n"
        f"sqrt(a_bg)={der['sqrt_a_bg']:.6f}\n"
        f"r_free={inp['r_wdelta_free']:.6f}\n"
        f"Δχ²(sqrt a_bg - free)={der['delta_chi2_fixed_sqrt_abg_minus_free']:.6f}"
    )
    ax1.text(0.98, 0.97, txt, transform=ax1.transAxes, ha="right", va="top", fontsize=9,
             bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.85))
    fig.tight_layout()
    fig.savefig(outdir / "fig_P6_T1_rWdelta_sqrt_abg.png", dpi=180)
    fig.savefig(outdir / "fig_P6_T1_rWdelta_sqrt_abg.pdf")
    plt.close(fig)


def main() -> None:
    repo = Path.cwd()
    outdir = repo / "outputs" / "paperVI_EG_theory"
    outdir.mkdir(parents=True, exist_ok=True)

    inp = LockedInputs()
    result = compute(inp)

    write_csvs(outdir, result)
    write_md(outdir, result)
    make_fig(outdir, result)
    with (outdir / "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json").open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(result["status"]["gate_status"])
    print("Wrote:")
    for p in sorted(outdir.glob("paperVI_EG_T1_*")):
        print(" -", p)
    print(" -", outdir / "fig_P6_T1_rWdelta_sqrt_abg.png")
    print(" -", outdir / "fig_P6_T1_rWdelta_sqrt_abg.pdf")


if __name__ == "__main__":
    main()
