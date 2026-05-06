#!/usr/bin/env python3
"""
Paper VI / E_G T10b — density-weighted HR/BPB mixing LaTeX insertion pack.

Non-destructive by default. It upgrades the T10/T9A insertion pack to the
T9D/T9E wording:
  - T9D: conditional density-weighted mixing derivation
  - T9E: D2 source-match not found internally, so D2 remains an explicit guardrail

No likelihood run. Does not touch NS paths.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Dict, List


OUTDIR = Path("outputs/paperVI_EG_theory")
OUTDIR.mkdir(parents=True, exist_ok=True)

T9D_VALUES = OUTDIR / "paperVI_EG_T9D_key_values.csv"
T9D_MD = OUTDIR / "paperVI_EG_T9D_density_weighted_mixing_closure.md"
T9E_MD = OUTDIR / "paperVI_EG_T9E_source_match_density_weighted_mixing.md"
T9E_DECISIONS = OUTDIR / "paperVI_EG_T9E_decisions.csv"

SECTION_TEX = OUTDIR / "paperVI_EG_T10b_density_weighted_closure_section.tex"
TABLE_TEX = OUTDIR / "paperVI_EG_T10b_density_weighted_values_table.tex"
SOURCE_LOCK = OUTDIR / "paperVI_EG_T10b_source_lock_final.csv"
DECISIONS = OUTDIR / "paperVI_EG_T10b_decisions.csv"
MANIFEST_MD = OUTDIR / "paperVI_EG_T10b_manifest.md"
MANIFEST_JSON = OUTDIR / "paperVI_EG_T10b_manifest.json"
KEY_VALUES_OUT = OUTDIR / "paperVI_EG_T10b_key_values.csv"


def read_key_values(path: Path) -> Dict[str, str]:
    vals: Dict[str, str] = {}
    if not path.exists():
        return vals
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = row.get("quantity") or row.get("key") or row.get("name")
            v = row.get("value")
            if k is not None and v is not None:
                vals[k] = v
    return vals


def fval(vals: Dict[str, str], key: str, default: str) -> str:
    return vals.get(key, default)


def write_csv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def maybe_apply_to_paper(section_path: Path, paper_tex: str, apply: bool) -> Dict[str, str]:
    result = {
        "applied": "False",
        "paper_tex": paper_tex or "",
        "backup": "",
        "apply_status": "NOT_REQUESTED",
    }

    if not apply:
        return result

    if not paper_tex:
        result["apply_status"] = "ERROR_NO_TARGET"
        return result

    target = Path(paper_tex)
    if not target.exists():
        result["apply_status"] = "ERROR_TARGET_NOT_FOUND"
        return result

    insertion = r"\input{outputs/paperVI_EG_theory/paperVI_EG_T10b_density_weighted_closure_section.tex}"
    text = target.read_text(encoding="utf-8", errors="replace")

    if insertion in text:
        result["apply_status"] = "ALREADY_PRESENT"
        result["applied"] = "False"
        return result

    backup = target.with_suffix(target.suffix + ".p6t10b.bak")
    shutil.copy2(target, backup)

    marker = r"\end{document}"
    if marker in text:
        text = text.replace(marker, insertion + "\n\n" + marker, 1)
    else:
        text = text + "\n\n" + insertion + "\n"

    target.write_text(text, encoding="utf-8")
    result["applied"] = "True"
    result["backup"] = str(backup)
    result["apply_status"] = "APPLIED"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper-tex", default="", help="Optional Paper VI .tex target.")
    parser.add_argument("--apply", action="store_true", help="Apply insertion to --paper-tex with backup.")
    args = parser.parse_args()

    vals = read_key_values(T9D_VALUES)

    # Fallbacks from frozen T9D/T9E console output.
    a_bg = fval(vals, "a_bg", "0.403687948")
    one_minus_a_bg = fval(vals, "1_minus_a_bg", "0.596312052")
    sqrt_a_bg = fval(vals, "sqrt_a_bg", "0.6353644214149861")
    sqrt_one_minus = fval(vals, "sqrt_one_minus_a_bg", "0.7722124396822418")
    rho_ratio = fval(vals, "rho_f_over_rho_g", "1.4771608985463198")
    theta_deg = fval(vals, "theta_deg", "50.55298183517175")
    cos2_minus = fval(vals, "cos2_theta_minus_a_bg", "0.0")
    free_r = fval(vals, "free_r_Wdelta", "0.630864746")
    delta_r = fval(vals, "delta_r_free_minus_cos_theta", "-0.004499675414986082")
    chi2_free = fval(vals, "chi2_free_r", "1.358943")
    chi2_fixed = fval(vals, "chi2_fixed_density_mixing", "1.367488")
    delta_chi2 = fval(vals, "delta_chi2_fixed_minus_free", "0.008545000000000025")
    t8_fraction = fval(vals, "T8_fraction_XX_over_AA", "0.3979903273355348")
    t8_sqrt = fval(vals, "T8_sqrt_fraction", "0.6308647456749622")
    t8_delta = fval(vals, "T8_fraction_minus_a_bg", "-0.005697620664465219")
    t8_rel = fval(vals, "T8_relative_fraction_delta", "-0.014113923124762742")
    t8_chi2pt = fval(vals, "T8_all_chi2_per_point", "0.22649042932465788")

    table_tex = rf"""\begin{{table}}
\centering
\caption{{Density-weighted projected Weyl--matter closure values for Paper VI. The closure is conditional on the density-weighted HR/BPB scalar Weyl mixing relation $ \tan^2\theta_W=\rho_f^{{\rm eff}}/\rho_g^{{\rm eff}} $, whose explicit perturbation-equation source match remains open.}}
\label{{tab:paperVI-eg-density-weighted-rwdelta}}
\begin{{tabular}}{{lc}}
\hline
Quantity & Value \\
\hline
$a_{{\rm bg}}$ & {a_bg} \\
$1-a_{{\rm bg}}$ & {one_minus_a_bg} \\
$\sqrt{{a_{{\rm bg}}}}$ & {sqrt_a_bg} \\
$\sqrt{{1-a_{{\rm bg}}}}$ & {sqrt_one_minus} \\
$\rho_f^{{\rm eff}}/\rho_g^{{\rm eff}}=(1-a_{{\rm bg}})/a_{{\rm bg}}$ & {rho_ratio} \\
$\theta_W$ & {theta_deg} deg \\
$\cos^2\theta_W-a_{{\rm bg}}$ & {cos2_minus} \\
free $r_{{W\delta}}$ & {free_r} \\
$r_{{W\delta}}^{{\rm free}}-\sqrt{{a_{{\rm bg}}}}$ & {delta_r} \\
$\chi^2_{{\rm free}}$ & {chi2_free} \\
$\chi^2_{{\sqrt{{a_{{\rm bg}}}}}}$ & {chi2_fixed} \\
$\Delta\chi^2$ fixed-free & {delta_chi2} \\
$\langle X,X\rangle/\langle A,A\rangle$ & {t8_fraction} \\
$\sqrt{{\langle X,X\rangle/\langle A,A\rangle}}$ & {t8_sqrt} \\
$\langle X,X\rangle/\langle A,A\rangle-a_{{\rm bg}}$ & {t8_delta} \\
all-row support $\chi^2/N$ & {t8_chi2pt} \\
\hline
\end{{tabular}}
\end{{table}}
"""
    TABLE_TEX.write_text(table_tex, encoding="utf-8")

    section_tex = rf"""% BEGIN P6_T10B_DENSITY_WEIGHTED_WEYL_MATTER_CLOSURE
\subsection{{Density-weighted projected Weyl--matter closure}}
\label{{sec:paperVI-eg-density-weighted-rwdelta}}

The compressed $E_G$ auto/cross diagnostic introduced an effective Weyl--matter correlation coefficient $r_{{W\delta}}$. The projected diagnostic result is that the free value,
\begin{{equation}}
  r_{{W\delta}}^{{\rm free}} = {free_r},
\end{{equation}}
is reproduced by the background-dictionary closure
\begin{{equation}}
  r_{{W\delta}}=\sqrt{{a_{{\rm bg}}}}={sqrt_a_bg},
  \qquad
  a_{{\rm bg}}={a_bg},
\end{{equation}}
with a negligible diagnostic penalty,
\begin{{equation}}
  \Delta\chi^2 = {delta_chi2}.
\end{{equation}}

A useful way to interpret this closure is through a two-sector Weyl mixing angle. Let the effective sector weights entering the scalar Weyl response at the biface transition be
\begin{{equation}}
  \rho_g^{{\rm eff}}:\rho_f^{{\rm eff}}
  =
  a_{{\rm bg}}:(1-a_{{\rm bg}}).
\end{{equation}}
If the two scalar Weyl modes obey the density-weighted mixing relation
\begin{{equation}}
  \tan^2\theta_W =
  \frac{{\rho_f^{{\rm eff}}}}{{\rho_g^{{\rm eff}}}},
  \label{{eq:paperVI-D2-density-weighted-mixing}}
\end{{equation}}
then
\begin{{equation}}
  \tan^2\theta_W
  =
  \frac{{1-a_{{\rm bg}}}}{{a_{{\rm bg}}}}.
\end{{equation}}
It follows algebraically that
\begin{{equation}}
  \cos^2\theta_W
  =
  \frac{{1}}{{1+\tan^2\theta_W}}
  =
  \frac{{1}}{{1+(1-a_{{\rm bg}})/a_{{\rm bg}}}}
  =
  a_{{\rm bg}},
\end{{equation}}
and therefore
\begin{{equation}}
  r_{{W\delta}}=\cos\theta_W=\sqrt{{a_{{\rm bg}}}}.
\end{{equation}}

This closure is non-circular at the projected diagnostic level: $a_{{\rm bg}}$ is inherited from the BPB background dictionary, not fitted to $E_G$. The $E_G$ auto/cross diagnostic then tests the consequence. In the amplitude-space artefact, define
\begin{{equation}}
  A_i={\tt unit\_ACTSigma\_cross}_i,
  \qquad
  X_i={\tt prediction\_global\_auto\_cross}_i .
\end{{equation}}
The row-level relation is $X_i=r_{{\rm global}}A_i$, hence for positive bin weights
\begin{{equation}}
  \frac{{\langle X,X\rangle}}{{\langle A,A\rangle}}
  =
  r_{{\rm global}}^2.
\end{{equation}}
The audited value is
\begin{{equation}}
  \frac{{\langle X,X\rangle}}{{\langle A,A\rangle}}
  =
  {t8_fraction},
  \qquad
  \sqrt{{\frac{{\langle X,X\rangle}}{{\langle A,A\rangle}}}}
  =
  {t8_sqrt},
\end{{equation}}
which differs from $a_{{\rm bg}}$ by {t8_delta}, i.e. {t8_rel} in relative units. The all-row compressed support is $\chi^2/N={t8_chi2pt}$.

Equation~\eqref{{eq:paperVI-D2-density-weighted-mixing}} is the remaining field-equation guardrail. The internal source-match audit did not find a non-circular derivation of this density-weighted scalar/Weyl mixing relation in the existing HR/BPB perturbation files. Thus the result should be stated as a conditional density-weighted HR/BPB mixing closure: it upgrades the purely projected diagnostic relation, but it does not yet constitute a fully source-matched action-level perturbation proof.

\input{{outputs/paperVI_EG_theory/paperVI_EG_T10b_density_weighted_values_table.tex}}
% END P6_T10B_DENSITY_WEIGHTED_WEYL_MATTER_CLOSURE
"""
    SECTION_TEX.write_text(section_tex, encoding="utf-8")

    key_rows = [
        {"quantity": "a_bg", "value": a_bg},
        {"quantity": "1_minus_a_bg", "value": one_minus_a_bg},
        {"quantity": "sqrt_a_bg", "value": sqrt_a_bg},
        {"quantity": "sqrt_one_minus_a_bg", "value": sqrt_one_minus},
        {"quantity": "rho_f_over_rho_g", "value": rho_ratio},
        {"quantity": "theta_deg", "value": theta_deg},
        {"quantity": "free_r_Wdelta", "value": free_r},
        {"quantity": "delta_chi2_fixed_minus_free", "value": delta_chi2},
        {"quantity": "T8_fraction_XX_over_AA", "value": t8_fraction},
        {"quantity": "T8_sqrt_fraction", "value": t8_sqrt},
        {"quantity": "T8_fraction_minus_a_bg", "value": t8_delta},
        {"quantity": "T8_relative_fraction_delta", "value": t8_rel},
        {"quantity": "T8_all_chi2_per_point", "value": t8_chi2pt},
    ]
    write_csv(KEY_VALUES_OUT, key_rows, ["quantity", "value"])

    source_rows = [
        {
            "claim_id": "P6.T10B.01",
            "claim": "Density-weighted two-sector mixing gives r_Wdelta=sqrt(a_bg) if tan^2(theta_W)=rho_f_eff/rho_g_eff.",
            "source": "outputs/paperVI_EG_theory/paperVI_EG_T9D_density_weighted_mixing_closure.md/csv",
            "status": "SOURCE_MATCHED_CONDITIONAL_DERIVATION",
        },
        {
            "claim_id": "P6.T10B.02",
            "claim": "a_bg is inherited from the BPB background dictionary and not fitted to E_G.",
            "source": "outputs/paperVI_EG_theory/paperVI_EG_T9D_density_weighted_mixing_closure.md",
            "status": "NON_CIRCULARITY_LOCKED_AT_DIAGNOSTIC_LEVEL",
        },
        {
            "claim_id": "P6.T10B.03",
            "claim": "E_G auto/cross amplitude space verifies r_global^2 close to a_bg with negligible fixed-free penalty.",
            "source": "outputs/paperVI_EG_theory/paperVI_EG_T8_auto_cross_amplitude_space.md/csv",
            "status": "SOURCE_MATCHED_PROJECTED_DIAGNOSTIC",
        },
        {
            "claim_id": "P6.T10B.04",
            "claim": "Internal source-match did not close D2: tan^2(theta_W)=rho_f_eff/rho_g_eff remains a labelled guardrail.",
            "source": "outputs/paperVI_EG_theory/paperVI_EG_T9E_source_match_density_weighted_mixing.md/csv",
            "status": "OPEN_GUARDRAIL_LOCKED",
        },
        {
            "claim_id": "P6.T10B.05",
            "claim": "LaTeX v0.2 insertion pack generated for Paper VI.",
            "source": "outputs/paperVI_EG_theory/paperVI_EG_T10b_density_weighted_closure_section.tex",
            "status": "LATEX_INSERTION_READY",
        },
    ]
    write_csv(SOURCE_LOCK, source_rows, ["claim_id", "claim", "source", "status"])

    decisions = [
        {
            "item": "insertion_pack_v02",
            "decision": "READY",
            "basis": "T9D/T9E wording converted to LaTeX section and table.",
        },
        {
            "item": "claim_strength",
            "decision": "CONDITIONAL_DENSITY_WEIGHTED_HR_BPB_MIXING_CLOSURE",
            "basis": "D2 is explicit and remains the field-equation source-match guardrail.",
        },
        {
            "item": "full_action_level_proof",
            "decision": "NOT_CLAIMED",
            "basis": "T9E did not source-match tan^2(theta_W)=rho_f_eff/rho_g_eff internally.",
        },
        {
            "item": "paper_use",
            "decision": "USE_FOR_PAPER_VI_V02",
            "basis": "Stronger than T9A but still honest about D2.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T11_INSERT_AND_COMPILE",
            "basis": "Insert section/table into Paper VI and check PDF.",
        },
    ]
    write_csv(DECISIONS, decisions, ["item", "decision", "basis"])

    apply_result = maybe_apply_to_paper(SECTION_TEX, args.paper_tex, args.apply)

    manifest = {
        "status": "P6_T10B_DENSITY_WEIGHTED_LATEX_INSERTION_PACK_DONE",
        "likelihood_run": "NO",
        "ns_touched": "NO",
        "files_written": [
            str(SECTION_TEX),
            str(TABLE_TEX),
            str(SOURCE_LOCK),
            str(DECISIONS),
            str(MANIFEST_MD),
            str(MANIFEST_JSON),
            str(KEY_VALUES_OUT),
        ],
        "apply_result": apply_result,
        "guardrail": "D2 remains open: tan^2(theta_W)=rho_f_eff/rho_g_eff must still be source-matched in explicit HR/BPB perturbation conventions.",
    }
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    manifest_md = rf"""# Paper VI / E_G T10b - density-weighted LaTeX insertion pack

## Status

**P6_T10B_DENSITY_WEIGHTED_LATEX_INSERTION_PACK_DONE**

Non-destructive by default. No likelihood run. NS untouched.

## Purpose

T10b upgrades the T10/T9A projected-diagnostic insertion to the T9D/T9E wording:

- T9D gives a conditional density-weighted HR/BPB mixing derivation.
- T9E confirms that the remaining D2 relation is not internally source-matched in existing non-generated files.

## Apply status

| item | value |
|---|---|
| applied to paper tex | {apply_result["applied"]} |
| target paper tex | `{apply_result["paper_tex"]}` |
| apply status | {apply_result["apply_status"]} |
| backup | `{apply_result["backup"]}` |

## Files written

- `paperVI_EG_T10b_density_weighted_closure_section.tex`
- `paperVI_EG_T10b_density_weighted_values_table.tex`
- `paperVI_EG_T10b_source_lock_final.csv`
- `paperVI_EG_T10b_decisions.csv`
- `paperVI_EG_T10b_manifest.md`
- `paperVI_EG_T10b_manifest.json`
- `paperVI_EG_T10b_key_values.csv`

## Recommended use

In the Paper VI LaTeX file, add:

```tex
\input{{outputs/paperVI_EG_theory/paperVI_EG_T10b_density_weighted_closure_section.tex}}
```

Recommended insertion point: after the compressed E_G diagnostic results and before the final discussion/conclusion.

## Claim strength

This section states a conditional density-weighted HR/BPB mixing closure:

```tex
r_{{W\delta}}=\sqrt{{a_{{\rm bg}}}}
```

under the explicit D2 condition:

```tex
\tan^2\theta_W=\frac{{\rho_f^{{\rm eff}}}}{{\rho_g^{{\rm eff}}}}.
```

## Guardrail

T9E did **not** source-match D2 internally. Therefore the section must not claim a completed action-level HR/BPB perturbation proof.

The open field-equation target remains to derive D2 from explicit HR/BPB scalar/Weyl perturbation conventions.
"""
    MANIFEST_MD.write_text(manifest_md, encoding="utf-8")

    print("P6_T10B_DENSITY_WEIGHTED_LATEX_INSERTION_PACK_DONE")
    print("Wrote:", SECTION_TEX)
    print("Wrote:", TABLE_TEX)
    print("Wrote:", SOURCE_LOCK)
    print("Apply:", apply_result)


if __name__ == "__main__":
    main()
