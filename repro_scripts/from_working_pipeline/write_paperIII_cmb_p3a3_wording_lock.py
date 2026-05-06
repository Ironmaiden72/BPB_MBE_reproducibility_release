#!/usr/bin/env python3
"""
Paper III CMB P3-A3 — wording lock.

Non-destructive gate: reads the P3-A1b/P3-A2 audit artefacts and writes the
paper-facing wording lock to outputs/cmb_paperIII. It does not call ACT,
Planck, CAMB, Cobaya, or any NS code.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path.cwd()
OUT = ROOT / "outputs" / "cmb_paperIII"
OUT.mkdir(parents=True, exist_ok=True)

STATUS = "P3_A3_WORDING_LOCK_DONE"

# Canonical numerical backbone, locked from P3-A1b/P3-A2.
NUMBERS = {
    "ACT": {
        "Planck_PR4_chi2_ACT": 13.951140258982589,
        "A_vanilla_late_amp1_chi2_ACT": 227.5207684179947,
        "Kprimary_amp1_chi2_ACT": 17.470020679940834,
        "Kprimary_CMB8_amp_chi2_ACT": 14.156104008823837,
        "Kprimary_scan_best_amp_PhiPhi": 1.043,
        "Kprimary_scan_best_Sigma_eff": 1.0212737145349429,
        "Kprimary_scan_best_chi2_ACT": 14.15359945719008,
    },
    "PlanckLiteHighEll": {
        "Planck_PR4_best_full_chi2": 765.8691646821364,
        "Planck_PR4_best_full_ns": 0.972,
        "Planck_PR4_best_full_As_scale": 0.986,
        "C_s_surv0_best_full_chi2": 749.191767798833,
        "C_s_surv0_best_full_ns": 0.982,
        "C_s_surv0_best_full_As_scale": 0.976,
        "delta_vs_Planck": -16.67739688330346,
        "A_vanilla_best_full_chi2": 8851.940268986189,
        "A_vanilla_delta_vs_Planck": 8086.071104304053,
    },
    "PlanckLiteHighLowTau": {
        "C_s_surv0_total_chi2": 1164.6451180840354,
        "Planck_PR4_total_chi2": 1182.5539175383526,
        "delta_total": -17.908799454317204,
        "C_s_surv0_highl_chi2": 747.5775462736663,
        "Planck_PR4_highl_chi2": 764.9643487205054,
        "delta_highl": -17.38680244683917,
        "delta_lowT_lowE": -0.521997007478034,
    },
    "PlanckLensing": {
        "BPB_best_amp_phi": 1.075,
        "BPB_Sigma_eff": 1.0368220676663862,
        "BPB_best_chi2_lensing": 8.742460151124908,
        "Planck_PR4_best_chi2": 9.18193886999318,
        "Vanilla_best_chi2": 20.54451419576416,
        "ThreeBand_BPB_chi2": 2.235354820530869,
        "ThreeBand_Planck_chi2": 2.3432253250195854,
        "ThreeBand_Vanilla_chi2": 5.042760792159987,
    },
    "JointLensing": {
        "ACT_only_amp_phi": 1.043,
        "ACT_only_Sigma_eff": 1.0212737145349429,
        "ACT_only_chi2_ACT": 14.15359945719008,
        "ACT_only_chi2_Planck_lensing": 10.198317000740625,
        "ACT_only_joint": 24.351916457930702,
        "Planck_only_amp_phi": 1.075,
        "Planck_only_Sigma_eff": 1.0368220676663862,
        "Planck_only_chi2_ACT": 15.90870191352052,
        "Planck_only_chi2_Planck_lensing": 8.742460151124908,
        "Planck_only_joint": 24.651162064645423,
        "Joint_amp_phi": 1.057,
        "Joint_Sigma_eff": 1.02810505299799,
        "Joint_chi2_ACT": 14.478027412095692,
        "Joint_chi2_Planck_lensing": 9.243264395408051,
        "Joint_chi2_joint": 23.721291807503743,
        "separate_minima_sum": 22.896059608314985,
        "joint_scalar_penalty": 0.8252321991887577,
    },
}

WORDING_RULES = [
    {
        "id": "WR.ACT.01",
        "object": "ACT DR6 lensing",
        "locked_wording": "ACT is treated with the full ACT DR6 lenslike diagnostic used by the pipeline. The strict face-minus/Kprimary branch is close to the ACT baseline after a small positive Weyl-auto amplitude profiling, whereas the vanilla late-density CMB dictionary is catastrophically disfavoured by the same diagnostic.",
        "allowed_strength": "STRONG_DIAGNOSTIC",
        "do_not_say": "Do not call this a standalone full cosmological MCMC or a final combined CMB likelihood.",
    },
    {
        "id": "WR.PLANCKLITE.01",
        "object": "Planck-lite TTTEEE high ell",
        "locked_wording": "The Planck-lite TTTEEE high-ell diagnostic is a controlled proxy/lite comparison, not the official full Planck nuisance-foreground likelihood. Within this proxy, the strict face-minus C_s_surv0 dictionary improves the profiled high-ell chi2 relative to the Planck sanity point, while the vanilla dictionary fails catastrophically.",
        "allowed_strength": "STRONG_PROXY_WITH_EXPLICIT_LIMIT",
        "do_not_say": "Do not call it official full Planck, do not imply foreground/nuisance marginalization, and do not present it as a final Planck model-selection result.",
    },
    {
        "id": "WR.LOWELL.01",
        "object": "Planck-lite high+low tau/As",
        "locked_wording": "The native lowT/lowE addition is a supporting consistency gate. It does not introduce a low-ell penalty against the strict branch at the audited point; the high+low tau/As grid remains favourable to C_s_surv0 relative to the Planck sanity point in this diagnostic setup.",
        "allowed_strength": "SUPPORTING_DIAGNOSTIC",
        "do_not_say": "Do not overstate the low-ell module as a full Planck low-ell+nuisance production run.",
    },
    {
        "id": "WR.PLANCKLENS.01",
        "object": "Planck CMB-marginalized lensing",
        "locked_wording": "Planck CMB-marginalized lensing is used as an official-formula lensing diagnostic with scalar amplitude profiling. The strict branch prefers a modest positive lensing/Weyl-auto amplitude and is not penalized relative to the Planck sanity point; the vanilla mapping is again the poor control.",
        "allowed_strength": "CANONICAL_LENSING_DIAGNOSTIC",
        "do_not_say": "Do not call the three-band forensic scan the primary physical kernel fit.",
    },
    {
        "id": "WR.JOINTLENS.01",
        "object": "ACT+Planck scalar lensing closure",
        "locked_wording": "The joint ACT+Planck lensing comparison is a scalar-amplitude closure diagnostic. A single positive Weyl-auto amplitude, amp_phi=1.057 (Sigma_eff=1.028105), jointly fits ACT and Planck lensing with only Delta chi2=0.825 relative to the sum of their separate minima.",
        "allowed_strength": "SCALAR_CLOSURE_STRONG",
        "do_not_say": "Do not call it an official combined ACT+Planck likelihood, and do not infer a full window-dependent physical kernel from the scalar closure alone.",
    },
    {
        "id": "WR.CORECLAIM.01",
        "object": "Paper III CMB core claim",
        "locked_wording": "The CMB branch does not kill MAP1782. It kills the vanilla late-density CMB dictionary and selects the strict face-minus Weyl-sector hierarchy: primary-CMB dictionary, matter growth response, and Weyl/lensing response must be kept distinct.",
        "allowed_strength": "CORE_CLAIM",
        "do_not_say": "Do not weaken this into 'inconclusive'; the correct limit is that the claim is diagnostic/lite where appropriate, not null.",
    },
]

CLAIM_WORDING = [
    {
        "claim_id": "P3.ACT.01",
        "paper_wording": "In the ACT DR6 lenslike diagnostic, the vanilla late-density projection is ruled out as a viable CMB dictionary: chi2_ACT=227.520768 versus 13.951140 for the Planck sanity point. The strict face-minus/Kprimary branch instead sits close to the ACT baseline, with chi2_ACT=17.470021 at unit amplitude and 14.153599 after scalar Weyl-auto profiling at amp_PhiPhi=1.043 (Sigma_eff=1.021274).",
        "strength": "strong diagnostic",
        "source_gate": "P3-A1/P3-A2",
    },
    {
        "claim_id": "P3.PLANCKLITE.01",
        "paper_wording": "In the Planck-lite high-ell TTTEEE proxy, the strict C_s_surv0 dictionary gives best_full_chi2=749.191768, improving over the Planck_PR4 sanity point by Delta chi2=-16.677397, while the vanilla dictionary gives chi2=8851.940269 and is catastrophically excluded within the same proxy.",
        "strength": "strong proxy, explicitly not official full Planck",
        "source_gate": "P3-A1/P3-A2",
    },
    {
        "claim_id": "P3.PLANCKLOW.01",
        "paper_wording": "Adding native lowT/lowE and profiling tau/As does not generate a compensating low-ell failure: C_s_surv0 gives total chi2=1164.645118 versus 1182.553918 for the Planck_PR4 sanity point, Delta chi2=-17.908799, with the lowT+lowE contribution slightly favourable by about -0.522.",
        "strength": "supporting diagnostic",
        "source_gate": "P3-A2",
    },
    {
        "claim_id": "P3.PLANCKLENS.01",
        "paper_wording": "In the Planck CMB-marginalized lensing diagnostic, the strict branch is compatible with a modest positive Weyl-auto/lensing amplitude: best_amp_phi=1.075, Sigma_eff=1.036822, chi2=8.742460, compared with 9.181939 for the Planck sanity point and 20.544514 for the vanilla control.",
        "strength": "canonical lensing diagnostic",
        "source_gate": "P3-A1/P3-A2",
    },
    {
        "claim_id": "P3.JOINTLENS.01",
        "paper_wording": "The ACT and Planck lensing amplitudes close with one common scalar parameter: the joint scalar best fit is amp_phi=1.057, Sigma_eff=1.028105, chi2_joint=23.721292. The penalty relative to the sum of the separate ACT-only and Planck-only minima is only Delta chi2=0.825232.",
        "strength": "joint scalar-amplitude closure diagnostic",
        "source_gate": "P3-A1b/P3-A2",
    },
    {
        "claim_id": "P3.CORE.01",
        "paper_wording": "Across ACT lensing, Planck-lite high-ell/low-ell diagnostics, and CMB-marginalized lensing, the consistent lesson is not that BPB/MBE fails in the CMB. The failure mode is the naive vanilla dictionary. The strict face-minus branch remains viable at diagnostic level and selects a response hierarchy in which primary CMB, matter growth, and Weyl lensing are not collapsed into a single amplitude.",
        "strength": "synthesis claim",
        "source_gate": "P3-A3",
    },
]

TABLE_CAPTIONS = [
    {
        "table_id": "T_CMB_ACT_DR6_lenslike",
        "caption": "ACT DR6 lenslike diagnostic for the MAP1782 CMB dictionaries. The vanilla late-density mapping is shown as the failure control; the strict face-minus/Kprimary branch is close to the Planck sanity point and admits a small positive scalar Weyl-auto amplitude profile.",
    },
    {
        "table_id": "T_CMB_PlanckLite_high_ell",
        "caption": "Planck-lite high-ell TTTEEE proxy comparison. This is not the official full Planck nuisance/foreground likelihood. It is used to test whether each CMB dictionary preserves the high-ell acoustic structure under the same lite/proxy assumptions.",
    },
    {
        "table_id": "T_CMB_PlanckLite_high_low_tau",
        "caption": "Supporting high-ell plus native lowT/lowE tau/As grid. The low-ell modules are used as a consistency gate and do not introduce a compensating penalty against the strict face-minus branch at the audited point.",
    },
    {
        "table_id": "T_CMB_Planck_lensing_CMBmarged",
        "caption": "Planck CMB-marginalized lensing official-formula diagnostic with scalar amplitude profiling. The strict branch prefers a modest positive lensing/Weyl-auto amplitude; the vanilla mapping remains the poor control.",
    },
    {
        "table_id": "T_CMB_Joint_ACT_Planck_scalar_lensing",
        "caption": "Joint ACT+Planck scalar lensing closure. This is a scalar-amplitude diagnostic, not an official combined likelihood. The small penalty relative to separate ACT-only and Planck-only minima tests whether one common Weyl-auto response can serve both windows.",
    },
]

FIGURE_CAPTIONS = [
    {
        "figure_id": "Fig_P3_1_ACT_lensing_gate",
        "caption": "ACT DR6 lenslike CMB dictionary gate. The vanilla late-density CMB mapping is catastrophically disfavoured, while the strict face-minus/Kprimary dictionary remains close to the Planck sanity point and improves under a small positive Weyl-auto amplitude profile.",
        "must_show": "Planck_PR4, A_vanilla_late_amp1, Kprimary_amp1, Kprimary scan best; annotate chi2_ACT values.",
    },
    {
        "figure_id": "Fig_P3_2_PlanckLite_high_ell",
        "caption": "Planck-lite high-ell TTTEEE proxy. The strict C_s_surv0 dictionary improves the profiled high-ell chi2 relative to the Planck sanity point in this lite diagnostic, while the vanilla dictionary fails by thousands in chi2.",
        "must_show": "Planck_PR4, C_s_surv0, optional C_Zbest, vanilla; visibly separate vanilla due to scale or broken axis/log inset.",
    },
    {
        "figure_id": "Fig_P3_3_joint_lensing_scalar_closure",
        "caption": "Scalar-amplitude closure for ACT DR6 and Planck CMB-marginalized lensing. ACT-only and Planck-only windows prefer nearby positive amplitudes, and one common amp_phi=1.057 gives a joint penalty of only Delta chi2=0.825 relative to separate minima.",
        "must_show": "ACT-only, Planck-only, joint scalar best; amp_phi, Sigma_eff, chi2_ACT, chi2_Planck_lensing, joint penalty.",
    },
    {
        "figure_id": "Fig_P3_4_optional_proxy_contact_sheet",
        "caption": "Optional proxy contact sheet collecting TT/TE/EE shape, peak-position, and lensing/Weyl proxy diagnostics. This figure is supporting only and should not replace the three canonical table-driven figures.",
        "must_show": "Only if needed in appendix/supplement; label clearly as proxy/diagnostic.",
    },
]

FORBIDDEN_PHRASES = [
    {"phrase": "official full Planck likelihood", "replacement": "Planck-lite/proxy high-ell diagnostic, unless referring to a future required gate"},
    {"phrase": "full Planck nuisance marginalization", "replacement": "not performed here; future official pipeline gate"},
    {"phrase": "official combined ACT+Planck likelihood", "replacement": "joint scalar-amplitude diagnostic closure"},
    {"phrase": "Planck proves BPB", "replacement": "Planck-lite/lensing diagnostics do not kill the strict branch and reject the vanilla dictionary"},
    {"phrase": "inconclusive CMB result", "replacement": "diagnostic-level but strong dictionary-selection result"},
    {"phrase": "the Weyl response is the matter response", "replacement": "primary CMB, matter growth, and Weyl/lensing response are distinct"},
    {"phrase": "three-band kernel is the final physical kernel", "replacement": "three-band scan is supporting forensic/kernel diagnostic"},
]

ABSTRACT_PARAGRAPH = (
    "The CMB sector provides a dictionary-selection test for the MAP1782 branch. "
    "Across ACT DR6 lensing, Planck-lite high-ell/low-ell diagnostics, and Planck CMB-marginalized lensing, the naive vanilla late-density projection fails, while the strict face-minus hierarchy remains viable at diagnostic level. "
    "The ACT lenslike gate rejects the vanilla mapping with chi2_ACT=227.520768 versus 13.951140 for the Planck sanity point, whereas the strict Kprimary branch reaches chi2_ACT=14.153599 after scalar Weyl-auto profiling. "
    "The Planck-lite high-ell proxy similarly favours C_s_surv0 over the Planck sanity point by Delta chi2=-16.677397 and gives a catastrophic chi2=8851.940269 for the vanilla control. "
    "Finally, ACT and Planck lensing close with one common scalar Weyl-auto amplitude, amp_phi=1.057 (Sigma_eff=1.028105), with a joint penalty of only Delta chi2=0.825232 relative to separate minima. "
    "These results do not constitute a full official Planck nuisance-likelihood claim; they show that the CMB does not kill MAP1782, but kills the vanilla CMB dictionary and selects distinct primary-CMB, matter-growth, and Weyl/lensing responses."
)

CONCLUSION_PARAGRAPH = (
    "The CMB result is therefore sharper than a generic consistency statement. "
    "The vanilla late-density dictionary is the failed branch. The strict face-minus dictionary is the surviving branch across the audited diagnostics. "
    "Because the lensing closure requires a modest positive Weyl-auto amplitude and because ACT, Planck-lite, and Planck lensing probe different windows, the paper should state the response hierarchy explicitly: the primary-CMB dictionary, the low-redshift matter response, and the Weyl/lensing response are not interchangeable amplitudes. "
    "The remaining caveat is procedural rather than conceptual: a future full official Planck nuisance/foreground MCMC and a physical Weyl-window kernel are required before claiming final CMB Bayesian model selection."
)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def md_table(rows: list[dict], columns: list[str]) -> str:
    out = []
    out.append("| " + " | ".join(columns) + " |")
    out.append("|" + "|".join(["---"] * len(columns)) + "|")
    for row in rows:
        out.append("| " + " | ".join(str(row.get(c, "")).replace("\n", " ") for c in columns) + " |")
    return "\n".join(out)


def main() -> None:
    # Passive existence checks only. Missing inputs do not block creation because the
    # lock is intentionally based on P3-A1b/P3-A2 numerical decisions.
    expected_inputs = [
        OUT / "paperIII_cmb_p3a2_canonical_tables.csv",
        OUT / "paperIII_cmb_p3a2_figure_plan.csv",
        OUT / "paperIII_cmb_joint_lens_v4_patch.json",
    ]
    input_status = [
        {"path": str(p), "exists": p.exists(), "role": "expected upstream P3-A1b/P3-A2 artefact"}
        for p in expected_inputs
    ]

    open_items = [
        {
            "id": "P3A3.O01",
            "item": "Generate publication-grade figures from locked table/figure captions",
            "status": "OPEN_NEXT",
            "note": "Recommended order: ACT gate, Planck-lite high-ell, joint scalar closure.",
        },
        {
            "id": "P3A3.O02",
            "item": "Use P3-A1b joint patch as authoritative source if old source_tables summary still says joint_best=None",
            "status": "OPEN_LOW_RISK",
            "note": "Do not propagate the legacy extraction bug into paper tables.",
        },
        {
            "id": "P3A3.O03",
            "item": "Future full official Planck nuisance/foreground MCMC",
            "status": "FUTURE_GATE_NOT_REQUIRED_FOR_CURRENT_WORDING",
            "note": "Explicit caveat retained; does not block Paper III diagnostic claim.",
        },
    ]

    summary = {
        "status": STATUS,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "safety": {
            "read_scope": "outputs/cmb_paperIII audit artefacts only, if present",
            "write_scope": "outputs/cmb_paperIII only",
            "likelihood_run": "NO",
            "NS_touched": "NO",
            "ACT_or_Planck_code_run": "NO",
        },
        "numbers": NUMBERS,
        "wording_rules": WORDING_RULES,
        "claim_wording": CLAIM_WORDING,
        "table_captions": TABLE_CAPTIONS,
        "figure_captions": FIGURE_CAPTIONS,
        "forbidden_phrases": FORBIDDEN_PHRASES,
        "abstract_paragraph": ABSTRACT_PARAGRAPH,
        "conclusion_paragraph": CONCLUSION_PARAGRAPH,
        "input_status": input_status,
        "open_items": open_items,
        "gate_decisions": {
            "P3-A3 wording lock": "PASS",
            "ACT wording": "LOCKED_FULL_LENSLIKE_DIAGNOSTIC",
            "Planck-lite wording": "LOCKED_PROXY_LITE_NOT_OFFICIAL_FULL_PLANCK",
            "Planck lensing wording": "LOCKED_OFFICIAL_FORMULA_CMBMARGED_DIAGNOSTIC",
            "Joint lensing wording": "LOCKED_SCALAR_AMPLITUDE_CLOSURE_NOT_OFFICIAL_COMBINED_LIKELIHOOD",
            "Core claim": "LOCKED_CMB_KILLS_VANILLA_NOT_STRICT_BRANCH",
            "next_gate": "P3-A4 figure generation",
        },
    }

    write_csv(OUT / "paperIII_cmb_p3a3_wording_rules.csv", WORDING_RULES)
    write_csv(OUT / "paperIII_cmb_p3a3_claim_wording.csv", CLAIM_WORDING)
    write_csv(OUT / "paperIII_cmb_p3a3_table_captions.csv", TABLE_CAPTIONS)
    write_csv(OUT / "paperIII_cmb_p3a3_figure_captions.csv", FIGURE_CAPTIONS)
    write_csv(OUT / "paperIII_cmb_p3a3_forbidden_phrases.csv", FORBIDDEN_PHRASES)
    write_csv(OUT / "paperIII_cmb_p3a3_input_status.csv", input_status)
    write_csv(OUT / "paperIII_cmb_p3a3_open_items.csv", open_items)
    (OUT / "paperIII_cmb_p3a3_wording_lock.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    md = []
    md.append("# Paper III CMB P3-A3 — wording lock")
    md.append("")
    md.append("## Status")
    md.append("")
    md.append(f"**{STATUS}**")
    md.append("")
    md.append("Safety: read-only on upstream CMB audit artefacts; writes only to `outputs/cmb_paperIII`; no ACT/Planck likelihood run; no NS touch.")
    md.append("")
    md.append("## One-sentence locked claim")
    md.append("")
    md.append("> The CMB does not kill MAP1782; it kills the vanilla late-density CMB dictionary and selects the strict face-minus Weyl-sector hierarchy, with primary-CMB, matter-growth, and Weyl/lensing responses kept distinct.")
    md.append("")
    md.append("## Locked abstract paragraph")
    md.append("")
    md.append(ABSTRACT_PARAGRAPH)
    md.append("")
    md.append("## Wording rules")
    md.append("")
    md.append(md_table(WORDING_RULES, ["id", "object", "locked_wording", "allowed_strength", "do_not_say"]))
    md.append("")
    md.append("## Claim wording")
    md.append("")
    md.append(md_table(CLAIM_WORDING, ["claim_id", "paper_wording", "strength", "source_gate"]))
    md.append("")
    md.append("## Table captions")
    md.append("")
    md.append(md_table(TABLE_CAPTIONS, ["table_id", "caption"]))
    md.append("")
    md.append("## Figure captions")
    md.append("")
    md.append(md_table(FIGURE_CAPTIONS, ["figure_id", "caption", "must_show"]))
    md.append("")
    md.append("## Forbidden / replacement phrases")
    md.append("")
    md.append(md_table(FORBIDDEN_PHRASES, ["phrase", "replacement"]))
    md.append("")
    md.append("## Locked conclusion paragraph")
    md.append("")
    md.append(CONCLUSION_PARAGRAPH)
    md.append("")
    md.append("## Gate decisions")
    md.append("")
    decisions = [{"decision": k, "value": v} for k, v in summary["gate_decisions"].items()]
    md.append(md_table(decisions, ["decision", "value"]))
    md.append("")
    md.append("## Open items")
    md.append("")
    md.append(md_table(open_items, ["id", "item", "status", "note"]))
    md.append("")
    md.append("## Outputs")
    md.append("")
    for name in [
        "paperIII_cmb_p3a3_wording_lock.md",
        "paperIII_cmb_p3a3_wording_lock.json",
        "paperIII_cmb_p3a3_wording_rules.csv",
        "paperIII_cmb_p3a3_claim_wording.csv",
        "paperIII_cmb_p3a3_table_captions.csv",
        "paperIII_cmb_p3a3_figure_captions.csv",
        "paperIII_cmb_p3a3_forbidden_phrases.csv",
        "paperIII_cmb_p3a3_open_items.csv",
    ]:
        md.append(f"- `outputs/cmb_paperIII/{name}`")
    md.append("")
    (OUT / "paperIII_cmb_p3a3_wording_lock.md").write_text("\n".join(md), encoding="utf-8")

    print(f"Saved: {OUT / 'paperIII_cmb_p3a3_wording_lock.md'}")
    print(f"Saved: {OUT / 'paperIII_cmb_p3a3_wording_lock.json'}")
    print(f"Status: {STATUS}")


if __name__ == "__main__":
    main()
