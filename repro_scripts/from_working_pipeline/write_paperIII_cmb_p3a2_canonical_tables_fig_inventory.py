#!/usr/bin/env python3
"""
Paper III CMB P3-A2 — canonical table choice + figure inventory.

Non-destructive gate:
- reads outputs/cmb_paperIII source-match artefacts and CMB output inventory;
- writes only to outputs/cmb_paperIII;
- does not run ACT/Planck likelihood code;
- does not touch outputs/ns, src/ns, or NS checkpoints.

Run from repo root:
    python runs/cmb/write_paperIII_cmb_p3a2_canonical_tables_fig_inventory.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

OUTDIR = Path("outputs/cmb_paperIII")
SOURCE_TABLES = OUTDIR / "paperIII_cmb_source_tables.csv"
RELEVANT_FILES = OUTDIR / "paperIII_cmb_relevant_files.csv"
JOINT_PATCH = OUTDIR / "paperIII_cmb_joint_lens_v4_patch.csv"
JOINT_PATCH_JSON = OUTDIR / "paperIII_cmb_joint_lens_v4_patch.json"

STATUS = "P3_A2_CANONICAL_TABLES_FIGURE_INVENTORY_DONE"


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path}")
    return pd.read_csv(path)


def _json_row(row: pd.Series) -> Dict[str, Any]:
    raw = row.get("row", "{}")
    if pd.isna(raw):
        return {}
    try:
        return json.loads(str(raw))
    except Exception:
        return {"_raw": str(raw)}


def _find(source_tables: pd.DataFrame, table: str, key: str) -> Dict[str, Any]:
    m = source_tables[(source_tables["table"] == table) & (source_tables["key"] == key)]
    if m.empty:
        raise KeyError(f"Missing source table row: table={table}, key={key}")
    row = m.iloc[0]
    payload = _json_row(row)
    payload["_source"] = row.get("source", "")
    payload["_table"] = table
    payload["_key"] = key
    return payload


def _fmt(x: Any, nd: int = 6) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    source_tables = _read_csv(SOURCE_TABLES)
    relevant_files = _read_csv(RELEVANT_FILES)
    joint_patch = _read_csv(JOINT_PATCH)

    # Canonical table choices.
    act_planck = _find(source_tables, "ACT_refined_cases", "Planck_PR4_amp1")
    act_vanilla = _find(source_tables, "ACT_refined_cases", "A_vanilla_late_amp1")
    act_kprimary = _find(source_tables, "ACT_refined_cases", "Kprimary_amp1")
    act_kcmb8 = _find(source_tables, "ACT_refined_cases", "Kprimary_CMB8_amp")
    act_scan = _find(source_tables, "ACT_refined_scan_best", "Kprimary_scan_best")

    high_planck = _find(source_tables, "Planck_ns_As_extended_best", "Planck_PR4")
    high_strict = _find(source_tables, "Planck_ns_As_extended_best", "C_s_surv0")
    high_vanilla = _find(source_tables, "Planck_ns_As_best", "A_vanilla")
    high_zbest = _find(source_tables, "Planck_ns_As_extended_best", "C_Zbest")

    total_planck = _find(source_tables, "Planck_highl_lowl_tau_best", "Planck_PR4")
    total_strict = _find(source_tables, "Planck_highl_lowl_tau_best", "C_s_surv0")

    lens_planck = _find(source_tables, "Planck_lensing_official_formula", "Planck_PR4_12F")
    lens_strict = _find(source_tables, "Planck_lensing_official_formula", "BPB_C_s_surv0_12F")
    lens_vanilla = _find(source_tables, "Planck_lensing_official_formula", "BPB_vanilla_control")

    kern_planck = _find(source_tables, "Planck_lensing_kernel_models", "Planck_PR4")
    kern_strict = _find(source_tables, "Planck_lensing_kernel_models", "BPB_C_s_surv0")
    kern_vanilla = _find(source_tables, "Planck_lensing_kernel_models", "BPB_vanilla")

    j = {r["label"]: r.to_dict() for _, r in joint_patch.iterrows()}

    canonical_rows: List[Dict[str, Any]] = [
        {
            "paper_table": "T_CMB_ACT_DR6_lenslike",
            "status": "CANONICAL",
            "role": "Full ACT DR6 lenslike diagnostic table",
            "source": act_planck["_source"] + "; " + act_scan["_source"],
            "entries": "Planck_PR4_amp1; A_vanilla_late_amp1; Kprimary_amp1; Kprimary_CMB8_amp; Kprimary_scan_best",
            "key_numbers": (
                f"Planck chi2_ACT={_fmt(act_planck['chi2_ACT'])}; "
                f"vanilla chi2_ACT={_fmt(act_vanilla['chi2_ACT'])}; "
                f"Kprimary amp1 chi2_ACT={_fmt(act_kprimary['chi2_ACT'])}; "
                f"Kprimary CMB8 amp chi2_ACT={_fmt(act_kcmb8['chi2_ACT'])}; "
                f"scan best amp_PhiPhi={_fmt(act_scan['amp_PhiPhi'],3)}, "
                f"Sigma_eff={_fmt(act_scan['Sigma_auto_eff'],9)}, chi2_ACT={_fmt(act_scan['chi2_ACT'])}"
            ),
            "wording_level": "Full ACT DR6 lenslike diagnostic; not a standalone full cosmological MCMC.",
        },
        {
            "paper_table": "T_CMB_PlanckLite_high_ell",
            "status": "CANONICAL",
            "role": "Primary Planck-lite TTTEEE high-ell table",
            "source": high_planck["_source"] + "; " + high_vanilla["_source"],
            "entries": "Planck_PR4; C_s_surv0; C_Zbest; A_vanilla",
            "key_numbers": (
                f"Planck best_full_chi2={_fmt(high_planck['best_full_chi2'])}; "
                f"C_s_surv0 best_full_chi2={_fmt(high_strict['best_full_chi2'])}; "
                f"delta={_fmt(high_strict['delta_best_full_chi2_vs_best_Planck'])}; "
                f"C_Zbest best_full_chi2={_fmt(high_zbest['best_full_chi2'])}; "
                f"vanilla best_full_chi2={_fmt(high_vanilla['best_full_chi2'])}, "
                f"vanilla delta={_fmt(high_vanilla['delta_best_full_chi2_vs_best_Planck'])}"
            ),
            "wording_level": "Planck-lite/proxy high-ell TTTEEE; not official full Planck nuisance likelihood.",
        },
        {
            "paper_table": "T_CMB_PlanckLite_high_low_tau",
            "status": "CANONICAL_SUPPORTING",
            "role": "High-ell + native lowT/lowE + tau/As extension table",
            "source": total_planck["_source"],
            "entries": "Planck_PR4; C_s_surv0",
            "key_numbers": (
                f"Planck total={_fmt(total_planck['best_total_chi2'])}; "
                f"C_s_surv0 total={_fmt(total_strict['best_total_chi2'])}; "
                f"delta_total={_fmt(total_strict['delta_best_total_chi2_vs_Planck'])}; "
                f"delta_highl={_fmt(total_strict['delta_best_highl_chi2_vs_Planck'])}; "
                f"delta_lowT_lowE={_fmt(total_strict['delta_best_lowT_lowE_chi2_vs_Planck'])}"
            ),
            "wording_level": "Planck-lite plus low-ell proxy/supporting closure.",
        },
        {
            "paper_table": "T_CMB_Planck_lensing_CMBmarged",
            "status": "CANONICAL",
            "role": "Planck CMB-marginalized lensing official-formula diagnostic",
            "source": lens_planck["_source"],
            "entries": "Planck_PR4_12F; BPB_C_s_surv0_12F; BPB_vanilla_control",
            "key_numbers": (
                f"Planck best_chi2={_fmt(lens_planck['best_chi2_lensing'])}; "
                f"BPB strict best_amp_phi={_fmt(lens_strict['best_amp_phi'],3)}, "
                f"Sigma_eff={_fmt(lens_strict['best_Sigma_eff'],9)}, "
                f"best_chi2={_fmt(lens_strict['best_chi2_lensing'])}; "
                f"vanilla best_chi2={_fmt(lens_vanilla['best_chi2_lensing'])}"
            ),
            "wording_level": "Official Planck lensing CMBmarged formula diagnostic; scalar amplitude profiled.",
        },
        {
            "paper_table": "T_CMB_Planck_lensing_kernel_three_band",
            "status": "SUPPORTING_NOT_PRIMARY",
            "role": "Kernel-shape diagnostic beyond scalar amplitude",
            "source": kern_planck["_source"],
            "entries": "Planck_PR4; BPB_C_s_surv0; BPB_vanilla",
            "key_numbers": (
                f"BPB three-band chi2={_fmt(kern_strict['chi2'])}; "
                f"Planck three-band chi2={_fmt(kern_planck['chi2'])}; "
                f"vanilla three-band chi2={_fmt(kern_vanilla['chi2'])}"
            ),
            "wording_level": "Supporting forensic/kernel diagnostic; do not oversell as fitted physical kernel.",
        },
        {
            "paper_table": "T_CMB_Joint_ACT_Planck_scalar_lensing",
            "status": "CANONICAL_PATCHED_P3A1B",
            "role": "Joint scalar Weyl-auto amplitude closure table",
            "source": str(JOINT_PATCH) + "; outputs/cmb/cmb15a_v4_joint_act_planck_scalar_amp_map1782.grid.csv; outputs/cmb/cmb15a_v4_joint_act_planck_scalar_amp_map1782.refs.csv",
            "entries": "ACT-only best; Planck-only best; Joint scalar best; separate-minima sum; scalar penalty",
            "key_numbers": (
                f"ACT-only amp_phi={_fmt(j['act_best']['amp_phi'],3)}, chi2_ACT={_fmt(j['act_best']['chi2_ACT'])}; "
                f"Planck-only amp_phi={_fmt(j['planck_best']['amp_phi'],3)}, chi2_Planck_lensing={_fmt(j['planck_best']['chi2_Planck_lensing'])}; "
                f"Joint amp_phi={_fmt(j['joint_best']['amp_phi'],3)}, chi2_joint={_fmt(j['joint_best']['chi2_joint'])}; "
                f"separate_minima_sum={_fmt(j['separate_minima_sum']['chi2_joint'])}; "
                f"penalty={_fmt(j['joint_scalar_penalty']['chi2_joint'])}"
            ),
            "wording_level": "Scalar-amplitude diagnostic closure; not a full official combined ACT+Planck likelihood.",
        },
    ]

    # Existing image/figure inventory from A0 scan.
    suffixes = {".pdf", ".png", ".jpg", ".jpeg", ".svg"}
    fig_inventory = relevant_files[
        relevant_files["suffix"].str.lower().isin(suffixes)
        & relevant_files["path"].str.contains("cmb|CMB|paperIII|PaperIII", case=False, na=False)
    ].copy()
    if fig_inventory.empty:
        fig_inventory = pd.DataFrame(columns=["path", "suffix", "size_bytes", "mtime", "paper_use", "status"])
    else:
        def classify(path: str) -> tuple[str, str]:
            p = path.lower()
            if "tt_shape" in p:
                return "Optional proxy figure: TT shape sanity / face-minus shape non-crash", "EXISTING_PROXY_OPTIONAL"
            if "polarization" in p:
                return "Optional proxy figure: TT/EE/TE polarization sanity", "EXISTING_PROXY_OPTIONAL"
            if "lensing_weyl_proxy" in p:
                return "Optional proxy figure: CMB lensing/Weyl proxy bridge", "EXISTING_PROXY_OPTIONAL"
            if "peak_proxy" in p:
                return "Optional proxy figure: acoustic peak budget", "EXISTING_PROXY_OPTIONAL"
            return "Existing CMB diagnostic figure", "EXISTING_DIAGNOSTIC"
        classified = fig_inventory["path"].map(classify)
        fig_inventory["paper_use"] = [x[0] for x in classified]
        fig_inventory["status"] = [x[1] for x in classified]
        fig_inventory = fig_inventory.sort_values("path")

    figure_plan_rows = [
        {
            "figure_id": "Fig_P3_1_ACT_lensing_gate",
            "priority": "PRIMARY",
            "recommended_content": "ACT DR6 lenslike chi2 comparison: Planck_PR4, vanilla late-density, Kprimary amp1, Kprimary scan best.",
            "source_tables": "T_CMB_ACT_DR6_lenslike",
            "generation_status": "TO_GENERATE_FROM_TABLES",
            "reason": "This is the cleanest visual proof that vanilla crashes while face-minus/Kprimary survives ACT with mild positive Weyl-auto amplitude.",
        },
        {
            "figure_id": "Fig_P3_2_PlanckLite_high_ell",
            "priority": "PRIMARY",
            "recommended_content": "Planck-lite high-ell TTTEEE bar/waterfall: Planck_PR4, C_s_surv0, C_Zbest, vanilla.",
            "source_tables": "T_CMB_PlanckLite_high_ell; T_CMB_PlanckLite_high_low_tau",
            "generation_status": "TO_GENERATE_FROM_TABLES",
            "reason": "This locks the core paper claim: strict face-minus has no high-ell shape crash; vanilla is catastrophically excluded.",
        },
        {
            "figure_id": "Fig_P3_3_joint_lensing_scalar_closure",
            "priority": "PRIMARY",
            "recommended_content": "ACT-only, Planck-only and joint scalar amp_phi minima with joint penalty annotation.",
            "source_tables": "T_CMB_Joint_ACT_Planck_scalar_lensing",
            "generation_status": "TO_GENERATE_FROM_PATCHED_P3A1B_TABLE",
            "reason": "This visually supports response hierarchy and coherent scalar Weyl-auto closure.",
        },
        {
            "figure_id": "Fig_P3_4_optional_proxy_contact_sheet",
            "priority": "OPTIONAL_APPENDIX_OR_SUPPLEMENT",
            "recommended_content": "Existing TT/EE/TE, peak and lensing proxy PNGs as a compact appendix/contact-sheet if desired.",
            "source_tables": "Existing proxy PNG inventory",
            "generation_status": "OPTIONAL_EXISTING_FIGURES",
            "reason": "Useful audit trail, but too proxy-heavy for main paper unless space allows.",
        },
    ]

    open_items = [
        {
            "id": "P3A2.O01",
            "item": "Generate publication-grade figures from canonical tables",
            "status": "OPEN_NEXT",
            "next_gate": "P3-A3 or P3-FIG1",
            "note": "Recommended first figures: ACT gate, Planck-lite high-ell waterfall, joint scalar closure.",
        },
        {
            "id": "P3A2.O02",
            "item": "Patch source_tables joint summary or explicitly supersede it",
            "status": "OPEN_LOW_RISK",
            "next_gate": "P3-A3 wording/source-lock cleanup",
            "note": "paperIII_cmb_source_tables.csv still contains old Joint_ACT_Planck_v4_summary with joint_best=None; P3-A1b patch is authoritative.",
        },
        {
            "id": "P3A2.O03",
            "item": "Freeze wording levels",
            "status": "OPEN_NEXT",
            "next_gate": "P3-A3 wording table",
            "note": "ACT full lenslike diagnostic; Planck-lite proxy/lite; Planck lensing official-formula diagnostic; joint lensing scalar-amplitude diagnostic closure.",
        },
    ]

    canonical_df = pd.DataFrame(canonical_rows)
    fig_plan_df = pd.DataFrame(figure_plan_rows)
    open_df = pd.DataFrame(open_items)

    canonical_csv = OUTDIR / "paperIII_cmb_p3a2_canonical_tables.csv"
    fig_inventory_csv = OUTDIR / "paperIII_cmb_p3a2_figure_inventory.csv"
    fig_plan_csv = OUTDIR / "paperIII_cmb_p3a2_figure_plan.csv"
    open_csv = OUTDIR / "paperIII_cmb_p3a2_open_items.csv"
    json_out = OUTDIR / "paperIII_cmb_p3a2_canonical_tables_fig_inventory.json"
    md_out = OUTDIR / "paperIII_cmb_p3a2_canonical_tables_fig_inventory.md"

    canonical_df.to_csv(canonical_csv, index=False)
    fig_inventory.to_csv(fig_inventory_csv, index=False)
    fig_plan_df.to_csv(fig_plan_csv, index=False)
    open_df.to_csv(open_csv, index=False)

    payload = {
        "status": STATUS,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "source_tables": str(SOURCE_TABLES),
            "relevant_files": str(RELEVANT_FILES),
            "joint_patch": str(JOINT_PATCH),
        },
        "safety": {
            "likelihood_run": "NO",
            "NS_touched": "NO",
            "write_scope": "outputs/cmb_paperIII only",
        },
        "canonical_tables": canonical_rows,
        "figure_plan": figure_plan_rows,
        "n_existing_cmb_figures_found": int(len(fig_inventory)),
        "open_items": open_items,
        "gate_decisions": {
            "P3-A2 canonical table choice": "PASS",
            "P3-A2 figure inventory": "PASS",
            "P3-A1b joint patch incorporated": "PASS",
            "claims paper-ready": "CANONICAL_TABLES_SELECTED_FIGURES_PENDING",
            "next_gate": "P3-A3 wording lock or figure generation",
        },
    }
    json_out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    lines: List[str] = []
    lines.append("# Paper III CMB P3-A2 — canonical table choice and figure inventory")
    lines.append("")
    lines.append("## Status")
    lines.append("")
    lines.append(f"**{STATUS}**")
    lines.append("")
    lines.append("Safety: read-only on `outputs/cmb` and `outputs/cmb_paperIII` inputs; writes only to `outputs/cmb_paperIII`; no ACT/Planck likelihood run; no NS touch.")
    lines.append("")
    lines.append("## Canonical paper tables")
    lines.append("")
    lines.append("| paper table | status | role | key numbers | wording level |")
    lines.append("|---|---|---|---|---|")
    for row in canonical_rows:
        lines.append(f"| `{row['paper_table']}` | {row['status']} | {row['role']} | {row['key_numbers']} | {row['wording_level']} |")
    lines.append("")
    lines.append("## Recommended figure plan")
    lines.append("")
    lines.append("| figure id | priority | content | generation status |")
    lines.append("|---|---|---|---|")
    for row in figure_plan_rows:
        lines.append(f"| `{row['figure_id']}` | {row['priority']} | {row['recommended_content']} | {row['generation_status']} |")
    lines.append("")
    lines.append(f"Existing CMB/proxy image files found by A0 inventory: **{len(fig_inventory)}**")
    lines.append("")
    if len(fig_inventory):
        lines.append("### Existing proxy/diagnostic figures")
        lines.append("")
        lines.append("| path | status | paper use |")
        lines.append("|---|---|---|")
        for _, row in fig_inventory.iterrows():
            lines.append(f"| `{row['path']}` | {row['status']} | {row['paper_use']} |")
        lines.append("")
    lines.append("## Gate decisions")
    lines.append("")
    for k, v in payload["gate_decisions"].items():
        lines.append(f"- **{k}:** {v}")
    lines.append("")
    lines.append("## Open items")
    lines.append("")
    lines.append("| id | item | status | note |")
    lines.append("|---|---|---|---|")
    for row in open_items:
        lines.append(f"| {row['id']} | {row['item']} | {row['status']} | {row['note']} |")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    for p in [canonical_csv, fig_inventory_csv, fig_plan_csv, open_csv, json_out, md_out]:
        lines.append(f"- `{p}`")
    lines.append("")
    md_out.write_text("\n".join(lines), encoding="utf-8")

    print(f"Saved: {md_out}")
    print(f"Saved: {json_out}")
    print(f"Saved: {canonical_csv}")
    print(f"Saved: {fig_inventory_csv}")
    print(f"Saved: {fig_plan_csv}")
    print(f"Saved: {open_csv}")
    print(STATUS)


if __name__ == "__main__":
    main()
