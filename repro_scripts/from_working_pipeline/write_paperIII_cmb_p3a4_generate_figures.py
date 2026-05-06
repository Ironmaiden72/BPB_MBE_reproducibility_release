#!/usr/bin/env python3
"""
Paper III CMB P3-A4 — canonical figure generation.

Non-destructive gate:
- Reads only audited CMB artefacts already present in outputs/cmb_paperIII.
- Does not run ACT, Planck, CAMB, Cobaya, clik, or any likelihood.
- Does not touch outputs/ns, src/ns, or checkpoints.
- Writes only under outputs/cmb_paperIII.

Generated primary figures:
- Fig_P3_1_ACT_lensing_gate
- Fig_P3_2_PlanckLite_high_ell
- Fig_P3_3_joint_lensing_scalar_closure
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATUS = "P3_A4_CANONICAL_FIGURES_GENERATED"

# Audited values fixed by P3-A1/P3-A1b/P3-A2/P3-A3.
ACT_VALUES = [
    ("Planck_PR4\namp=1", 13.951140258982589),
    ("Kprimary\namp=1", 17.470020679940834),
    ("Kprimary scan\nbest", 14.15359945719008),
    ("Vanilla late\ndensity", 227.5207684179947),
]
ACT_SCAN_AMP = 1.043
ACT_SCAN_SIGMA = 1.0212737145349429

PLANCKLITE_VALUES = [
    ("C_s_surv0\nstrict", 749.191767798833),
    ("Planck_PR4\nsanity", 765.8691646821364),
    ("C_Zbest", 777.286828),
    ("Vanilla\ncontrol", 8851.940268986189),
]
PLANCKLITE_DELTA_STRICT = -16.67739688330346
PLANCKLITE_VANILLA_DELTA = 8086.071104304053

DEFAULT_JOINT_ROWS = [
    {
        "label": "act_best",
        "amp_phi": 1.043,
        "Sigma_eff": 1.0212737145349429,
        "chi2_ACT": 14.15359945719008,
        "chi2_Planck_lensing": 10.198317000740625,
        "chi2_joint": 24.351916457930702,
    },
    {
        "label": "planck_best",
        "amp_phi": 1.0750000000000002,
        "Sigma_eff": 1.0368220676663862,
        "chi2_ACT": 15.90870191352052,
        "chi2_Planck_lensing": 8.742460151124908,
        "chi2_joint": 24.651162064645423,
    },
    {
        "label": "joint_best",
        "amp_phi": 1.057,
        "Sigma_eff": 1.02810505299799,
        "chi2_ACT": 14.478027412095692,
        "chi2_Planck_lensing": 9.243264395408051,
        "chi2_joint": 23.721291807503743,
    },
    {
        "label": "separate_minima_sum",
        "amp_phi": np.nan,
        "Sigma_eff": np.nan,
        "chi2_ACT": 14.15359945719008,
        "chi2_Planck_lensing": 8.742460151124908,
        "chi2_joint": 22.896059608314985,
    },
    {
        "label": "joint_scalar_penalty",
        "amp_phi": np.nan,
        "Sigma_eff": np.nan,
        "chi2_ACT": np.nan,
        "chi2_Planck_lensing": np.nan,
        "chi2_joint": 0.8252321991887577,
    },
]

FIG_CAPTIONS = {
    "Fig_P3_1_ACT_lensing_gate": "ACT DR6 lenslike CMB dictionary gate. The vanilla late-density CMB mapping is catastrophically disfavoured, while the strict face-minus/Kprimary dictionary remains close to the Planck sanity point and improves under a small positive Weyl-auto amplitude profile.",
    "Fig_P3_2_PlanckLite_high_ell": "Planck-lite high-ell TTTEEE proxy. The strict C_s_surv0 dictionary improves the profiled high-ell chi2 relative to the Planck sanity point in this lite diagnostic, while the vanilla dictionary fails by thousands in chi2.",
    "Fig_P3_3_joint_lensing_scalar_closure": "Scalar-amplitude closure for ACT DR6 and Planck CMB-marginalized lensing. ACT-only and Planck-only windows prefer nearby positive amplitudes, and one common amp_phi=1.057 gives a joint penalty of only Delta chi2=0.825 relative to separate minima.",
}


def repo_root() -> Path:
    return Path.cwd()


def out_paths(root: Path) -> Tuple[Path, Path]:
    out = root / "outputs" / "cmb_paperIII"
    figdir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    return out, figdir


def load_joint_patch(out: Path) -> Tuple[pd.DataFrame, str]:
    p = out / "paperIII_cmb_joint_lens_v4_patch.csv"
    if p.exists():
        df = pd.read_csv(p)
        required = {"label", "amp_phi", "Sigma_eff", "chi2_ACT", "chi2_Planck_lensing", "chi2_joint"}
        missing = sorted(required - set(df.columns))
        if missing:
            raise RuntimeError(f"Joint patch exists but is missing columns: {missing}")
        return df, str(p)
    return pd.DataFrame(DEFAULT_JOINT_ROWS), "embedded audited P3-A1b values"


def load_captions(out: Path) -> Tuple[Dict[str, str], str]:
    p = out / "paperIII_cmb_p3a3_figure_captions.csv"
    if not p.exists():
        return dict(FIG_CAPTIONS), "embedded P3-A3 captions"
    df = pd.read_csv(p)
    caps = dict(FIG_CAPTIONS)
    if {"figure_id", "caption"}.issubset(df.columns):
        for _, row in df.iterrows():
            if isinstance(row.get("figure_id"), str) and isinstance(row.get("caption"), str):
                caps[row["figure_id"]] = row["caption"]
    return caps, str(p)


def save_figure(fig: plt.Figure, figdir: Path, stem: str) -> Dict[str, str]:
    pdf = figdir / f"{stem}.pdf"
    png = figdir / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return {"pdf": str(pdf), "png": str(png), "pdf_size": pdf.stat().st_size, "png_size": png.stat().st_size}


def annotate_bars(ax, bars, fmt="{:.3f}", xlog: bool = False):
    for bar in bars:
        val = bar.get_width() if hasattr(bar, "get_width") else bar.get_height()
        if hasattr(bar, "get_width"):
            y = bar.get_y() + bar.get_height() / 2
            if xlog:
                x = val * 1.05
            else:
                x = val + max(1.0, val * 0.01)
            ax.text(x, y, fmt.format(val), va="center", ha="left", fontsize=9)
        else:
            x = bar.get_x() + bar.get_width() / 2
            y = val + max(0.2, val * 0.01)
            ax.text(x, y, fmt.format(val), va="bottom", ha="center", fontsize=9)


def fig_act(figdir: Path) -> Dict[str, str]:
    labels = [x[0] for x in ACT_VALUES]
    vals = [x[1] for x in ACT_VALUES]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    bars = ax.barh(y, vals)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("ACT lenslike chi2 (log scale)")
    ax.set_title("Fig. P3-1 - ACT DR6 lenslike dictionary gate")
    ax.grid(True, axis="x", alpha=0.25)
    annotate_bars(ax, bars, xlog=True)
    ax.text(
        0.02, 0.02,
        f"Kprimary scan best: amp_phi={ACT_SCAN_AMP:.3f}, Sigma_eff={ACT_SCAN_SIGMA:.6f}\nVanilla control is the failed CMB dictionary.",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    return save_figure(fig, figdir, "fig_P3_1_ACT_lensing_gate")


def fig_plancklite(figdir: Path) -> Dict[str, str]:
    labels = [x[0] for x in PLANCKLITE_VALUES]
    vals = [x[1] for x in PLANCKLITE_VALUES]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    bars = ax.barh(y, vals)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Planck-lite high-ell TTTEEE chi2 (log scale)")
    ax.set_title("Fig. P3-2 - Planck-lite high-ell dictionary gate")
    ax.grid(True, axis="x", alpha=0.25)
    annotate_bars(ax, bars, xlog=True)
    ax.text(
        0.02, 0.02,
        f"Strict C_s_surv0 vs Planck sanity: Delta chi2={PLANCKLITE_DELTA_STRICT:.6f}\nVanilla control: Delta chi2=+{PLANCKLITE_VANILLA_DELTA:.3f}",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    return save_figure(fig, figdir, "fig_P3_2_PlanckLite_high_ell")


def row_by_label(df: pd.DataFrame, label: str) -> pd.Series:
    hit = df[df["label"] == label]
    if hit.empty:
        raise RuntimeError(f"Missing joint patch row: {label}")
    return hit.iloc[0]


def fig_joint(figdir: Path, joint_df: pd.DataFrame) -> Dict[str, str]:
    rows = [row_by_label(joint_df, x) for x in ["act_best", "planck_best", "joint_best"]]
    labels = ["ACT-only\nmin", "Planck-only\nmin", "Joint scalar\nmin"]
    joint_vals = [float(r["chi2_joint"]) for r in rows]
    act_vals = [float(r["chi2_ACT"]) for r in rows]
    planck_vals = [float(r["chi2_Planck_lensing"]) for r in rows]
    amps = [float(r["amp_phi"]) for r in rows]
    sigmas = [float(r["Sigma_eff"]) for r in rows]

    sep = float(row_by_label(joint_df, "separate_minima_sum")["chi2_joint"])
    penalty = float(row_by_label(joint_df, "joint_scalar_penalty")["chi2_joint"])

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    width = 0.62
    bars_act = ax.bar(x, act_vals, width=width, label="ACT contribution")
    bars_planck = ax.bar(x, planck_vals, width=width, bottom=act_vals, label="Planck lensing contribution")
    ax.axhline(sep, linestyle="--", linewidth=1.2, label=f"Separate minima sum = {sep:.6f}")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Joint chi2 = chi2_ACT + chi2_Planck_lensing")
    ax.set_title("Fig. P3-3 - Joint ACT+Planck scalar lensing closure", pad=14)
    ax.grid(True, axis="y", alpha=0.25)
    ax.set_ylim(0, max(joint_vals) + 4.1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=8, frameon=True)
    for i, (total, amp, sigma) in enumerate(zip(joint_vals, amps, sigmas)):
        ax.text(i, total + 0.28, f"{total:.3f}\namp={amp:.3f}\nSigma={sigma:.6f}", ha="center", va="bottom", fontsize=7.6)
    ax.text(
        0.02, 0.02,
        f"Joint scalar penalty relative to separate minima: Delta chi2={penalty:.6f}",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    return save_figure(fig, figdir, "fig_P3_3_joint_lensing_scalar_closure")


def write_manifest(out: Path, figdir: Path, fig_outputs: Dict[str, Dict[str, str]], joint_source: str, caption_source: str, captions: Dict[str, str]):
    timestamp = datetime.now(timezone.utc).isoformat()
    rows = []
    for fig_id, paths in fig_outputs.items():
        rows.append({
            "figure_id": fig_id,
            "status": "GENERATED_PRIMARY",
            "pdf": paths["pdf"],
            "png": paths["png"],
            "pdf_size_bytes": paths["pdf_size"],
            "png_size_bytes": paths["png_size"],
            "caption": captions.get(fig_id, ""),
        })
    fig_csv = out / "paperIII_cmb_p3a4_figure_outputs.csv"
    pd.DataFrame(rows).to_csv(fig_csv, index=False)

    checks = [
        {"check": "status", "value": STATUS, "pass": True},
        {"check": "figure_count_primary", "value": len(fig_outputs), "pass": len(fig_outputs) == 3},
        {"check": "joint_patch_source", "value": joint_source, "pass": True},
        {"check": "caption_source", "value": caption_source, "pass": True},
        {"check": "likelihood_run", "value": "NO", "pass": True},
        {"check": "NS_touched", "value": "NO", "pass": True},
        {"check": "write_scope", "value": "outputs/cmb_paperIII only", "pass": True},
    ]
    checks_csv = out / "paperIII_cmb_p3a4_figure_checks.csv"
    pd.DataFrame(checks).to_csv(checks_csv, index=False)

    payload = {
        "status": STATUS,
        "timestamp_utc": timestamp,
        "safety": {
            "likelihood_run": "NO",
            "NS_touched": "NO",
            "read_scope": ["outputs/cmb_paperIII audited artefacts", "embedded audited P3 values if needed"],
            "write_scope": "outputs/cmb_paperIII only",
        },
        "inputs": {
            "joint_patch": joint_source,
            "captions": caption_source,
        },
        "figures": rows,
        "next_gate": "P3-A5 LaTeX figure insertion / table rendering",
    }
    json_path = out / "paperIII_cmb_p3a4_figure_generation.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_path = out / "paperIII_cmb_p3a4_figure_generation.md"
    lines = [
        "# Paper III CMB P3-A4 - canonical figure generation",
        "",
        "## Status",
        "",
        f"**{STATUS}**",
        "",
        "Safety: no ACT/Planck likelihood run; no NS touch; writes only to `outputs/cmb_paperIII`.",
        "",
        "## Generated primary figures",
        "",
        "| figure | PDF | PNG | caption lock |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| `{row['figure_id']}` | `{row['pdf']}` | `{row['png']}` | {row['caption']} |")
    lines += [
        "",
        "## Gate decisions",
        "",
        "| decision | value |",
        "|---|---|",
        "| P3-A4 figure generation | PASS |",
        "| Primary figures generated | 3/3 |",
        "| P3-A1b joint patch used | YES |",
        "| NS touched | NO |",
        "| likelihood run | NO |",
        "| next gate | P3-A5 LaTeX figure insertion / table rendering |",
        "",
        "## Outputs",
        "",
        f"- `{fig_csv}`",
        f"- `{checks_csv}`",
        f"- `{json_path}`",
        f"- `{md_path}`",
    ]
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return md_path, json_path, fig_csv, checks_csv


def main() -> None:
    root = repo_root()
    out, figdir = out_paths(root)
    joint_df, joint_source = load_joint_patch(out)
    captions, caption_source = load_captions(out)

    fig_outputs: Dict[str, Dict[str, str]] = {}
    fig_outputs["Fig_P3_1_ACT_lensing_gate"] = fig_act(figdir)
    fig_outputs["Fig_P3_2_PlanckLite_high_ell"] = fig_plancklite(figdir)
    fig_outputs["Fig_P3_3_joint_lensing_scalar_closure"] = fig_joint(figdir, joint_df)

    md_path, json_path, fig_csv, checks_csv = write_manifest(out, figdir, fig_outputs, joint_source, caption_source, captions)
    print(f"Status: {STATUS}")
    print(f"Saved: {md_path}")
    print(f"Saved: {json_path}")
    print(f"Saved: {fig_csv}")
    print(f"Saved: {checks_csv}")
    for fig_id, paths in fig_outputs.items():
        print(f"Saved: {fig_id}: {paths['pdf']} ; {paths['png']}")


if __name__ == "__main__":
    main()
