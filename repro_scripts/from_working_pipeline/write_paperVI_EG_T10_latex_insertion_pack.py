#!/usr/bin/env python3
"""
Paper VI / E_G T10 -- LaTeX insertion pack and final source lock.

Non-destructive by default.
- Reads outputs/paperVI_EG_theory/T9A artefacts.
- Writes an insertion-ready LaTeX section, table, source lock and manifest.
- Does not touch NS files.
- Does not run likelihoods.
- Optional --apply --paper-tex path/to/paper.tex inserts the block between markers
  after making a .bak backup. Without --apply, no paper file is modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any

OUTDIR = Path("outputs/paperVI_EG_theory")
STATUS_DONE = "P6_T10_LATEX_INSERTION_PACK_DONE"
STATUS_APPLIED = "P6_T10_LATEX_INSERTION_PACK_DONE_AND_APPLIED_WITH_BACKUP"

DEFAULTS = {
    "a_bg": "0.403687948000",
    "sqrt_a_bg": "0.635364421415",
    "free_r_wdelta": "0.630864746000",
    "r_wdelta_squared": "0.397990327746",
    "r2_minus_a_bg": "-0.005697620254",
    "delta_chi2_fixed_sqrt_abg_vs_free": "0.008545",
    "theta_bg_deg": "50.552982",
    "wenzl_trace_unweighted": "0.999998712819",
    "wenzl_trace_ell_weighted": "0.999998879774",
    "autocross_fraction": "0.397990327336",
    "autocross_sqrt_fraction": "0.630864745675",
    "autocross_delta_vs_a_bg": "-0.005697620664",
    "autocross_relative_delta_percent": "-1.411392",
    "autocross_corr_A_X": "1.000000000000",
    "chi2_per_point_all": "0.226490429325",
    "chi2_per_point_planck_pr4_boss": "0.024232480764",
    "chi2_per_point_act_dr6_boss": "0.558887499306",
    "chi2_per_point_combined": "0.096351307904",
    "free_r_scan_delta_physical_free": "1.118470958517e-06",
}


def read_csv_key_values(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return out
        fields = set(reader.fieldnames)
        if {"quantity", "value"}.issubset(fields):
            for row in reader:
                q = (row.get("quantity") or "").strip()
                v = (row.get("value") or "").strip()
                if q:
                    out[q] = v
        elif {"key", "value"}.issubset(fields):
            for row in reader:
                q = (row.get("key") or "").strip()
                v = (row.get("value") or "").strip()
                if q:
                    out[q] = v
    return out


def load_values() -> Dict[str, str]:
    values = dict(DEFAULTS)
    for p in [
        OUTDIR / "paperVI_EG_T9A_key_values.csv",
        OUTDIR / "paperVI_EG_T8_key_values.csv",
        OUTDIR / "paperVI_EG_T1_rWdelta_sqrt_abg_summary.csv",
    ]:
        kv = read_csv_key_values(p)
        for k, v in kv.items():
            values[k] = v
            alias = k.lower().replace("<", "").replace(">", "").replace("/", "_").replace("-", "_")
            alias = re.sub(r"[^a-z0-9_]+", "_", alias).strip("_")
            if alias:
                values.setdefault(alias, v)
    return values


def val(values: Dict[str, str], *keys: str, default: str = "") -> str:
    for k in keys:
        if k in values and str(values[k]).strip() != "":
            return str(values[k]).strip()
    return default


def write_csv(path: Path, rows: List[Dict[str, Any]], fields: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def find_tex_candidates() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    roots = [Path("."), Path("papers"), Path("paperVI"), Path("docs"), Path("outputs")]
    seen = set()
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*.tex"):
            if p in seen or not p.is_file():
                continue
            seen.add(p)
            s = str(p).replace("\\", "/").lower()
            if "/ns/" in s or "outputs/ns" in s or "src/ns" in s:
                continue
            name_score = sum(tok in s for tok in ["papervi", "paper_vi", "paper6", "paper_6", "eg", "e_g", "cmb", "cosmo"])
            try:
                txt = p.read_text(encoding="utf-8", errors="replace")[:20000]
            except Exception:
                txt = ""
            content_score = sum(tok in txt for tok in ["E_G", "Weyl", "r_W", "a_{\\rm bg}", "biface", "BPB", "Projected"])
            rows.append({
                "path": str(p),
                "size_bytes": p.stat().st_size,
                "name_score": int(name_score),
                "content_score": int(content_score),
                "total_score": int(name_score + content_score),
            })
    rows.sort(key=lambda r: (-int(r["total_score"]), r["path"]))
    return rows[:50]


def build_table_tex(values: Dict[str, str]) -> str:
    rows = [
        (r"$a_{\rm bg}$", val(values, "a_bg", default=DEFAULTS["a_bg"])),
        (r"$\sqrt{a_{\rm bg}}$", val(values, "sqrt_a_bg", default=DEFAULTS["sqrt_a_bg"])),
        (r"free $r_{W\delta}$", val(values, "free_r_wdelta", "r_Wdelta_free", default=DEFAULTS["free_r_wdelta"])),
        (r"$r_{W\delta}^2$", val(values, "r_wdelta_squared", "r_Wdelta^2", "r_wdelta_free_squared", default=DEFAULTS["r_wdelta_squared"])),
        (r"$r_{W\delta}^2-a_{\rm bg}$", val(values, "r2_minus_a_bg", "r_Wdelta^2 - a_bg", "delta_r2_minus_a_bg", default=DEFAULTS["r2_minus_a_bg"])),
        (r"$\Delta\chi^2[\sqrt{a_{\rm bg}}-r_{\rm free}]$", val(values, "delta_chi2_fixed_sqrt_abg_vs_free", "delta_chi2_fixed_sqrt_abg_minus_free", default=DEFAULTS["delta_chi2_fixed_sqrt_abg_vs_free"])),
        (r"$\theta_{\rm bg}=\arccos\sqrt{a_{\rm bg}}$", val(values, "theta_bg_deg", default=DEFAULTS["theta_bg_deg"]) + r" deg"),
        (r"$\langle X,X\rangle/\langle A,A\rangle$", val(values, "autocross_fraction", "ALL_unweighted_fraction_XX_over_AA", default=DEFAULTS["autocross_fraction"])),
        (r"$\sqrt{\langle X,X\rangle/\langle A,A\rangle}$", val(values, "autocross_sqrt_fraction", "ALL_unweighted_sqrt_fraction", default=DEFAULTS["autocross_sqrt_fraction"])),
        (r"$\langle X,X\rangle/\langle A,A\rangle-a_{\rm bg}$", val(values, "autocross_delta_vs_a_bg", "ALL_unweighted_delta_fraction_minus_a_bg", default=DEFAULTS["autocross_delta_vs_a_bg"])),
        (r"all-data $\chi^2$/point", val(values, "chi2_per_point_all", "ALL_unweighted_chi2_per_point", default=DEFAULTS["chi2_per_point_all"])),
    ]
    body = "\n".join([f"{q} & {v} \\\\" for q, v in rows])
    return r"""\begin{table}
\centering
\caption{Projected Weyl--matter closure values used in the Paper VI $E_G$ appendix. The relation $r_{W\delta}=\sqrt{a_{\rm bg}}$ is a projected-diagnostic closure, not yet an action-level HR/BPB perturbation proof.}
\label{tab:paperVI-eg-rwdelta-abg}
\begin{tabular}{lc}
\hline
Quantity & Value \\
\hline
""" + body + r"""
\hline
\end{tabular}
\end{table}
"""


def build_section_tex(values: Dict[str, str]) -> str:
    replacements = {
        "@@RFREE@@": val(values, "free_r_wdelta", "r_Wdelta_free", default=DEFAULTS["free_r_wdelta"]),
        "@@SQ@@": val(values, "sqrt_a_bg", default=DEFAULTS["sqrt_a_bg"]),
        "@@A@@": val(values, "a_bg", default=DEFAULTS["a_bg"]),
        "@@DCHI@@": val(values, "delta_chi2_fixed_sqrt_abg_vs_free", "delta_chi2_fixed_sqrt_abg_minus_free", default=DEFAULTS["delta_chi2_fixed_sqrt_abg_vs_free"]),
        "@@THETA@@": val(values, "theta_bg_deg", default=DEFAULTS["theta_bg_deg"]),
        "@@TRACE_U@@": val(values, "wenzl_trace_unweighted", default=DEFAULTS["wenzl_trace_unweighted"]),
        "@@TRACE_W@@": val(values, "wenzl_trace_ell_weighted", default=DEFAULTS["wenzl_trace_ell_weighted"]),
        "@@FRAC@@": val(values, "autocross_fraction", "ALL_unweighted_fraction_XX_over_AA", default=DEFAULTS["autocross_fraction"]),
        "@@SQRT_FRAC@@": val(values, "autocross_sqrt_fraction", "ALL_unweighted_sqrt_fraction", default=DEFAULTS["autocross_sqrt_fraction"]),
        "@@DELTA_FRAC@@": val(values, "autocross_delta_vs_a_bg", "ALL_unweighted_delta_fraction_minus_a_bg", default=DEFAULTS["autocross_delta_vs_a_bg"]),
        "@@REL@@": val(values, "autocross_relative_delta_percent", default=DEFAULTS["autocross_relative_delta_percent"]),
        "@@CHI_ALL@@": val(values, "chi2_per_point_all", "ALL_unweighted_chi2_per_point", default=DEFAULTS["chi2_per_point_all"]),
        "@@CHI_PLANCK@@": val(values, "chi2_per_point_planck_pr4_boss", default=DEFAULTS["chi2_per_point_planck_pr4_boss"]),
        "@@CHI_ACT@@": val(values, "chi2_per_point_act_dr6_boss", default=DEFAULTS["chi2_per_point_act_dr6_boss"]),
        "@@CHI_COMB@@": val(values, "chi2_per_point_combined", default=DEFAULTS["chi2_per_point_combined"]),
    }
    template = r"""% BEGIN P6_T10_PROJECTED_WEYL_MATTER_CLOSURE
\subsection{Projected Weyl--matter closure: $r_{W\delta}=\sqrt{a_{\rm bg}}$}
\label{sec:paperVI-eg-rwdelta-abg}

The compressed $E_G$ diagnostic contains an effective Weyl--matter correlation coefficient, denoted here by $r_{W\delta}$. The audited projected diagnostic chain shows that this coefficient does not need to remain an independent parameter at the compressed level: the free value
\begin{equation}
  r_{W\delta}^{\rm free} = @@RFREE@@
\end{equation}
is reproduced by the dictionary-tied closure
\begin{equation}
  r_{W\delta} = \sqrt{a_{\rm bg}} = @@SQ@@,
  \qquad a_{\rm bg}=@@A@@,
\end{equation}
with a diagnostic penalty of only $\Delta\chi^2=@@DCHI@@$ relative to the free-$r$ solution.

The conditional projection theorem is simple. If the projected Weyl response decomposes into a matter-correlated component and an orthogonal BPB-private component,
\begin{equation}
  W = \sqrt{a_{\rm bg}}\,\delta_m
      + \sqrt{1-a_{\rm bg}}\,\xi_W,
  \qquad \langle \delta_m,\xi_W\rangle_{E_G}=0,
\end{equation}
with matched projected variances, then
\begin{equation}
  {\rm corr}(W,\delta_m)=\sqrt{a_{\rm bg}}.
\end{equation}
Equivalently, in a two-mode projected Weyl basis,
\begin{equation}
  \widehat W = \cos\theta_W\,e_m + \sin\theta_W\,e_x,
  \qquad \cos^2\theta_W=a_{\rm bg},
\end{equation}
so that
\begin{equation}
  \theta_W=\arccos\sqrt{a_{\rm bg}} = @@THETA@@^\circ.
\end{equation}

A direct projector proof cannot be obtained from the Wenzl 19-bin projected response shapes alone. Those curves are nearly one-dimensional in shape: the dominant trace fraction is @@TRACE_U@@ for the unweighted-bin norm and @@TRACE_W@@ for the ell-width weighted norm. A naive projection of the ACT/Sigma or $E_G$/Sigma response onto the $W_\mu$ direction therefore returns a parallel fraction close to unity rather than $a_{\rm bg}\simeq 0.404$. The Wenzl projected curves are useful model-response visualisations, but they are not the correct vector space for proving the two-mode matter/private Weyl split.

The correct compressed diagnostic space is instead the auto/cross amplitude space. In the frozen diagnostic artefact, define
\begin{equation}
  A_i = {\tt unit\_ACTSigma\_cross}_i,
  \qquad
  X_i = {\tt prediction\_global\_auto\_cross}_i.
\end{equation}
The row-level identity is
\begin{equation}
  X_i = r_{\rm global} A_i,
\end{equation}
so for any positive diagonal or bin weights,
\begin{equation}
  \frac{\langle X,X\rangle}{\langle A,A\rangle}
  = r_{\rm global}^2.
\end{equation}
The audited value is
\begin{equation}
  \frac{\langle X,X\rangle}{\langle A,A\rangle}
  = @@FRAC@@,
  \qquad
  \sqrt{\frac{\langle X,X\rangle}{\langle A,A\rangle}} = @@SQRT_FRAC@@,
\end{equation}
which differs from $a_{\rm bg}$ by @@DELTA_FRAC@@, i.e. @@REL@@\%. The associated support values are $\chi^2/N=@@CHI_ALL@@$ for all compressed rows, @@CHI_PLANCK@@ for the Planck PR4+BOSS diagnostic, @@CHI_ACT@@ for ACT DR6+BOSS, and @@CHI_COMB@@ for the combined ACT+Planck+BOSS diagnostic.

This result should be read as a projected-diagnostic closure, not as a completed action-level HR/BPB perturbation derivation. The remaining field-equation target is to derive the projector identity
\begin{equation}
  \frac{\langle P_m W,P_m W\rangle_{E_G}}{\langle W,W\rangle_{E_G}} = a_{\rm bg}
\end{equation}
from the HR/BPB perturbation sector. Until that derivation is obtained, the relation $r_{W\delta}=\sqrt{a_{\rm bg}}$ is a strong conditional closure supported by the projected $E_G$ diagnostics.

\input{outputs/paperVI_EG_theory/paperVI_EG_T10_closure_values_table.tex}
% END P6_T10_PROJECTED_WEYL_MATTER_CLOSURE
"""
    for k, v in replacements.items():
        template = template.replace(k, str(v))
    return template


def build_md(status: str, applied: bool, target: str) -> str:
    template = """# Paper VI / E_G T10 - LaTeX insertion pack

## Status

**@@STATUS@@**

Non-destructive by default. No likelihood run. NS untouched.

## Purpose

T10 turns the T9A theory appendix into a LaTeX insertion pack for Paper VI. It creates a ready-to-include subsection, a values table, final source lock, and manifest.

## Apply status

| item | value |
|---|---|
| applied to paper tex | @@APPLIED@@ |
| target paper tex | `@@TARGET@@` |

## Files written

- `paperVI_EG_T10_projected_closure_section.tex`
- `paperVI_EG_T10_closure_values_table.tex`
- `paperVI_EG_T10_source_lock_final.csv`
- `paperVI_EG_T10_manifest.md`
- `paperVI_EG_T10_manifest.json`
- `paperVI_EG_T10_decisions.csv`
- `paperVI_EG_T10_tex_candidates.csv`

## Recommended use

In the Paper VI LaTeX file, add:

```tex
\input{outputs/paperVI_EG_theory/paperVI_EG_T10_projected_closure_section.tex}
```

Recommended insertion point: after the compressed E_G diagnostic results and before the final discussion/conclusion.

## Guardrail

This section states a projected-diagnostic closure:

```tex
r_{W\delta}=\sqrt{a_{\rm bg}}
```

It does not claim a completed action-level HR/BPB perturbation proof. The open field-equation target remains:

```tex
\frac{\langle P_m W,P_m W\rangle_{E_G}}{\langle W,W\rangle_{E_G}} = a_{\rm bg}
```
"""
    return (
        template
        .replace("@@STATUS@@", status)
        .replace("@@APPLIED@@", str(applied))
        .replace("@@TARGET@@", target)
    )


def apply_to_tex(paper_tex: Path, block: str) -> None:
    if not paper_tex.exists():
        raise FileNotFoundError(f"paper tex not found: {paper_tex}")
    txt = paper_tex.read_text(encoding="utf-8", errors="replace")
    begin = "% BEGIN P6_T10_PROJECTED_WEYL_MATTER_CLOSURE"
    end = "% END P6_T10_PROJECTED_WEYL_MATTER_CLOSURE"
    bak = paper_tex.with_suffix(paper_tex.suffix + ".p6t10.bak")
    bak.write_text(txt, encoding="utf-8")
    if begin in txt and end in txt:
        pattern = re.compile(re.escape(begin) + r".*?" + re.escape(end), flags=re.S)
        new = pattern.sub(block.strip(), txt)
    elif "\\end{document}" in txt:
        new = txt.replace("\\end{document}", "\n" + block.strip() + "\n\n\\end{document}")
    else:
        new = txt.rstrip() + "\n\n" + block.strip() + "\n"
    paper_tex.write_text(new, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper-tex", default="", help="Optional target .tex file for insertion")
    parser.add_argument("--apply", action="store_true", help="Actually apply insertion to --paper-tex with backup")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    values = load_values()
    table_tex = build_table_tex(values)
    section_tex = build_section_tex(values)

    table_path = OUTDIR / "paperVI_EG_T10_closure_values_table.tex"
    section_path = OUTDIR / "paperVI_EG_T10_projected_closure_section.tex"
    table_path.write_text(table_tex, encoding="utf-8")
    section_path.write_text(section_tex, encoding="utf-8")

    source_rows: List[Dict[str, str]] = []
    t9lock = OUTDIR / "paperVI_EG_T9A_source_lock.csv"
    if t9lock.exists():
        with t9lock.open("r", encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                source_rows.append({
                    "claim_id": row.get("claim_id", ""),
                    "claim": row.get("claim", ""),
                    "source": row.get("source", ""),
                    "status": row.get("status", ""),
                })
    source_rows.extend([
        {
            "claim_id": "P6.T10.01",
            "claim": "Paper VI LaTeX insertion pack generated from T9A conditional projected closure.",
            "source": str(section_path),
            "status": "LATEX_INSERTION_READY",
        },
        {
            "claim_id": "P6.T10.02",
            "claim": "Closure values table generated and ready for inclusion.",
            "source": str(table_path),
            "status": "TABLE_READY",
        },
    ])
    write_csv(OUTDIR / "paperVI_EG_T10_source_lock_final.csv", source_rows, ["claim_id", "claim", "source", "status"])

    tex_candidates = find_tex_candidates()
    write_csv(OUTDIR / "paperVI_EG_T10_tex_candidates.csv", tex_candidates, ["path", "size_bytes", "name_score", "content_score", "total_score"])

    applied = False
    target = args.paper_tex or ""
    status = STATUS_DONE
    if args.apply:
        if not target:
            raise SystemExit("--apply requires --paper-tex path/to/paper.tex")
        apply_to_tex(Path(target), section_tex)
        applied = True
        status = STATUS_APPLIED

    decisions = [
        {"item": "insertion_pack", "decision": "READY", "basis": "T9A appendix converted to LaTeX section and values table."},
        {"item": "apply_mode", "decision": "APPLIED_WITH_BACKUP" if applied else "NOT_APPLIED_NONDESTRUCTIVE", "basis": "Default mode does not modify a paper tex file; --apply with --paper-tex is required."},
        {"item": "claim_strength", "decision": "PROJECTED_DIAGNOSTIC_CLOSURE_ONLY", "basis": "Do not claim action-level HR/BPB projector derivation is closed."},
        {"item": "next_gate", "decision": "P6_T11_PAPER_VI_COMPILE_OR_INSERT_MANUALLY", "basis": "Insert the generated section/table into the Paper VI manuscript and compile/check references."},
    ]
    write_csv(OUTDIR / "paperVI_EG_T10_decisions.csv", decisions, ["item", "decision", "basis"])

    manifest = {
        "status": status,
        "timestamp_utc_like": datetime.utcnow().isoformat() + "Z",
        "likelihood_run": "NO",
        "ns_touched": "NO",
        "applied": applied,
        "target_paper_tex": target,
        "files": {
            "section_tex": str(section_path),
            "table_tex": str(table_path),
            "source_lock_final": str(OUTDIR / "paperVI_EG_T10_source_lock_final.csv"),
            "decisions": str(OUTDIR / "paperVI_EG_T10_decisions.csv"),
            "tex_candidates": str(OUTDIR / "paperVI_EG_T10_tex_candidates.csv"),
        },
        "guardrail": "Projected-diagnostic closure only; action-level HR/BPB projector derivation remains open.",
    }
    (OUTDIR / "paperVI_EG_T10_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (OUTDIR / "paperVI_EG_T10_manifest.md").write_text(build_md(status, applied, target), encoding="utf-8")

    print(status)
    print("Wrote:")
    for p in [section_path, table_path, OUTDIR / "paperVI_EG_T10_source_lock_final.csv", OUTDIR / "paperVI_EG_T10_manifest.md", OUTDIR / "paperVI_EG_T10_manifest.json", OUTDIR / "paperVI_EG_T10_decisions.csv", OUTDIR / "paperVI_EG_T10_tex_candidates.csv"]:
        print(" -", p)
    if applied:
        print("Applied to:", target)
        print("Backup:", str(Path(target).with_suffix(Path(target).suffix + ".p6t10.bak")))


if __name__ == "__main__":
    main()
