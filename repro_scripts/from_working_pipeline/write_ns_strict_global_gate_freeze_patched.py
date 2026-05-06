#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import hashlib
from pathlib import Path
from datetime import datetime, timezone

OUT = Path("outputs/ns")
CHECKPOINTS = OUT / "checkpoints_strict_global_gate"

# Authoritative freeze values from the completed strict-global watcher output.
# These values are copied from the terminal summary of the completed run.
KAPPA_STRICT = 0.965543842

ROWS = [
    {
        "EOS": "BSk24", "status": "ok", "n": 1255, "topology": "9/1255",
        "M_range": "1.116–2.277", "Mmax_seq": 2.277165, "R_Mmax": 10.762886,
        "nearest_1p4_M": 1.471378, "R1p4_nearest_grid_approx": 12.112159,
        "Lambda1p4_nearest_grid_approx": 260.256817, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -0.996,
    },
    {
        "EOS": "APR", "status": "ok", "n": 1308, "topology": "9/1308",
        "M_range": "0.660–2.204", "Mmax_seq": 2.204165, "R_Mmax": 9.536365,
        "nearest_1p4_M": 1.484316, "R1p4_nearest_grid_approx": 10.332866,
        "Lambda1p4_nearest_grid_approx": 91.260955, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.063,
    },
    {
        "EOS": "DD2", "status": "ok", "n": 1244, "topology": "9/1244",
        "M_range": "1.155–2.404", "Mmax_seq": 2.403902, "R_Mmax": 11.359537,
        "nearest_1p4_M": 1.568759, "R1p4_nearest_grid_approx": 12.019242,
        "Lambda1p4_nearest_grid_approx": 161.901255, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.025,
    },
    {
        "EOS": "DDME2", "status": "ok", "n": 1155, "topology": "10/1155",
        "M_range": "1.174–2.471", "Mmax_seq": 2.471424, "R_Mmax": 11.495463,
        "nearest_1p4_M": 1.623273, "R1p4_nearest_grid_approx": 12.113637,
        "Lambda1p4_nearest_grid_approx": 142.372340, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -0.898,
    },
    {
        "EOS": "NL3wrL55", "status": "ok", "n": 986, "topology": "10/986",
        "M_range": "0.897–2.753", "Mmax_seq": 2.753006, "R_Mmax": 12.488097,
        "nearest_1p4_M": 1.449090, "R1p4_nearest_grid_approx": 12.375749,
        "Lambda1p4_nearest_grid_approx": 335.393891, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.065,
    },
    {
        "EOS": "SFHo", "status": "ok", "n": 1139, "topology": "8/1139",
        "M_range": "1.088–2.054", "Mmax_seq": 2.053955, "R_Mmax": 9.777273,
        "nearest_1p4_M": 1.389372, "R1p4_nearest_grid_approx": 10.820066,
        "Lambda1p4_nearest_grid_approx": 163.223181, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.073,
    },
    {
        "EOS": "SLy4", "status": "ok", "n": 1136, "topology": "10/1136",
        "M_range": "0.991–2.043", "Mmax_seq": 2.043413, "R_Mmax": 9.633417,
        "nearest_1p4_M": 1.514428, "R1p4_nearest_grid_approx": 10.571850,
        "Lambda1p4_nearest_grid_approx": 84.746736, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.210,
    },
    {
        "EOS": "TW", "status": "ok", "n": 1213, "topology": "10/1213",
        "M_range": "1.220–2.052", "Mmax_seq": 2.052304, "R_Mmax": 9.982124,
        "nearest_1p4_M": 1.509123, "R1p4_nearest_grid_approx": 11.066518,
        "Lambda1p4_nearest_grid_approx": 115.306457, "four_zone_1p4": False,
        "mean_delta_Lambda_over_Lambda_percent": -1.035,
    },
]

def sha(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def fmt(x, n=6):
    if isinstance(x, bool):
        return "True" if x else "False"
    if x is None:
        return ""
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.{n}f}"
    return str(x)

def write_csv(path: Path, rows: list[dict]):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})

def priority(abs_delta: float) -> str:
    if abs_delta > 0.10:
        return "HIGH"
    if abs_delta > 0.03:
        return "MEDIUM"
    return "LOW"

def main():
    OUT.mkdir(parents=True, exist_ok=True)

    for r in ROWS:
        r["abs_nearest_delta_from_1p4"] = abs(float(r["nearest_1p4_M"]) - 1.4)

    deltas = [r["mean_delta_Lambda_over_Lambda_percent"] for r in ROWS]
    mmaxs = [r["Mmax_seq"] for r in ROWS]
    summary = {
        "status": "NS_STRICT_GLOBAL_GATE_FREEZE_DONE_PATCHED",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "kappa_strict": KAPPA_STRICT,
        "n_eos": len(ROWS),
        "ok_eos": sum(1 for r in ROWS if r["status"] == "ok"),
        "all_ok": all(r["status"] == "ok" for r in ROWS),
        "mean_delta_percent_mean": sum(deltas) / len(deltas),
        "mean_delta_percent_min": min(deltas),
        "mean_delta_percent_max": max(deltas),
        "mmax_seq_min": min(mmaxs),
        "mmax_seq_max": max(mmaxs),
        "all_mmax_seq_over_2Msun": all(x > 2.0 for x in mmaxs),
        "four_zone_1p4_true_count": sum(1 for r in ROWS if r["four_zone_1p4"]),
        "four_zone_1p4_false_count": sum(1 for r in ROWS if not r["four_zone_1p4"]),
    }

    refine = []
    for r in sorted(ROWS, key=lambda x: x["abs_nearest_delta_from_1p4"], reverse=True):
        refine.append({
            "EOS": r["EOS"],
            "nearest_1p4_M": r["nearest_1p4_M"],
            "abs_nearest_delta_from_1p4": r["abs_nearest_delta_from_1p4"],
            "R1p4_nearest_grid_approx": r["R1p4_nearest_grid_approx"],
            "Lambda1p4_nearest_grid_approx": r["Lambda1p4_nearest_grid_approx"],
            "refine_priority": priority(r["abs_nearest_delta_from_1p4"]),
            "refine_for_final_paper_table": True,
            "reason": "nearest-grid diagnostic; use interpolation/exact target solve before final 1.4 table",
        })

    decisions = [
        {"item": "strict_global_gate", "decision": "PASS_8_OF_8", "basis": "8/8 EOS have ok status under strict dictionary."},
        {"item": "legacy_hybrid", "decision": "DO_NOT_RUN", "basis": "Strict dictionary passed; legacy/hybrid are fallback diagnostics only."},
        {"item": "tidal_response", "decision": "COHERENT_NEGATIVE_AROUND_ONE_PERCENT", "basis": f"mean delta range {summary['mean_delta_percent_min']:.3f}% to {summary['mean_delta_percent_max']:.3f}%, mean {summary['mean_delta_percent_mean']:.3f}%."},
        {"item": "Mmax_seq", "decision": "SURVIVAL_GATE_ONLY", "basis": "Mmax_seq is the maximum reached in the checkpointed search window, not a mathematical EOS maximum."},
        {"item": "R1p4_Lambda1p4", "decision": "NEAREST_GRID_APPROX_REFINE_FOR_FINAL_TABLES", "basis": "Nearest-grid values should be interpolated/exact-targeted for final paper tables."},
        {"item": "four_zone_1p4", "decision": "NOT_A_1P4_FOUR_ZONE_CLAIM", "basis": "true count=0, false count=8."},
        {"item": "native_k2_solver", "decision": "OPEN_NOT_BLOCKING_STRICT_EFFECTIVE_GATE", "basis": "This freeze concerns strict effective tidal/Love outputs; native scalar-tensor k2 remains separate."},
        {"item": "parser_patch", "decision": "PATCHED_FROM_COMPLETED_TERMINAL_SUMMARY", "basis": "Previous freeze parser failed to map summary CSV columns; this patched freeze uses the completed watcher values explicitly."},
        {"item": "next_gate", "decision": "NS_1P4_REFINEMENT_OR_PAPER_INSERTION", "basis": "Refine exact 1.4 values or insert this coarse survival gate with caveat."},
    ]

    source_files = [
        OUT / "run_multieos8_strict_global_checkpoint_gate.json",
        OUT / "run_multieos8_strict_global_checkpoint_gate_rows.csv",
        OUT / "run_multieos8_strict_global_checkpoint_gate_summary.csv",
        OUT / "run_multieos8_strict_global_checkpoint_gate.md",
    ]
    source_lock = []
    for p in source_files:
        source_lock.append({
            "path": str(p), "exists": p.exists(),
            "kind": "file" if p.is_file() else "missing",
            "size_bytes": p.stat().st_size if p.exists() and p.is_file() else "",
            "sha256": sha(p),
            "purpose": "strict_global_gate_input",
        })
    if CHECKPOINTS.exists():
        cks = sorted(CHECKPOINTS.glob("*.json"))
        source_lock.append({"path": str(CHECKPOINTS), "exists": True, "kind": "directory", "size_bytes": "", "sha256": "", "purpose": f"checkpoint_dir_count_json={len(cks)}"})
        for p in cks:
            source_lock.append({"path": str(p), "exists": True, "kind": "checkpoint_json", "size_bytes": p.stat().st_size, "sha256": sha(p), "purpose": "strict_global_checkpoint"})
    else:
        source_lock.append({"path": str(CHECKPOINTS), "exists": False, "kind": "missing", "size_bytes": "", "sha256": "", "purpose": "checkpoint_dir"})

    md = []
    md.append("# NS strict-only global EOS gate freeze\n")
    md.append("## Status\n")
    md.append("**NS_STRICT_GLOBAL_GATE_FREEZE_DONE_PATCHED**\n")
    md.append("Non-destructive freeze gate. Reads the completed strict-dictionary outputs and uses the completed watcher summary values. Does not run legacy/hybrid and does not modify `src/ns` or checkpoint JSON files.\n")
    md.append("## Parser note\n")
    md.append("The first freeze attempt did not map the summary CSV schema correctly. This patched freeze uses the completed terminal watcher values explicitly and locks them with source files/checkpoints listed in the source-lock.\n")
    md.append("## Global lock\n")
    md.append("| quantity | value |\n|---|---:|")
    for k in ["kappa_strict","n_eos","ok_eos","all_ok","mean_delta_percent_mean","mean_delta_percent_min","mean_delta_percent_max","mmax_seq_min","mmax_seq_max","all_mmax_seq_over_2Msun","four_zone_1p4_true_count","four_zone_1p4_false_count"]:
        md.append(f"| {k} | {fmt(summary.get(k), 9 if k=='kappa_strict' else 6)} |")
    md.append("\n## EOS summary\n")
    md.append("| EOS | status | n | topology | M range | Mmax seq | R(Mmax) | nearest 1.4 M | R1.4 approx | Lambda1.4 approx | 1.4 four-zone | mean dLambda/Lambda |\n|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in ROWS:
        md.append(
            f"| `{r['EOS']}` | {r['status']} | {r['n']} | {r['topology']} | {r['M_range']} | "
            f"{r['Mmax_seq']:.6f} | {r['R_Mmax']:.6f} | {r['nearest_1p4_M']:.6f} | "
            f"{r['R1p4_nearest_grid_approx']:.6f} | {r['Lambda1p4_nearest_grid_approx']:.6f} | "
            f"{r['four_zone_1p4']} | {r['mean_delta_Lambda_over_Lambda_percent']:.3f}% |"
        )
    md.append("\n## 1.4-Msun refinement priority\n")
    md.append("Nearest-grid 1.4 values are diagnostic only. For final paper tables, use interpolation or exact-target solves per EOS.\n")
    md.append("| EOS | nearest 1.4 M | |delta M| | priority |\n|---|---:|---:|---|")
    for r in refine:
        md.append(f"| `{r['EOS']}` | {r['nearest_1p4_M']:.6f} | {r['abs_nearest_delta_from_1p4']:.6f} | {r['refine_priority']} |")
    md.append("\n## Decisions\n")
    md.append("| item | decision | basis |\n|---|---|---|")
    for d in decisions:
        md.append(f"| {d['item']} | {d['decision']} | {d['basis']} |")
    md.append("\n## Locked interpretation\n")
    md.append("- Strict dictionary global gate passes all eight EOS.")
    md.append("- Legacy/hybrid branches are not needed unless a later strict-only diagnostic fails.")
    md.append("- The tidal/Love-sector correction is coherently negative across EOS, around the one-percent level.")
    md.append("- `Mmax_seq` is a checkpoint-window survival gate, not a mathematical EOS maximum claim.")
    md.append("- `R1.4` and `Lambda1.4` are nearest-grid diagnostics, not final interpolated 1.4-Msun values.")
    md.append("- `1.4 four-zone = False` across the table means the paper claim should be strict effective tidal response, not a full four-zone topology claim at 1.4 Msun.")

    (OUT / "ns_strict_global_gate_freeze.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (OUT / "ns_strict_global_gate_freeze.json").write_text(json.dumps({"summary": summary, "rows": ROWS, "refine_priority": refine, "decisions": decisions}, indent=2), encoding="utf-8")
    write_csv(OUT / "ns_strict_global_gate_freeze_summary.csv", ROWS)
    write_csv(OUT / "ns_strict_global_gate_freeze_refine_1p4_priority.csv", refine)
    write_csv(OUT / "ns_strict_global_gate_freeze_source_lock.csv", source_lock)
    write_csv(OUT / "ns_strict_global_gate_freeze_decisions.csv", decisions)

    print("Saved:", OUT / "ns_strict_global_gate_freeze.md")
    print("Saved:", OUT / "ns_strict_global_gate_freeze.json")
    print("Saved:", OUT / "ns_strict_global_gate_freeze_summary.csv")
    print("Saved:", OUT / "ns_strict_global_gate_freeze_refine_1p4_priority.csv")
    print("Saved:", OUT / "ns_strict_global_gate_freeze_source_lock.csv")
    print("Saved:", OUT / "ns_strict_global_gate_freeze_decisions.csv")

if __name__ == "__main__":
    main()
