#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import hashlib
from pathlib import Path
from datetime import datetime, timezone
from typing import Any

DEFAULT_CHECKPOINT_DIR = Path("outputs/ns/checkpoints_strict_global_gate")
DEFAULT_ROWS_CSV = Path("outputs/ns/run_multieos8_strict_global_checkpoint_gate_rows.csv")
DEFAULT_OUT_PREFIX = Path("outputs/ns/ns_checkpoint_table_audit")

MASS_TOKENS = {"m", "mass", "m_msun", "mmsun", "m_msun", "m_grav", "mgrav", "mgravitational", "m_sun", "msun", "m_msun"}
RADIUS_TOKENS = {"r", "radius", "r_km", "rkm", "rcirc", "rcirc_km"}
EOS_NAMES = ["BSk24", "APR", "DD2", "DDME2", "NL3wrL55", "SFHo", "SLy4", "TW"]

def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).lower()).strip("_")

def parse_float(x: Any) -> float | None:
    if x is None or isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        v = float(x)
        return v if math.isfinite(v) else None
    s = str(x).strip().replace("−", "-").replace(",", ".")
    if not s or s.lower() in {"nan", "none", "null", "true", "false"}:
        return None
    m = re.search(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        v = float(m.group(0))
    except Exception:
        return None
    return v if math.isfinite(v) else None

def sha256(p: Path) -> str:
    if not p.exists() or not p.is_file():
        return ""
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})

def read_csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def all_cols(rows: list[dict[str, Any]]) -> list[str]:
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    return cols

def find_mass_cols(cols: list[str]) -> list[str]:
    out = []
    for c in cols:
        nc = norm(c)
        if nc in MASS_TOKENS or nc.startswith("m_") or "mass" in nc or "msun" in nc:
            out.append(c)
    return out

def find_radius_cols(cols: list[str]) -> list[str]:
    out = []
    for c in cols:
        nc = norm(c)
        if nc in RADIUS_TOKENS or "radius" in nc or nc == "r_km" or nc.endswith("_km"):
            # avoid unrelated km-ish columns if possible
            if "lambda" not in nc and "mass" not in nc:
                out.append(c)
    return out

def find_lambda_cols(cols: list[str]) -> list[str]:
    return [c for c in cols if "lambda" in norm(c) or norm(c) in {"lam", "tidal_lambda"}]

def find_k2_cols(cols: list[str]) -> list[str]:
    return [c for c in cols if "k2" in norm(c)]

def eos_from_path_or_rows(path: Path, rows: list[dict[str, Any]]) -> str:
    s = str(path)
    for eos in EOS_NAMES:
        if eos.lower() in s.lower():
            return eos
    for r in rows[:20]:
        for k, v in r.items():
            if "eos" in norm(k) or "name" == norm(k):
                for eos in EOS_NAMES:
                    if eos.lower() in str(v).lower():
                        return eos
    return ""

def summarize_table(path: Path, table_path: str, rows: list[dict[str, Any]], source_kind: str) -> dict[str, Any]:
    cols = all_cols(rows)
    mass_cols = find_mass_cols(cols)
    radius_cols = find_radius_cols(cols)
    lambda_cols = find_lambda_cols(cols)
    k2_cols = find_k2_cols(cols)

    best_mass = ""
    best_unique = 0
    mass_min = ""
    mass_max = ""
    near_14 = ""
    bracket_14 = False
    bracket_width = ""
    for c in mass_cols:
        vals = [parse_float(r.get(c)) for r in rows]
        vals = [v for v in vals if v is not None and 0.0 < v < 5.0]
        if not vals:
            continue
        uniq = sorted(set(round(v, 12) for v in vals))
        if len(uniq) > best_unique:
            best_unique = len(uniq)
            best_mass = c
            mass_min = min(vals)
            mass_max = max(vals)
            below = [v for v in uniq if v <= 1.4]
            above = [v for v in uniq if v >= 1.4]
            if below and above:
                lo = max(below); hi = min(above)
                bracket_14 = True
                bracket_width = hi - lo
                near_14 = min(uniq, key=lambda x: abs(x - 1.4))
            elif uniq:
                near_14 = min(uniq, key=lambda x: abs(x - 1.4))

    score = 0
    score += min(len(rows), 10000) / 100
    score += 50 if best_unique >= 100 else 0
    score += 20 if best_unique >= 30 else 0
    score += 30 if lambda_cols else 0
    score += 10 if radius_cols else 0
    score += 10 if bracket_14 else 0
    if isinstance(bracket_width, float) and bracket_width <= 0.03:
        score += 50

    verdict = "LOW_VALUE"
    if best_unique >= 100 and bracket_14 and lambda_cols:
        verdict = "DENSE_LAMBDA_SEQUENCE_CANDIDATE"
    elif best_unique >= 30 and bracket_14 and lambda_cols:
        verdict = "SPARSE_LAMBDA_SEQUENCE_CANDIDATE"
    elif best_unique >= 30 and bracket_14:
        verdict = "MASS_RADIUS_SEQUENCE_NO_LAMBDA"
    elif lambda_cols:
        verdict = "LAMBDA_COLUMNS_NO_GOOD_MASS_GRID"

    return {
        "source": str(path),
        "source_kind": source_kind,
        "table_path": table_path,
        "eos_hint": eos_from_path_or_rows(path, rows),
        "n_rows": len(rows),
        "n_columns": len(cols),
        "score": score,
        "verdict": verdict,
        "mass_col_best": best_mass,
        "n_unique_mass_best": best_unique,
        "mass_min": mass_min,
        "mass_max": mass_max,
        "nearest_mass_to_1p4": near_14,
        "brackets_1p4": bracket_14,
        "bracket_width": bracket_width,
        "radius_cols": ";".join(radius_cols[:20]),
        "lambda_cols": ";".join(lambda_cols[:20]),
        "k2_cols": ";".join(k2_cols[:20]),
        "columns_preview": ";".join(cols[:60]),
        "sha256": sha256(path),
    }

def recurse_json_tables(obj: Any, path: str = "$") -> list[tuple[str, list[dict[str, Any]]]]:
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            child = f"{path}.{k}"
            if isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
                out.append((child, v))
            out.extend(recurse_json_tables(v, child))
    elif isinstance(obj, list):
        for i, v in enumerate(obj[:20]):  # avoid huge path explosion; tables are caught at list level
            out.extend(recurse_json_tables(v, f"{path}[{i}]"))
    return out

def main():
    ap = argparse.ArgumentParser(description="Audit NS checkpoint JSON/CSV table schemas before 1.4 paper-table refinement.")
    ap.add_argument("--checkpoint-dir", default=str(DEFAULT_CHECKPOINT_DIR))
    ap.add_argument("--rows-csv", default=str(DEFAULT_ROWS_CSV))
    ap.add_argument("--out-prefix", default=str(DEFAULT_OUT_PREFIX))
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    checkpoint_dir = Path(args.checkpoint_dir)
    rows_csv = Path(args.rows_csv)

    tables = []
    source_lock = []

    csv_rows = read_csv_rows(rows_csv)
    if csv_rows:
        tables.append(summarize_table(rows_csv, "$CSV", csv_rows, "rows_csv"))
    source_lock.append({
        "path": str(rows_csv),
        "kind": "rows_csv",
        "exists": rows_csv.exists(),
        "sha256": sha256(rows_csv),
        "n_rows": len(csv_rows),
    })

    if checkpoint_dir.exists():
        for jp in sorted(checkpoint_dir.glob("*.json")):
            source_lock.append({
                "path": str(jp),
                "kind": "checkpoint_json",
                "exists": True,
                "sha256": sha256(jp),
                "n_rows": "",
            })
            try:
                obj = json.loads(jp.read_text(encoding="utf-8"))
            except Exception as exc:
                tables.append({
                    "source": str(jp),
                    "source_kind": "checkpoint_json_error",
                    "table_path": "",
                    "verdict": "JSON_READ_ERROR",
                    "error": repr(exc),
                    "sha256": sha256(jp),
                })
                continue
            for table_path, rows in recurse_json_tables(obj):
                tables.append(summarize_table(jp, table_path, rows, "checkpoint_json"))
    else:
        source_lock.append({
            "path": str(checkpoint_dir),
            "kind": "checkpoint_dir",
            "exists": False,
            "sha256": "",
            "n_rows": "",
        })

    tables.sort(key=lambda r: float(r.get("score") or 0), reverse=True)

    strong = [r for r in tables if r.get("verdict") == "DENSE_LAMBDA_SEQUENCE_CANDIDATE"]
    usable_lambda = [r for r in tables if "LAMBDA" in str(r.get("verdict"))]
    mass_radius = [r for r in tables if r.get("verdict") == "MASS_RADIUS_SEQUENCE_NO_LAMBDA"]

    decisions = [
        {
            "item": "dense_lambda_sequence",
            "decision": "FOUND" if strong else "NOT_FOUND",
            "basis": f"{len(strong)} dense bracketing tables with lambda columns and >=100 unique masses.",
        },
        {
            "item": "any_lambda_sequence",
            "decision": "FOUND" if usable_lambda else "NOT_FOUND",
            "basis": f"{len(usable_lambda)} candidate tables with lambda columns.",
        },
        {
            "item": "mass_radius_only",
            "decision": "FOUND" if mass_radius else "NOT_FOUND",
            "basis": f"{len(mass_radius)} mass/radius tables bracket 1.4 but have no detected Lambda columns.",
        },
        {
            "item": "next_gate",
            "decision": "USE_DETECTED_DENSE_TABLE_OR_TARGETED_SOLVE",
            "basis": "If no dense lambda sequence is found, final Paper IV 1.4 Lambda values require a targeted 1.4 solve rather than interpolation from checkpoint summaries.",
        },
    ]

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_csv(out_prefix.with_name(out_prefix.name + "_tables.csv"), tables)
    write_csv(out_prefix.with_name(out_prefix.name + "_source_lock.csv"), source_lock)
    write_csv(out_prefix.with_name(out_prefix.name + "_decisions.csv"), decisions)

    payload = {
        "status": "NS_CHECKPOINT_TABLE_AUDIT_DONE",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "n_tables": len(tables),
        "n_dense_lambda_sequence_candidates": len(strong),
        "n_any_lambda_candidates": len(usable_lambda),
        "top_tables": tables[:20],
        "decisions": decisions,
    }
    out_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = []
    md.append("# NS checkpoint table audit\n")
    md.append("## Status\n")
    md.append("**NS_CHECKPOINT_TABLE_AUDIT_DONE**\n")
    md.append("Purpose: locate the actual table(s) inside strict-global checkpoint JSON/CSV before final 1.4-Msun Paper IV interpolation.\n")
    md.append("## Counts\n")
    md.append(f"- tables found: {len(tables)}")
    md.append(f"- dense lambda sequence candidates: {len(strong)}")
    md.append(f"- any lambda candidates: {len(usable_lambda)}")
    md.append(f"- mass/radius-only bracketing candidates: {len(mass_radius)}\n")
    md.append("## Top table candidates\n")
    headers = ["score", "verdict", "eos_hint", "source", "table_path", "n_rows", "n_unique_mass_best", "mass_col_best", "bracket_width", "radius_cols", "lambda_cols"]
    md.append("| " + " | ".join(headers) + " |")
    md.append("|" + "|".join(["---"] * len(headers)) + "|")
    for r in tables[:25]:
        vals = []
        for h in headers:
            v = r.get(h, "")
            if isinstance(v, float):
                v = f"{v:.6g}"
            vals.append(str(v).replace("|", "/")[:140])
        md.append("| " + " | ".join(vals) + " |")
    md.append("\n## Decisions\n")
    md.append("| item | decision | basis |")
    md.append("|---|---|---|")
    for d in decisions:
        md.append(f"| {d['item']} | {d['decision']} | {d['basis']} |")
    md.append("\n## Interpretation rule\n")
    md.append("- If a `DENSE_LAMBDA_SEQUENCE_CANDIDATE` exists, use that table path for final local interpolation.")
    md.append("- If only `MASS_RADIUS_SEQUENCE_NO_LAMBDA` exists, the checkpoint has background rows but not final tidal Lambda rows.")
    md.append("- If only sparse lambda candidates exist, final Lambda at 1.4 requires a targeted solve/dense run around 1.4.")
    out_prefix.with_suffix(".md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("Saved:", out_prefix.with_suffix(".md"))
    print("Saved:", out_prefix.with_suffix(".json"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_tables.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_decisions.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_source_lock.csv"))

if __name__ == "__main__":
    main()
