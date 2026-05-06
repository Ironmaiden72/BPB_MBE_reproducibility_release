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

DEFAULT_ROWS_CSV = Path("outputs/ns/run_multieos8_strict_global_checkpoint_gate_rows.csv")
DEFAULT_SUMMARY_CSV = Path("outputs/ns/ns_strict_global_gate_freeze_summary.csv")
DEFAULT_CHECKPOINT_DIR = Path("outputs/ns/checkpoints_strict_global_gate")
DEFAULT_OUT_PREFIX = Path("outputs/ns/ns_1p4_refinement_from_strict_gate")
TARGET_DEFAULT = 1.4

EOS_NAMES = ["BSk24", "APR", "DD2", "DDME2", "NL3wrL55", "SFHo", "SLy4", "TW"]

MASS_ALIASES = [
    "M", "M_msun", "M_Msun", "mass", "mass_msun", "M_solar", "Mgrav",
    "M_grav", "M_G", "Mstar", "mstar", "m_msun", "gravitational_mass",
]
RADIUS_ALIASES = [
    "R", "R_km", "radius", "radius_km", "Rcirc", "Rcirc_km", "R_kilometer",
]
EOS_ALIASES = ["EOS", "eos", "eos_name", "name", "template", "eos_label"]

LAMBDA_GR_HINTS = [
    "lambda_gr", "lambda_tov", "lambda_gr_tidal", "lambda_lcdm",
    "lambda_background", "lambda_base", "lambda_ref",
]
LAMBDA_BPB_HINTS = [
    "lambda_bpb", "lambda_st", "lambda_strict", "lambda_eff", "lambda_bpb_eff",
    "lambda_tidal_bpb", "lambda_tidal_st", "lambda_modified",
]
K2_GR_HINTS = ["k2_gr", "k2_tov", "k2_background", "k2_base", "k2_ref"]
K2_BPB_HINTS = ["k2_bpb", "k2_st", "k2_strict", "k2_eff", "k2_bpb_eff"]

def sha256(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower()).strip("_")

def parse_float(x: Any) -> float | None:
    if x is None:
        return None
    if isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        if math.isfinite(float(x)):
            return float(x)
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "null", "false", "true"}:
        return None
    s = s.replace("−", "-").replace(",", ".")
    # Strip units/percent but keep exponent notation.
    m = re.search(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        v = float(m.group(0))
    except ValueError:
        return None
    return v if math.isfinite(v) else None

def parse_bool(x: Any) -> bool | None:
    if isinstance(x, bool):
        return x
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"true", "t", "yes", "y", "1"}:
        return True
    if s in {"false", "f", "no", "n", "0"}:
        return False
    return None

def read_csv_dicts(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))

def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})

def find_col(cols: list[str], aliases: list[str]) -> str | None:
    nmap = {norm_name(c): c for c in cols}
    for a in aliases:
        na = norm_name(a)
        if na in nmap:
            return nmap[na]
    # fuzzy fallback: exact token containment
    for a in aliases:
        na = norm_name(a)
        for nc, c in nmap.items():
            if na and (nc == na or nc.endswith("_" + na) or nc.startswith(na + "_")):
                return c
    return None

def detect_eos_from_path(path: Path) -> str | None:
    s = str(path)
    for eos in EOS_NAMES:
        if eos.lower() in s.lower():
            return eos
    return None

def extract_lists_from_json(obj: Any, eos_hint: str | None = None) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []

    def rec(x: Any, ctx: dict[str, Any] | None = None):
        ctx = dict(ctx or {})
        if isinstance(x, dict):
            # Keep EOS-like context values.
            for k in EOS_ALIASES:
                if k in x and x[k]:
                    ctx["EOS"] = x[k]
            for key, val in x.items():
                if isinstance(val, list):
                    if val and all(isinstance(v, dict) for v in val):
                        for row in val:
                            rr = dict(ctx)
                            rr.update(row)
                            out.append(rr)
                    else:
                        rec(val, ctx)
                elif isinstance(val, dict):
                    rec(val, ctx)
        elif isinstance(x, list):
            if x and all(isinstance(v, dict) for v in x):
                for row in x:
                    rr = dict(ctx)
                    rr.update(row)
                    out.append(rr)
            else:
                for item in x:
                    rec(item, ctx)

    rec(obj, {"EOS": eos_hint} if eos_hint else {})
    return out

def load_rows(rows_csv: Path, checkpoint_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sources = []
    rows = read_csv_dicts(rows_csv)
    if rows:
        sources.append({
            "path": str(rows_csv),
            "exists": True,
            "used": True,
            "kind": "rows_csv",
            "n_rows_loaded": len(rows),
            "sha256": sha256(rows_csv),
        })
        return rows, sources

    sources.append({
        "path": str(rows_csv),
        "exists": rows_csv.exists(),
        "used": False,
        "kind": "rows_csv",
        "n_rows_loaded": 0,
        "sha256": sha256(rows_csv),
    })

    json_rows: list[dict[str, Any]] = []
    if checkpoint_dir.exists():
        for jp in sorted(checkpoint_dir.glob("*.json")):
            try:
                obj = json.loads(jp.read_text(encoding="utf-8"))
                extracted = extract_lists_from_json(obj, detect_eos_from_path(jp))
                json_rows.extend(extracted)
                sources.append({
                    "path": str(jp),
                    "exists": True,
                    "used": bool(extracted),
                    "kind": "checkpoint_json",
                    "n_rows_loaded": len(extracted),
                    "sha256": sha256(jp),
                })
            except Exception as exc:
                sources.append({
                    "path": str(jp),
                    "exists": True,
                    "used": False,
                    "kind": "checkpoint_json_error",
                    "n_rows_loaded": 0,
                    "error": repr(exc),
                    "sha256": sha256(jp),
                })
    else:
        sources.append({
            "path": str(checkpoint_dir),
            "exists": False,
            "used": False,
            "kind": "checkpoint_dir_missing",
            "n_rows_loaded": 0,
            "sha256": "",
        })

    return json_rows, sources

def all_columns(rows: list[dict[str, Any]]) -> list[str]:
    cols: list[str] = []
    for r in rows:
        for k in r.keys():
            if k not in cols:
                cols.append(k)
    return cols

def pick_numeric_columns(rows: list[dict[str, Any]], min_count: int = 2) -> list[str]:
    cols = all_columns(rows)
    nums = []
    for c in cols:
        count = sum(1 for r in rows if parse_float(r.get(c)) is not None)
        if count >= min_count:
            nums.append(c)
    return nums

def identify_lambda_columns(cols: list[str]) -> tuple[str | None, str | None]:
    norm = {c: norm_name(c) for c in cols}
    lambda_cols = [c for c in cols if "lambda" in norm[c] or norm[c] in {"lam", "lamb"}]
    if not lambda_cols:
        return None, None

    def score(c: str, hints: list[str]) -> int:
        nc = norm[c]
        s = 0
        for h in hints:
            if norm_name(h) in nc:
                s += 10
        if "gr" in nc or "lcdm" in nc or "base" in nc or "ref" in nc:
            s += 2 if hints is LAMBDA_GR_HINTS else -1
        if "bpb" in nc or "st" in nc or "strict" in nc or "eff" in nc:
            s += 2 if hints is LAMBDA_BPB_HINTS else -1
        return s

    gr = max(lambda_cols, key=lambda c: score(c, LAMBDA_GR_HINTS))
    bpb = max(lambda_cols, key=lambda c: score(c, LAMBDA_BPB_HINTS))
    if score(gr, LAMBDA_GR_HINTS) <= 0:
        gr = None
    if score(bpb, LAMBDA_BPB_HINTS) <= 0:
        bpb = None
    if gr == bpb:
        # If only one Lambda column exists, keep it as BPB/observed-like, no delta.
        bpb = gr
        gr = None
    return gr, bpb

def identify_k2_columns(cols: list[str]) -> tuple[str | None, str | None]:
    norm = {c: norm_name(c) for c in cols}
    k2_cols = [c for c in cols if "k2" in norm[c]]
    if not k2_cols:
        return None, None

    def score(c: str, hints: list[str]) -> int:
        nc = norm[c]
        s = 0
        for h in hints:
            if norm_name(h) in nc:
                s += 10
        if "gr" in nc or "base" in nc or "ref" in nc:
            s += 2 if hints is K2_GR_HINTS else -1
        if "bpb" in nc or "st" in nc or "strict" in nc or "eff" in nc:
            s += 2 if hints is K2_BPB_HINTS else -1
        return s

    gr = max(k2_cols, key=lambda c: score(c, K2_GR_HINTS))
    bpb = max(k2_cols, key=lambda c: score(c, K2_BPB_HINTS))
    if score(gr, K2_GR_HINTS) <= 0:
        gr = None
    if score(bpb, K2_BPB_HINTS) <= 0:
        bpb = None
    if gr == bpb:
        bpb = gr
        gr = None
    return gr, bpb

def unique_by_mass(points: list[dict[str, Any]], mass_col: str, numeric_cols: list[str]) -> list[dict[str, Any]]:
    bins: dict[float, list[dict[str, Any]]] = {}
    for r in points:
        m = parse_float(r.get(mass_col))
        if m is None:
            continue
        bins.setdefault(m, []).append(r)

    out = []
    for m, group in bins.items():
        rr = dict(group[0])
        for c in numeric_cols:
            vals = [parse_float(g.get(c)) for g in group]
            vals = [v for v in vals if v is not None]
            if vals:
                rr[c] = sum(vals) / len(vals)
        rr[mass_col] = m
        out.append(rr)
    out.sort(key=lambda r: float(r[mass_col]))
    return out

def interp_value(x0: float, y0: float, x1: float, y1: float, x: float) -> float:
    if abs(x1 - x0) < 1e-14:
        return y0
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

def bracket(points: list[dict[str, Any]], mass_col: str, target: float) -> tuple[str, dict[str, Any], dict[str, Any]]:
    pts = points
    if not pts:
        raise ValueError("no points")
    # exact
    for p in pts:
        m = parse_float(p.get(mass_col))
        if m is not None and abs(m - target) < 1e-8:
            return "exact", p, p

    below = [p for p in pts if parse_float(p.get(mass_col)) is not None and parse_float(p.get(mass_col)) <= target]
    above = [p for p in pts if parse_float(p.get(mass_col)) is not None and parse_float(p.get(mass_col)) >= target]
    if below and above:
        lo = max(below, key=lambda p: parse_float(p.get(mass_col)))
        hi = min(above, key=lambda p: parse_float(p.get(mass_col)))
        if lo is hi:
            return "exact", lo, hi
        return "bracketed_linear", lo, hi

    nearest = min(pts, key=lambda p: abs(parse_float(p.get(mass_col)) - target))
    return "nearest_only_no_bracket", nearest, nearest

def interpolate_eos(eos: str, points: list[dict[str, Any]], mass_col: str, numeric_cols: list[str], target: float) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pts = unique_by_mass(points, mass_col, numeric_cols)
    mode, lo, hi = bracket(pts, mass_col, target)
    mlo = parse_float(lo.get(mass_col))
    mhi = parse_float(hi.get(mass_col))

    out = {
        "EOS": eos,
        "target_M": target,
        "n_points": len(pts),
        "mass_min": min(parse_float(p.get(mass_col)) for p in pts),
        "mass_max": max(parse_float(p.get(mass_col)) for p in pts),
        "interp_mode": mode,
        "M_lo": mlo,
        "M_hi": mhi,
        "bracket_width": None if mlo is None or mhi is None else abs(mhi - mlo),
    }
    per_col = []
    for c in numeric_cols:
        vlo = parse_float(lo.get(c))
        vhi = parse_float(hi.get(c))
        if vlo is None:
            continue
        if mode in {"exact", "nearest_only_no_bracket"} or vhi is None:
            val = vlo
        else:
            val = interp_value(mlo, vlo, mhi, vhi, target)
        key = f"{norm_name(c)}_1p4"
        out[key] = val
        per_col.append({
            "EOS": eos,
            "column": c,
            "output_key": key,
            "interp_mode": mode,
            "M_lo": mlo,
            "value_lo": vlo,
            "M_hi": mhi,
            "value_hi": vhi,
            "value_1p4": val,
        })

    # Booleans: nearest if no clean interpolation.
    for c in all_columns(points):
        vals = [parse_bool(p.get(c)) for p in points]
        if any(v is not None for v in vals):
            blo = parse_bool(lo.get(c))
            bhi = parse_bool(hi.get(c))
            if mode == "exact" or blo == bhi:
                out[f"{norm_name(c)}_1p4_bool"] = blo

    return out, per_col

def fmt(x: Any, n: int = 6) -> str:
    if x is None or x == "":
        return ""
    if isinstance(x, bool):
        return "True" if x else "False"
    if isinstance(x, (int,)):
        return str(x)
    if isinstance(x, float):
        return f"{x:.{n}f}"
    return str(x)

def main():
    ap = argparse.ArgumentParser(description="Interpolate strict-global NS checkpoint rows at M=1.4 Msun without rerunning solvers.")
    ap.add_argument("--rows-csv", default=str(DEFAULT_ROWS_CSV))
    ap.add_argument("--summary-csv", default=str(DEFAULT_SUMMARY_CSV))
    ap.add_argument("--checkpoint-dir", default=str(DEFAULT_CHECKPOINT_DIR))
    ap.add_argument("--target-mass", type=float, default=TARGET_DEFAULT)
    ap.add_argument("--out-prefix", default=str(DEFAULT_OUT_PREFIX))
    ap.add_argument("--eos", nargs="*", default=EOS_NAMES)
    args = ap.parse_args()

    rows_csv = Path(args.rows_csv)
    checkpoint_dir = Path(args.checkpoint_dir)
    out_prefix = Path(args.out_prefix)
    target = args.target_mass

    rows, source_rows = load_rows(rows_csv, checkpoint_dir)
    if not rows:
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        msg = {
            "status": "NS_1P4_REFINEMENT_FAILED_NO_ROWS",
            "rows_csv": str(rows_csv),
            "checkpoint_dir": str(checkpoint_dir),
            "sources": source_rows,
        }
        (out_prefix.with_suffix(".json")).write_text(json.dumps(msg, indent=2), encoding="utf-8")
        raise SystemExit("No usable rows found. Check rows CSV or checkpoint JSON schema.")

    cols = all_columns(rows)
    eos_col = find_col(cols, EOS_ALIASES)
    mass_col = find_col(cols, MASS_ALIASES)
    radius_col = find_col(cols, RADIUS_ALIASES)

    if mass_col is None:
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        msg = {
            "status": "NS_1P4_REFINEMENT_FAILED_NO_MASS_COLUMN",
            "columns": cols,
            "sources": source_rows,
        }
        (out_prefix.with_suffix(".json")).write_text(json.dumps(msg, indent=2), encoding="utf-8")
        raise SystemExit("Could not detect a mass column. Inspect columns in JSON output.")

    # Add inferred EOS if absent.
    if eos_col is None:
        eos_col = "EOS"
        for r in rows:
            if not r.get("EOS"):
                r["EOS"] = r.get("eos") or r.get("template") or "UNKNOWN"

    numeric_cols = pick_numeric_columns(rows)
    # Ensure mass col is numeric and not double-output as a physical column.
    if mass_col not in numeric_cols:
        numeric_cols.append(mass_col)

    # Group by EOS, but keep only requested EOS when possible.
    groups: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        eos = str(r.get(eos_col, "")).strip()
        if not eos:
            # Try infer from a filename/path column.
            eos = "UNKNOWN"
        # Map partial names.
        mapped = None
        for name in args.eos:
            if name.lower() in eos.lower():
                mapped = name
                break
        eos = mapped or eos
        if args.eos and eos not in args.eos and eos != "UNKNOWN":
            continue
        groups.setdefault(eos, []).append(r)

    # If all unknown and requested EOS cannot be separated, stop transparently.
    if set(groups) == {"UNKNOWN"} and len(args.eos) > 1:
        status = "NS_1P4_REFINEMENT_FAILED_NO_EOS_COLUMN"
    else:
        status = "NS_1P4_REFINEMENT_DONE_INTERPOLATED"

    results: list[dict[str, Any]] = []
    per_column: list[dict[str, Any]] = []
    for eos in sorted(groups.keys(), key=lambda x: args.eos.index(x) if x in args.eos else 999):
        pts = groups[eos]
        try:
            res, pc = interpolate_eos(eos, pts, mass_col, numeric_cols, target)
            results.append(res)
            per_column.extend(pc)
        except Exception as exc:
            results.append({
                "EOS": eos,
                "target_M": target,
                "n_points": len(pts),
                "interp_mode": "FAILED",
                "error": repr(exc),
            })

    lambda_gr_col, lambda_bpb_col = identify_lambda_columns(cols)
    k2_gr_col, k2_bpb_col = identify_k2_columns(cols)
    r_col = radius_col

    compact_rows = []
    for r in results:
        cr = {
            "EOS": r.get("EOS"),
            "target_M": target,
            "n_points": r.get("n_points"),
            "mass_range": f"{fmt(r.get('mass_min'))}–{fmt(r.get('mass_max'))}",
            "interp_mode": r.get("interp_mode"),
            "M_lo": r.get("M_lo"),
            "M_hi": r.get("M_hi"),
            "bracket_width": r.get("bracket_width"),
        }
        if r_col:
            cr["R1p4_km"] = r.get(f"{norm_name(r_col)}_1p4")
        if lambda_gr_col:
            cr["Lambda_GR_1p4"] = r.get(f"{norm_name(lambda_gr_col)}_1p4")
        if lambda_bpb_col:
            cr["Lambda_BPB_or_ST_1p4"] = r.get(f"{norm_name(lambda_bpb_col)}_1p4")
        # If no pair, include first available lambda-ish numeric output.
        if not lambda_bpb_col:
            lambda_cols = [c for c in cols if "lambda" in norm_name(c)]
            if lambda_cols:
                c = lambda_cols[0]
                cr["Lambda_detected_1p4"] = r.get(f"{norm_name(c)}_1p4")
        if k2_gr_col:
            cr["k2_GR_1p4"] = r.get(f"{norm_name(k2_gr_col)}_1p4")
        if k2_bpb_col:
            cr["k2_BPB_or_ST_1p4"] = r.get(f"{norm_name(k2_bpb_col)}_1p4")
        # delta percent
        lg = cr.get("Lambda_GR_1p4")
        lb = cr.get("Lambda_BPB_or_ST_1p4")
        if lg not in (None, "", 0) and lb not in (None, ""):
            try:
                cr["delta_Lambda_over_Lambda_percent_1p4"] = 100.0 * (float(lb) / float(lg) - 1.0)
            except Exception:
                pass
        compact_rows.append(cr)

    decisions = [
        {
            "item": "scope",
            "decision": "INTERPOLATION_ONLY_NO_SOLVER_RERUN",
            "basis": "Uses completed strict-global rows/checkpoints; does not modify src/ns, checkpoint JSON, legacy or hybrid branches.",
        },
        {
            "item": "target_mass",
            "decision": f"{target:.6f}_MSUN",
            "basis": "Paper-table refinement target.",
        },
        {
            "item": "mass_column",
            "decision": mass_col or "NOT_FOUND",
            "basis": "Auto-detected from rows schema.",
        },
        {
            "item": "radius_column",
            "decision": radius_col or "NOT_FOUND",
            "basis": "Auto-detected from rows schema.",
        },
        {
            "item": "lambda_columns",
            "decision": f"GR={lambda_gr_col or 'NOT_FOUND'}; BPB/ST={lambda_bpb_col or 'NOT_FOUND'}",
            "basis": "Auto-detected from lambda-like column names; inspect per-column CSV if ambiguous.",
        },
        {
            "item": "native_k2_solver",
            "decision": "OPEN_NOT_TOUCHED",
            "basis": "This gate refines strict effective outputs only; it is not a native scalar-tensor k2 derivation.",
        },
    ]

    source_lock = list(source_rows)
    summary_csv = Path(args.summary_csv)
    source_lock.append({
        "path": str(summary_csv),
        "exists": summary_csv.exists(),
        "used": summary_csv.exists(),
        "kind": "freeze_summary_csv_context",
        "n_rows_loaded": len(read_csv_dicts(summary_csv)) if summary_csv.exists() else 0,
        "sha256": sha256(summary_csv),
    })

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_csv(out_prefix.with_name(out_prefix.name + "_compact.csv"), compact_rows)
    write_csv(out_prefix.with_name(out_prefix.name + "_all_interpolated_columns.csv"), results)
    write_csv(out_prefix.with_name(out_prefix.name + "_per_column_interp.csv"), per_column)
    write_csv(out_prefix.with_name(out_prefix.name + "_source_lock.csv"), source_lock)
    write_csv(out_prefix.with_name(out_prefix.name + "_decisions.csv"), decisions)

    payload = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target_mass": target,
        "detected_columns": {
            "eos_col": eos_col,
            "mass_col": mass_col,
            "radius_col": radius_col,
            "lambda_gr_col": lambda_gr_col,
            "lambda_bpb_or_st_col": lambda_bpb_col,
            "k2_gr_col": k2_gr_col,
            "k2_bpb_or_st_col": k2_bpb_col,
        },
        "n_input_rows": len(rows),
        "n_eos_results": len(results),
        "compact_rows": compact_rows,
        "decisions": decisions,
        "sources": source_lock,
    }
    out_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = []
    md.append("# NS 1.4 Msun refinement from strict-global gate\n")
    md.append("## Status\n")
    md.append(f"**{status}**\n")
    md.append("Non-destructive interpolation gate. It reads the completed strict-dictionary rows/checkpoints and does not run legacy/hybrid or modify `src/ns`.\n")
    md.append("## Detected schema\n")
    md.append("| role | column |\n|---|---|")
    md.append(f"| EOS | `{eos_col}` |")
    md.append(f"| mass | `{mass_col}` |")
    md.append(f"| radius | `{radius_col or ''}` |")
    md.append(f"| Lambda GR | `{lambda_gr_col or ''}` |")
    md.append(f"| Lambda BPB/ST | `{lambda_bpb_col or ''}` |")
    md.append(f"| k2 GR | `{k2_gr_col or ''}` |")
    md.append(f"| k2 BPB/ST | `{k2_bpb_col or ''}` |")
    md.append("\n## Compact 1.4-Msun table\n")
    if compact_rows:
        headers = list(compact_rows[0].keys())
        md.append("| " + " | ".join(headers) + " |")
        md.append("|" + "|".join(["---"] * len(headers)) + "|")
        for r in compact_rows:
            md.append("| " + " | ".join(fmt(r.get(h), 8 if "delta" in h.lower() else 6) for h in headers) + " |")
    else:
        md.append("No compact rows generated.")
    md.append("\n## Decisions\n")
    md.append("| item | decision | basis |\n|---|---|---|")
    for d in decisions:
        md.append(f"| {d['item']} | {d['decision']} | {d['basis']} |")
    md.append("\n## Guardrails\n")
    md.append("- This is an interpolation/refinement from the completed strict-global gate, not a native `k2_BPB` solver.")
    md.append("- If an EOS is marked `nearest_only_no_bracket`, the checkpoint rows do not bracket 1.4 Msun and a targeted solve is still required.")
    md.append("- If Lambda columns are ambiguous or missing, inspect `*_all_interpolated_columns.csv` and `*_per_column_interp.csv` before using the table in Paper IV.")

    out_prefix.with_suffix(".md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("Saved:", out_prefix.with_suffix(".md"))
    print("Saved:", out_prefix.with_suffix(".json"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_compact.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_all_interpolated_columns.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_per_column_interp.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_source_lock.csv"))
    print("Saved:", out_prefix.with_name(out_prefix.name + "_decisions.csv"))

if __name__ == "__main__":
    main()
