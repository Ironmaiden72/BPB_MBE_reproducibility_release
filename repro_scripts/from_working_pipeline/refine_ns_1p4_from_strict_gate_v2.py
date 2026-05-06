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

EOS_NAMES = ["BSk24", "APR", "DD2", "DDME2", "NL3wrL55", "SFHo", "SLy4", "TW"]

DEFAULT_SPARSE_ROWS_CSV = Path("outputs/ns/run_multieos8_strict_global_checkpoint_gate_rows.csv")
DEFAULT_CHECKPOINT_DIR = Path("outputs/ns/checkpoints_strict_global_gate")
DEFAULT_OUT_PREFIX = Path("outputs/ns/ns_1p4_refinement_from_strict_gate_v2")
TARGET_DEFAULT = 1.4

MASS_ALIASES = [
    "M", "M_msun", "M_Msun", "mass", "mass_msun", "M_solar", "Mgrav",
    "M_grav", "M_G", "m_msun", "gravitational_mass",
]
RADIUS_ALIASES = ["R", "R_km", "radius", "radius_km", "Rcirc", "Rcirc_km"]
EOS_ALIASES = ["EOS", "eos", "eos_name", "name", "template", "eos_label"]

LAMBDA_GR_HINTS = ["lambda_gr", "lambda_tov", "lambda_lcdm", "lambda_background", "lambda_base", "lambda_ref"]
LAMBDA_BPB_HINTS = ["lambda_bpb", "lambda_st", "lambda_strict", "lambda_eff", "lambda_bpb_eff", "lambda_modified"]

def norm_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower()).strip("_")

def sha256(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_float(x: Any) -> float | None:
    if x is None or isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        v = float(x)
        return v if math.isfinite(v) else None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "null", "true", "false"}:
        return None
    s = s.replace("−", "-").replace(",", ".")
    m = re.search(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        v = float(m.group(0))
    except Exception:
        return None
    return v if math.isfinite(v) else None

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

def all_columns(rows: list[dict[str, Any]]) -> list[str]:
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    return cols

def find_col(cols: list[str], aliases: list[str]) -> str | None:
    nmap = {norm_name(c): c for c in cols}
    for a in aliases:
        if norm_name(a) in nmap:
            return nmap[norm_name(a)]
    for a in aliases:
        na = norm_name(a)
        for nc, c in nmap.items():
            if na and (nc.endswith("_" + na) or nc.startswith(na + "_")):
                return c
    return None

def eos_from_text(s: str) -> str | None:
    for eos in EOS_NAMES:
        if eos.lower() in str(s).lower():
            return eos
    return None

def flatten_candidate_rows(obj: Any, eos_hint: str | None = None) -> list[dict[str, Any]]:
    """
    Recursively collect dict rows that look like physical sequence rows.
    This intentionally keeps broad candidates, then a later schema filter decides.
    """
    rows: list[dict[str, Any]] = []

    def rec(x: Any, ctx: dict[str, Any]):
        if isinstance(x, dict):
            ctx2 = dict(ctx)
            for k in EOS_ALIASES:
                if k in x and x[k]:
                    ctx2["eos"] = x[k]
            # A single dict can itself be a row.
            rows.append({**ctx2, **x})
            for v in x.values():
                if isinstance(v, (dict, list)):
                    rec(v, ctx2)
        elif isinstance(x, list):
            for item in x:
                rec(item, ctx)

    rec(obj, {"eos": eos_hint} if eos_hint else {})
    return rows

def schema_score(rows: list[dict[str, Any]]) -> dict[str, Any]:
    cols = all_columns(rows)
    mcol = find_col(cols, MASS_ALIASES)
    rcol = find_col(cols, RADIUS_ALIASES)
    ecol = find_col(cols, EOS_ALIASES)
    lam_cols = [c for c in cols if "lambda" in norm_name(c)]
    if not mcol:
        return {"usable": False, "score": 0, "mass_col": None, "radius_col": rcol, "eos_col": ecol, "lambda_cols": lam_cols, "usable_rows": 0}
    usable_rows = 0
    for r in rows:
        if parse_float(r.get(mcol)) is not None:
            # require at least one physics-like extra numeric value
            if rcol and parse_float(r.get(rcol)) is not None:
                usable_rows += 1
            elif any(parse_float(r.get(c)) is not None for c in lam_cols):
                usable_rows += 1
    score = usable_rows + 10 * int(bool(rcol)) + 10 * len(lam_cols)
    return {"usable": usable_rows > 0, "score": score, "mass_col": mcol, "radius_col": rcol, "eos_col": ecol, "lambda_cols": lam_cols, "usable_rows": usable_rows}

def load_sources(sparse_csv: Path, checkpoint_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """
    Load sparse CSV and checkpoint JSONs. Choose the source with the most usable rows.
    This prevents accidental interpolation from the sparse summary CSV when dense
    checkpoint JSON rows exist.
    """
    candidates: list[tuple[str, list[dict[str, Any]], dict[str, Any], Path | None]] = []
    source_lock: list[dict[str, Any]] = []

    csv_rows = read_csv_dicts(sparse_csv)
    sc = schema_score(csv_rows)
    candidates.append(("sparse_rows_csv", csv_rows, sc, sparse_csv))
    source_lock.append({
        "path": str(sparse_csv),
        "kind": "sparse_rows_csv",
        "exists": sparse_csv.exists(),
        "n_candidate_rows": len(csv_rows),
        "usable_rows": sc.get("usable_rows", 0),
        "sha256": sha256(sparse_csv),
    })

    merged_json_rows: list[dict[str, Any]] = []
    if checkpoint_dir.exists():
        for jp in sorted(checkpoint_dir.glob("*.json")):
            eos_hint = eos_from_text(jp.name)
            try:
                obj = json.loads(jp.read_text(encoding="utf-8"))
                rows = flatten_candidate_rows(obj, eos_hint=eos_hint)
                scj = schema_score(rows)
                merged_json_rows.extend(rows)
                source_lock.append({
                    "path": str(jp),
                    "kind": "checkpoint_json",
                    "exists": True,
                    "n_candidate_rows": len(rows),
                    "usable_rows": scj.get("usable_rows", 0),
                    "sha256": sha256(jp),
                })
            except Exception as exc:
                source_lock.append({
                    "path": str(jp),
                    "kind": "checkpoint_json_error",
                    "exists": True,
                    "error": repr(exc),
                    "n_candidate_rows": 0,
                    "usable_rows": 0,
                    "sha256": sha256(jp),
                })
    else:
        source_lock.append({
            "path": str(checkpoint_dir),
            "kind": "checkpoint_dir",
            "exists": False,
            "n_candidate_rows": 0,
            "usable_rows": 0,
            "sha256": "",
        })

    sc_json = schema_score(merged_json_rows)
    candidates.append(("merged_checkpoint_json", merged_json_rows, sc_json, checkpoint_dir))

    # Pick best usable source by score; if JSON has more usable rows, it wins.
    usable = [c for c in candidates if c[2].get("usable")]
    if not usable:
        return [], source_lock, {"selected_source": "NONE", "status": "NO_USABLE_SOURCE"}
    selected = max(usable, key=lambda c: c[2]["score"])
    selected_name, rows, schema, path = selected

    for s in source_lock:
        s["selected"] = (selected_name in s["kind"]) or (selected_name == "merged_checkpoint_json" and s["kind"].startswith("checkpoint_json"))

    meta = {
        "selected_source": selected_name,
        "selected_path": str(path) if path else "",
        "schema": schema,
        "candidate_summaries": [
            {"name": name, "path": str(p) if p else "", **schema_i}
            for name, _, schema_i, p in candidates
        ],
    }
    return rows, source_lock, meta

def filter_physical_rows(rows: list[dict[str, Any]], mass_col: str, radius_col: str | None, lambda_cols: list[str]) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        m = parse_float(r.get(mass_col))
        if m is None:
            continue
        has_radius = radius_col and parse_float(r.get(radius_col)) is not None
        has_lambda = any(parse_float(r.get(c)) is not None for c in lambda_cols)
        if has_radius or has_lambda:
            out.append(r)
    return out

def pick_numeric_cols(rows: list[dict[str, Any]]) -> list[str]:
    cols = all_columns(rows)
    nums = []
    for c in cols:
        if sum(parse_float(r.get(c)) is not None for r in rows) >= 2:
            nums.append(c)
    return nums

def identify_lambda_columns(cols: list[str]) -> tuple[str | None, str | None]:
    lambda_cols = [c for c in cols if "lambda" in norm_name(c)]
    if not lambda_cols:
        return None, None

    def score(c: str, hints: list[str]) -> int:
        nc = norm_name(c)
        s = 0
        for h in hints:
            if norm_name(h) in nc:
                s += 20
        if any(t in nc for t in ["gr", "lcdm", "base", "ref", "tov"]):
            s += 3 if hints is LAMBDA_GR_HINTS else -2
        if any(t in nc for t in ["bpb", "st", "strict", "eff", "modified"]):
            s += 3 if hints is LAMBDA_BPB_HINTS else -2
        return s

    gr = max(lambda_cols, key=lambda c: score(c, LAMBDA_GR_HINTS))
    bpb = max(lambda_cols, key=lambda c: score(c, LAMBDA_BPB_HINTS))
    if score(gr, LAMBDA_GR_HINTS) <= 0:
        gr = None
    if score(bpb, LAMBDA_BPB_HINTS) <= 0:
        bpb = None
    if gr == bpb:
        # Keep one detected Lambda but no delta pair.
        bpb = gr
        gr = None
    return gr, bpb

def unique_by_mass(points: list[dict[str, Any]], mass_col: str, numeric_cols: list[str]) -> list[dict[str, Any]]:
    grouped: dict[float, list[dict[str, Any]]] = {}
    for r in points:
        m = parse_float(r.get(mass_col))
        if m is None:
            continue
        grouped.setdefault(m, []).append(r)

    out = []
    for m, group in grouped.items():
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

def interp(x0: float, y0: float, x1: float, y1: float, x: float) -> float:
    if abs(x1 - x0) < 1e-14:
        return y0
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

def local_interp(points: list[dict[str, Any]], mass_col: str, numeric_cols: list[str], target: float) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pts = unique_by_mass(points, mass_col, numeric_cols)
    below = [p for p in pts if parse_float(p.get(mass_col)) is not None and parse_float(p.get(mass_col)) <= target]
    above = [p for p in pts if parse_float(p.get(mass_col)) is not None and parse_float(p.get(mass_col)) >= target]

    if not pts:
        raise ValueError("no unique mass points")

    if below and above:
        lo = max(below, key=lambda p: parse_float(p.get(mass_col)))
        hi = min(above, key=lambda p: parse_float(p.get(mass_col)))
        mode = "exact" if lo is hi or abs(parse_float(lo.get(mass_col)) - parse_float(hi.get(mass_col))) < 1e-12 else "bracketed_linear"
    else:
        lo = hi = min(pts, key=lambda p: abs(parse_float(p.get(mass_col)) - target))
        mode = "nearest_only_no_bracket"

    mlo = parse_float(lo.get(mass_col))
    mhi = parse_float(hi.get(mass_col))

    out: dict[str, Any] = {
        "n_points": len(pts),
        "mass_min": min(parse_float(p.get(mass_col)) for p in pts),
        "mass_max": max(parse_float(p.get(mass_col)) for p in pts),
        "interp_mode": mode,
        "M_lo": mlo,
        "M_hi": mhi,
        "bracket_width": abs(mhi - mlo) if mlo is not None and mhi is not None else "",
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
            val = interp(mlo, vlo, mhi, vhi, target)
        key = f"{norm_name(c)}_1p4"
        out[key] = val
        per_col.append({
            "column": c,
            "output_key": key,
            "interp_mode": mode,
            "M_lo": mlo,
            "value_lo": vlo,
            "M_hi": mhi,
            "value_hi": vhi,
            "value_1p4": val,
        })
    return out, per_col

def fmt(x: Any, n: int = 6) -> str:
    if x is None or x == "":
        return ""
    if isinstance(x, bool):
        return "True" if x else "False"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.{n}f}"
    return str(x)

def main():
    ap = argparse.ArgumentParser(description="Dense/source-aware NS 1.4 interpolation from strict-global checkpoints.")
    ap.add_argument("--sparse-rows-csv", default=str(DEFAULT_SPARSE_ROWS_CSV))
    ap.add_argument("--checkpoint-dir", default=str(DEFAULT_CHECKPOINT_DIR))
    ap.add_argument("--target-mass", type=float, default=TARGET_DEFAULT)
    ap.add_argument("--out-prefix", default=str(DEFAULT_OUT_PREFIX))
    ap.add_argument("--max-final-bracket-width", type=float, default=0.03, help="If bracket wider than this Msun, mark as coarse-not-final.")
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    rows_raw, source_lock, meta = load_sources(Path(args.sparse_rows_csv), Path(args.checkpoint_dir))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not rows_raw:
        payload = {
            "status": "NS_1P4_REFINEMENT_V2_FAILED_NO_USABLE_SOURCE",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "source_meta": meta,
            "source_lock": source_lock,
        }
        out_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
        raise SystemExit("No usable source found.")

    schema = meta["schema"]
    mass_col = schema["mass_col"]
    radius_col = schema.get("radius_col")
    eos_col = schema.get("eos_col") or "eos"
    lambda_cols = schema.get("lambda_cols") or []
    rows = filter_physical_rows(rows_raw, mass_col, radius_col, lambda_cols)
    cols = all_columns(rows)
    numeric_cols = pick_numeric_cols(rows)
    lambda_gr_col, lambda_bpb_col = identify_lambda_columns(cols)

    groups: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        eos_val = r.get(eos_col) or r.get("eos") or r.get("EOS") or ""
        eos = eos_from_text(str(eos_val)) or eos_from_text(json.dumps(r)[:300]) or "UNKNOWN"
        groups.setdefault(eos, []).append(r)

    results = []
    per_col_all = []
    for eos in EOS_NAMES:
        pts = groups.get(eos, [])
        if not pts:
            results.append({"EOS": eos, "status": "NO_ROWS_FOR_EOS"})
            continue
        try:
            res, pc = local_interp(pts, mass_col, numeric_cols, args.target_mass)
            res["EOS"] = eos
            res["status"] = "OK"
            for p in pc:
                p["EOS"] = eos
            per_col_all.extend(pc)
            results.append(res)
        except Exception as exc:
            results.append({"EOS": eos, "status": "FAILED", "error": repr(exc), "n_raw_points": len(pts)})

    compact = []
    for r in results:
        cr = {
            "EOS": r.get("EOS"),
            "status": r.get("status"),
            "target_M": args.target_mass,
            "n_points": r.get("n_points", ""),
            "mass_range": f"{fmt(r.get('mass_min'))}–{fmt(r.get('mass_max'))}" if r.get("mass_min") not in (None, "") else "",
            "interp_mode": r.get("interp_mode", ""),
            "M_lo": r.get("M_lo", ""),
            "M_hi": r.get("M_hi", ""),
            "bracket_width": r.get("bracket_width", ""),
            "final_table_quality": "",
        }
        bw = parse_float(r.get("bracket_width"))
        if r.get("interp_mode") == "bracketed_linear":
            cr["final_table_quality"] = "FINAL_READY_LOCAL_INTERP" if bw is not None and bw <= args.max_final_bracket_width else "COARSE_BRACKET_REQUIRES_TARGETED_SOLVE"
        elif r.get("interp_mode") == "exact":
            cr["final_table_quality"] = "FINAL_READY_EXACT_GRID_POINT"
        elif r.get("interp_mode"):
            cr["final_table_quality"] = "NOT_FINAL"

        if radius_col:
            cr["R1p4_km"] = r.get(f"{norm_name(radius_col)}_1p4", "")
        if lambda_gr_col:
            cr["Lambda_GR_1p4"] = r.get(f"{norm_name(lambda_gr_col)}_1p4", "")
        if lambda_bpb_col:
            cr["Lambda_BPB_or_ST_1p4"] = r.get(f"{norm_name(lambda_bpb_col)}_1p4", "")
        lg = parse_float(cr.get("Lambda_GR_1p4"))
        lb = parse_float(cr.get("Lambda_BPB_or_ST_1p4"))
        if lg not in (None, 0) and lb is not None:
            cr["delta_Lambda_over_Lambda_percent_1p4"] = 100.0 * (lb / lg - 1.0)
        compact.append(cr)

    n_ready = sum(1 for r in compact if str(r.get("final_table_quality")).startswith("FINAL_READY"))
    status = "NS_1P4_REFINEMENT_V2_DONE"
    if n_ready < len(EOS_NAMES):
        status += "_COARSE_FOR_SOME_EOS_NOT_FINAL"

    decisions = [
        {
            "item": "source_selection",
            "decision": meta.get("selected_source"),
            "basis": "V2 compares sparse rows CSV and merged checkpoint JSON rows, then selects the source with the most usable physical rows.",
        },
        {
            "item": "paper_table_status",
            "decision": "FINAL_READY" if n_ready == len(EOS_NAMES) else "NOT_FINAL_FOR_ALL_EOS",
            "basis": f"{n_ready}/8 EOS have local bracket width <= {args.max_final_bracket_width} Msun or exact target.",
        },
        {
            "item": "coarse_v1_warning",
            "decision": "V1_VALUES_ARE_COARSE_IF_NPOINTS_MATCH_SPARSE_CSV",
            "basis": "The previous run used the sparse CSV when n_points were 30-61 instead of the full checkpoint sequence counts; wide brackets make Lambda interpolation non-final.",
        },
        {
            "item": "native_k2_solver",
            "decision": "OPEN_NOT_TOUCHED",
            "basis": "This remains strict effective interpolation, not a native scalar-tensor Love-number derivation.",
        },
    ]

    write_csv(out_prefix.with_name(out_prefix.name + "_compact.csv"), compact)
    write_csv(out_prefix.with_name(out_prefix.name + "_all_interpolated_columns.csv"), results)
    write_csv(out_prefix.with_name(out_prefix.name + "_per_column_interp.csv"), per_col_all)
    write_csv(out_prefix.with_name(out_prefix.name + "_source_lock.csv"), source_lock)
    write_csv(out_prefix.with_name(out_prefix.name + "_decisions.csv"), decisions)

    payload = {
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target_mass": args.target_mass,
        "source_meta": meta,
        "detected_columns": {
            "eos_col": eos_col,
            "mass_col": mass_col,
            "radius_col": radius_col,
            "lambda_gr_col": lambda_gr_col,
            "lambda_bpb_or_st_col": lambda_bpb_col,
        },
        "n_raw_rows_selected": len(rows_raw),
        "n_physical_rows_selected": len(rows),
        "n_final_ready": n_ready,
        "compact_rows": compact,
        "decisions": decisions,
        "source_lock": source_lock,
    }
    out_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = []
    md.append("# NS 1.4 Msun refinement from strict-global gate — v2 source-aware\n")
    md.append("## Status\n")
    md.append(f"**{status}**\n")
    md.append("Non-destructive interpolation gate. V2 prevents accidental use of sparse summary rows when denser checkpoint JSON data are available.\n")
    md.append("## Source selection\n")
    md.append(f"- selected source: `{meta.get('selected_source')}`")
    md.append(f"- selected path: `{meta.get('selected_path')}`")
    md.append(f"- selected usable rows: `{schema.get('usable_rows')}`")
    md.append("\n## Detected schema\n")
    md.append("| role | column |\n|---|---|")
    md.append(f"| EOS | `{eos_col}` |")
    md.append(f"| mass | `{mass_col}` |")
    md.append(f"| radius | `{radius_col or ''}` |")
    md.append(f"| Lambda GR | `{lambda_gr_col or ''}` |")
    md.append(f"| Lambda BPB/ST | `{lambda_bpb_col or ''}` |")
    md.append("\n## Compact 1.4-Msun table\n")
    headers = list(compact[0].keys()) if compact else []
    if headers:
        md.append("| " + " | ".join(headers) + " |")
        md.append("|" + "|".join(["---"] * len(headers)) + "|")
        for r in compact:
            md.append("| " + " | ".join(fmt(r.get(h), 8 if "delta" in h.lower() else 6) for h in headers) + " |")
    md.append("\n## Decisions\n")
    md.append("| item | decision | basis |\n|---|---|---|")
    for d in decisions:
        md.append(f"| {d['item']} | {d['decision']} | {d['basis']} |")
    md.append("\n## Guardrails\n")
    md.append("- Rows marked `COARSE_BRACKET_REQUIRES_TARGETED_SOLVE` must not be used as final Paper IV 1.4 tables.")
    md.append("- Wide brackets are especially unsafe for Lambda because Lambda(M) is highly nonlinear.")
    md.append("- If v2 still selects sparse CSV, the checkpoint JSON schema did not expose denser physical rows to this generic parser; use a native targeted 1.4 solve script instead.")
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
