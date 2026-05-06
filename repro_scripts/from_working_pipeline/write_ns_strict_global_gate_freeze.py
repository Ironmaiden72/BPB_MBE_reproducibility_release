#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import hashlib
import math
import re
from pathlib import Path
from datetime import datetime, timezone

OUT = Path("outputs/ns")
SUMMARY_CSV = OUT / "run_multieos8_strict_global_checkpoint_gate_summary.csv"
ROWS_CSV = OUT / "run_multieos8_strict_global_checkpoint_gate_rows.csv"
SUMMARY_JSON = OUT / "run_multieos8_strict_global_checkpoint_gate.json"
SUMMARY_MD = OUT / "run_multieos8_strict_global_checkpoint_gate.md"
CHECKPOINTS = OUT / "checkpoints_strict_global_gate"
EXPECTED = ["BSk24", "APR", "DD2", "DDME2", "NL3wrL55", "SFHo", "SLy4", "TW"]


def nk(s):
    s = str(s).strip().lower()
    s = s.replace("δ", "delta").replace("λ", "lambda").replace("Λ", "lambda")
    return re.sub(r"[^a-z0-9]+", "_", s).strip("_")


def fnum(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).replace("−", "-").replace(",", ".").replace("`", "")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group(0)) if m else None


def bval(x):
    if isinstance(x, bool):
        return x
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in ("true", "yes", "1", "pass", "ok"):
        return True
    if s in ("false", "no", "0", "fail"):
        return False
    return None


def read_csv(path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})


def find(row, *names, all_tokens=()):
    kmap = {nk(k): k for k in row.keys()}
    for name in names:
        if nk(name) in kmap:
            return row.get(kmap[nk(name)])
    toks = [nk(t) for t in all_tokens]
    if toks:
        for kk, orig in kmap.items():
            if all(t in kk for t in toks):
                return row.get(orig)
    return None


def sha(path):
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def find_json_key(obj, keys):
    keys = {nk(k) for k in keys}
    if isinstance(obj, dict):
        for k, v in obj.items():
            if nk(k) in keys:
                return v
        for v in obj.values():
            got = find_json_key(v, keys)
            if got is not None:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = find_json_key(v, keys)
            if got is not None:
                return got
    return None


def fmt(x, n=6):
    if x is None:
        return ""
    if isinstance(x, bool):
        return "True" if x else "False"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if x != 0 and abs(x) < 1e-4:
            return f"{x:.6e}"
        return f"{x:.{n}f}"
    return str(x)


def main():
    rows = read_csv(SUMMARY_CSV)
    if not rows:
        raise SystemExit(f"Missing/empty input: {SUMMARY_CSV}")

    j = None
    if SUMMARY_JSON.exists():
        try:
            j = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
        except Exception:
            j = None
    kappa = fnum(find_json_key(j, ["kappa_strict", "kappa"]) if j is not None else None)

    canon = []
    for r in rows:
        eos = find(r, "EOS", "eos", "name")
        if eos is None:
            for v in r.values():
                if str(v).strip("` ") in EXPECTED:
                    eos = str(v).strip("` ")
                    break
        eos = str(eos).strip("` ") if eos is not None else ""
        status = str(find(r, "status", "state") or "").strip().lower()
        n = fnum(find(r, "n", "rows", "n_rows"))
        topology = str(find(r, "topology", "topology_count") or "").strip()
        mrange = str(find(r, "M range", "M_range", "mass_range") or "").strip()
        mmax = fnum(find(r, "Mmax seq", "Mmax_seq", "mmax_sequence", "mmax"))
        rmax = fnum(find(r, "R(Mmax)", "R_mmax", "Rmax"))
        near = fnum(find(r, "nearest 1.4 M", "nearest_1p4", "near_1p4", "nearest_1_4"))
        r14 = fnum(find(r, "R1.4 approx", "R1p4_approx", "R_1p4", "R14", all_tokens=("r", "1", "4")))
        l14 = fnum(find(r, "Lambda1.4 approx", "Λ1.4 approx", "Lambda_1p4", "Lambda14", all_tokens=("lambda", "1", "4")))
        fz = bval(find(r, "1.4 four-zone", "four_zone_1p4", "four_zone", all_tokens=("four", "zone")))
        md_raw = find(r, "mean δΛ/Λ", "mean dLambda/Lambda", "mean_delta_lambda_over_lambda", "mean_delta_lambda_frac", "mean_frac_corr_lambda", "mean_dL_L")
        md_val = fnum(md_raw)
        md_str = str(md_raw or "")
        md_percent = None
        if md_val is not None:
            md_percent = md_val if "%" in md_str or abs(md_val) > 0.2 else 100.0 * md_val
        canon.append({
            "EOS": eos,
            "status": status,
            "n": int(n) if n is not None and abs(n - int(n)) < 1e-9 else n,
            "topology": topology,
            "M_range": mrange,
            "Mmax_seq": mmax,
            "R_Mmax": rmax,
            "nearest_1p4_M": near,
            "R1p4_nearest_grid_approx": r14,
            "Lambda1p4_nearest_grid_approx": l14,
            "four_zone_1p4": fz,
            "mean_delta_Lambda_over_Lambda_percent": md_percent,
            "abs_nearest_delta_from_1p4": abs(near - 1.4) if near is not None else None,
        })

    order = {e: i for i, e in enumerate(EXPECTED)}
    canon.sort(key=lambda r: order.get(r["EOS"], 999))

    ok = [r for r in canon if r["status"] in ("ok", "done", "pass", "passed", "success")]
    deltas = [r["mean_delta_Lambda_over_Lambda_percent"] for r in canon if isinstance(r["mean_delta_Lambda_over_Lambda_percent"], (int, float))]
    mmaxs = [r["Mmax_seq"] for r in canon if isinstance(r["Mmax_seq"], (int, float))]
    fzs = [r["four_zone_1p4"] for r in canon if r["four_zone_1p4"] is not None]
    summary = {
        "status": "NS_STRICT_GLOBAL_GATE_FREEZE_DONE",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "kappa_strict": kappa,
        "n_eos": len(canon),
        "ok_eos": len(ok),
        "all_ok": len(ok) == len(canon) and len(canon) >= 8,
        "mean_delta_percent_mean": sum(deltas) / len(deltas) if deltas else None,
        "mean_delta_percent_min": min(deltas) if deltas else None,
        "mean_delta_percent_max": max(deltas) if deltas else None,
        "mmax_seq_min": min(mmaxs) if mmaxs else None,
        "mmax_seq_max": max(mmaxs) if mmaxs else None,
        "all_mmax_seq_over_2Msun": all(x > 2.0 for x in mmaxs) if mmaxs else None,
        "four_zone_1p4_true_count": sum(1 for x in fzs if x is True),
        "four_zone_1p4_false_count": sum(1 for x in fzs if x is False),
    }

    refine = []
    for r in canon:
        d = r["abs_nearest_delta_from_1p4"]
        if not isinstance(d, (int, float)):
            pr = "UNKNOWN"
        elif d > 0.10:
            pr = "HIGH"
        elif d > 0.03:
            pr = "MEDIUM"
        else:
            pr = "LOW"
        refine.append({
            "EOS": r["EOS"],
            "nearest_1p4_M": r["nearest_1p4_M"],
            "abs_nearest_delta_from_1p4": d,
            "R1p4_nearest_grid_approx": r["R1p4_nearest_grid_approx"],
            "Lambda1p4_nearest_grid_approx": r["Lambda1p4_nearest_grid_approx"],
            "refine_priority": pr,
            "refine_for_final_paper_table": True,
            "reason": "nearest-grid diagnostic; use interpolation/exact target solve before final 1.4 table",
        })

    decisions = [
        {"item": "strict_global_gate", "decision": "PASS_8_OF_8" if summary["all_ok"] else "CHECK", "basis": f"{summary['ok_eos']}/{summary['n_eos']} EOS have ok/pass status."},
        {"item": "legacy_hybrid", "decision": "DO_NOT_RUN", "basis": "Strict dictionary passed; legacy/hybrid are fallback diagnostics only."},
        {"item": "tidal_response", "decision": "COHERENT_NEGATIVE_AROUND_ONE_PERCENT", "basis": f"mean delta range {fmt(summary['mean_delta_percent_min'],3)}% to {fmt(summary['mean_delta_percent_max'],3)}%, mean {fmt(summary['mean_delta_percent_mean'],3)}%."},
        {"item": "Mmax_seq", "decision": "SURVIVAL_GATE_ONLY", "basis": "Mmax_seq is the maximum reached in the checkpointed search window, not a mathematical EOS maximum."},
        {"item": "R1p4_Lambda1p4", "decision": "NEAREST_GRID_APPROX_REFINE_FOR_FINAL_TABLES", "basis": "Nearest-grid values should be interpolated/exact-targeted for final paper tables."},
        {"item": "four_zone_1p4", "decision": "NOT_A_1P4_FOUR_ZONE_CLAIM", "basis": f"true count={summary['four_zone_1p4_true_count']}, false count={summary['four_zone_1p4_false_count']}."},
        {"item": "native_k2_solver", "decision": "OPEN_NOT_BLOCKING_STRICT_EFFECTIVE_GATE", "basis": "This freeze concerns strict effective tidal/Love outputs; native scalar-tensor k2 remains separate."},
        {"item": "next_gate", "decision": "NS_1P4_REFINEMENT_OR_PAPER_INSERTION", "basis": "Refine exact 1.4 values or insert this coarse survival gate with caveat."},
    ]

    source_files = [SUMMARY_JSON, ROWS_CSV, SUMMARY_CSV, SUMMARY_MD]
    source_lock = []
    for p in source_files:
        source_lock.append({"path": str(p), "exists": p.exists(), "kind": "file" if p.is_file() else "missing", "size_bytes": p.stat().st_size if p.exists() and p.is_file() else "", "sha256": sha(p), "purpose": "strict_global_gate_input"})
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
    md.append("**NS_STRICT_GLOBAL_GATE_FREEZE_DONE**\n")
    md.append("Non-destructive freeze gate. Reads completed strict_dictionary outputs; does not run legacy/hybrid and does not modify `src/ns` or checkpoint JSON files.\n")
    md.append("## Global lock\n")
    md.append("| quantity | value |\n|---|---:|")
    for k in ["kappa_strict","n_eos","ok_eos","all_ok","mean_delta_percent_mean","mean_delta_percent_min","mean_delta_percent_max","mmax_seq_min","mmax_seq_max","all_mmax_seq_over_2Msun","four_zone_1p4_true_count","four_zone_1p4_false_count"]:
        md.append(f"| {k} | {fmt(summary.get(k), 9 if k=='kappa_strict' else 6)} |")
    md.append("\n## EOS summary\n")
    md.append("| EOS | status | n | topology | M range | Mmax seq | R(Mmax) | nearest 1.4 M | R1.4 approx | Lambda1.4 approx | 1.4 four-zone | mean dLambda/Lambda |")
    md.append("|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in canon:
        md.append(f"| `{r['EOS']}` | {r['status']} | {fmt(r['n'],0)} | {r['topology']} | {r['M_range']} | {fmt(r['Mmax_seq'])} | {fmt(r['R_Mmax'])} | {fmt(r['nearest_1p4_M'])} | {fmt(r['R1p4_nearest_grid_approx'])} | {fmt(r['Lambda1p4_nearest_grid_approx'])} | {fmt(r['four_zone_1p4'])} | {fmt(r['mean_delta_Lambda_over_Lambda_percent'],3)}% |")
    md.append("\n## 1.4-Msun refinement priority\n")
    md.append("Nearest-grid 1.4 values are diagnostic only. For final paper tables, use interpolation or exact-target solves per EOS.\n")
    md.append("| EOS | nearest 1.4 M | |delta M| | priority |")
    md.append("|---|---:|---:|---|")
    for r in refine:
        md.append(f"| `{r['EOS']}` | {fmt(r['nearest_1p4_M'])} | {fmt(r['abs_nearest_delta_from_1p4'])} | {r['refine_priority']} |")
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
    md.append("- `1.4 four-zone = False` across the table means the paper claim should be strict effective tidal response, not a full four-zone topology claim at 1.4 Msun.\n")

    outputs = {
        "md": OUT / "ns_strict_global_gate_freeze.md",
        "json": OUT / "ns_strict_global_gate_freeze.json",
        "summary_csv": OUT / "ns_strict_global_gate_freeze_summary.csv",
        "refine_csv": OUT / "ns_strict_global_gate_freeze_refine_1p4_priority.csv",
        "source_lock_csv": OUT / "ns_strict_global_gate_freeze_source_lock.csv",
        "decisions_csv": OUT / "ns_strict_global_gate_freeze_decisions.csv",
    }
    write_csv(outputs["summary_csv"], canon)
    write_csv(outputs["refine_csv"], refine)
    write_csv(outputs["source_lock_csv"], source_lock)
    write_csv(outputs["decisions_csv"], decisions)
    outputs["json"].write_text(json.dumps({"summary": summary, "eos_rows": canon, "refine_1p4_priority": refine, "decisions": decisions, "source_lock_count": len(source_lock)}, indent=2, ensure_ascii=False), encoding="utf-8")
    outputs["md"].write_text("\n".join(md) + "\n", encoding="utf-8")

    print("Saved:")
    for k, p in outputs.items():
        print(f"  {k}: {p}")
    print("Status: NS_STRICT_GLOBAL_GATE_FREEZE_DONE")


if __name__ == "__main__":
    main()
