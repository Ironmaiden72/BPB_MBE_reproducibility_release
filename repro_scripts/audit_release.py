#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release"
SUMS = RELEASE / "release_index" / "00_PACKAGE_SHA256SUMS.csv"
PDF_AUDIT = RELEASE / "audit" / "00_PDF_AUDIT_TABLE.csv"
MANIFEST = RELEASE / "release_index" / "00_MANIFEST_CORPUS.csv"

EXPECTED_PAPERS = {
    "Paper_0": "v0.9",
    "Paper_I": "v0.4.1",
    "Paper_II": "v0.4.1",
    "Paper_III": "v0.4.1",
    "Paper_IV": "v0.4.1",
    "Paper_V": "v0.3.2",
    "Paper_VI": "v0.2.1",
    "Paper_VII": "v0.3",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_checksums() -> list[str]:
    errors = []
    if not SUMS.exists():
        return [f"Missing checksum table: {SUMS}"]
    with SUMS.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        rel = row["path"]
        expected_bytes = int(row["bytes"])
        expected_sha = row["sha256"]
        path = RELEASE / rel
        if not path.exists():
            errors.append(f"MISSING {rel}")
            continue
        actual_bytes = path.stat().st_size
        actual_sha = sha256(path)
        if actual_bytes != expected_bytes:
            errors.append(f"SIZE {rel}: expected {expected_bytes}, got {actual_bytes}")
        if actual_sha != expected_sha:
            errors.append(f"SHA {rel}: expected {expected_sha}, got {actual_sha}")
    return errors


def audit_manifest() -> list[str]:
    errors = []
    if not MANIFEST.exists():
        return [f"Missing manifest: {MANIFEST}"]
    with MANIFEST.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    seen = {r["paper_id"]: r for r in rows}
    for paper_id, version in EXPECTED_PAPERS.items():
        if paper_id not in seen:
            errors.append(f"MANIFEST missing {paper_id}")
            continue
        if seen[paper_id]["version"] != version:
            errors.append(f"VERSION {paper_id}: expected {version}, got {seen[paper_id]['version']}")
        for key in ["pdf_path", "source_zip_path"]:
            p = RELEASE / seen[paper_id][key]
            if not p.exists():
                errors.append(f"MANIFEST file missing for {paper_id}: {key}={p}")
    return errors


def audit_pdf_table() -> list[str]:
    errors = []
    if not PDF_AUDIT.exists():
        return [f"Missing PDF audit table: {PDF_AUDIT}"]
    with PDF_AUDIT.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 8:
        errors.append(f"PDF audit row count expected 8, got {len(rows)}")
    for row in rows:
        status = row.get("status") or row.get("audit_status") or ""
        if status and status.upper() not in {"PASS", "OK"}:
            errors.append(f"PDF audit non-pass row: {row}")
    return errors


def main() -> None:
    errors = []
    errors += audit_checksums()
    errors += audit_manifest()
    errors += audit_pdf_table()
    if errors:
        print("CHECKSUM_AUDIT: FAIL")
        for e in errors:
            print(" -", e)
        raise SystemExit(1)
    print("CHECKSUM_AUDIT: PASS")
    print("PAPERS_EXPECTED:", len(EXPECTED_PAPERS))
    print("RELEASE_DIR:", RELEASE)

if __name__ == "__main__":
    main()
