from pathlib import Path
import shutil, hashlib, json, csv, textwrap, zipfile, os
from datetime import date

BASE = Path('/mnt/data')
SRC_RELEASE = BASE / 'BPB_MBE_HAL_CORPUS_FINAL_20260506'
if not SRC_RELEASE.exists():
    raise SystemExit(f'Missing {SRC_RELEASE}')

REPO = BASE / 'BPB_MBE_reproducibility_release_v2026_05_06'
if REPO.exists():
    shutil.rmtree(REPO)
REPO.mkdir()

# Copy final audited release as the payload of the repo.
release_dst = REPO / 'release'
shutil.copytree(SRC_RELEASE, release_dst)

# Convenience top-level copies of root release files.
for name in ['00_README_CORPUS.md','00_README_CORPUS.pdf']:
    src = release_dst / name
    if src.exists(): shutil.copy2(src, REPO / name)
for name in ['00_MANIFEST_CORPUS.csv','00_CLAIMS_INDEX.csv','00_REVIEWER_ROUTE.md','00_REVIEWER_ROUTE.csv','00_RELEASE_NOTES.md']:
    src = release_dst / 'release_index' / name
    if src.exists(): shutil.copy2(src, REPO / name)
for name in ['00_FINAL_RELEASE_AUDIT.md','00_PDF_AUDIT_TABLE.csv']:
    src = release_dst / 'audit' / name
    if src.exists(): shutil.copy2(src, REPO / name)

# Copy license and minimal requirements.
for name in ['LICENSE','requirements.txt']:
    src = BASE / name
    if src.exists(): shutil.copy2(src, REPO / name)

# Directories.
for d in ['repro_scripts','docs','data_instructions','repro_levels','.github/workflows']:
    (REPO / d).mkdir(parents=True, exist_ok=True)

# README
readme = '''# BPB/MBE reproducibility release capsule

This repository is a **clean reproducibility capsule** for the BPB/MBE HAL corpus release dated 2026-05-06.

It is not the private working pipeline and it is not a full laboratory history. It contains the audited release artefacts needed to inspect the submitted corpus: final PDFs, LaTeX source ZIPs, release indices, claim-source locks, checksums, and lightweight audit scripts.

## What is included

- `release/` — the final audited HAL corpus package.
- `release/pdfs/` — the eight final paper PDFs.
- `release/sources/` — source-only LaTeX/source ZIPs for each paper.
- `release/release_index/` — manifest, claims index, reviewer route, SHA256 index.
- `release/audit/` — automated PDF/release audit outputs.
- `repro_scripts/` — lightweight scripts to verify checksums and summarize claims.
- `data_instructions/` — notes about external or heavy data that are not redistributed here.
- `repro_levels/` — explanation of what can be checked at L0–L3.

## What is intentionally excluded

The private working pipeline contains historical experiments, abandoned branches, debug scripts, intermediate heavy outputs, and exploratory runs. Those are intentionally excluded from this public capsule.

The release is designed to answer the question:

> Are the submitted papers, numerical locks, source packages, and claim-source maps internally consistent and reproducible from the frozen release artefacts?

It is not designed to rerun every exploratory MCMC, CAMB, ACT, SPARC, or neutron-star production run from scratch.

## Quick audit

From the repository root:

```bash
python repro_scripts/audit_release.py
python repro_scripts/claim_summary.py
```

Expected headline result:

```text
CHECKSUM_AUDIT: PASS
CLAIM_ROWS: 107
```

## Reproducibility levels

| Level | Purpose | Included here |
|---|---|---:|
| L0 | Verify file presence, SHA256 checksums, paper versions, release manifest | Yes |
| L1 | Rebuild manuscript tables/figures from frozen CSV/source ZIPs where included | Partially, via per-paper source packages |
| L2 | Rerun lightweight gate diagnostics from frozen source artefacts | Partially; future expansion module |
| L3 | Rerun all heavy production/exploratory pipelines from raw external data | No; documented as future work |

## Citation

Use `CITATION.cff` once this repository is pushed to GitHub and archived by Zenodo. The Zenodo DOI can then be inserted into the HAL deposits.

## Release status

See:

- `00_MANIFEST_CORPUS.csv`
- `00_CLAIMS_INDEX.csv`
- `00_REVIEWER_ROUTE.md`
- `00_FINAL_RELEASE_AUDIT.md`
'''
(REPO/'README.md').write_text(readme, encoding='utf-8')

# CITATION.cff
citation = '''cff-version: 1.2.0
message: "If you use this reproducibility capsule, please cite the archived Zenodo release and the associated HAL corpus papers."
title: "BPB/MBE HAL corpus reproducibility release capsule"
authors:
  - family-names: "Pelletier"
    given-names: "F."
    affiliation: "Independent Research"
version: "v2026.05.06"
date-released: "2026-05-06"
repository-code: "https://github.com/<OWNER>/BPB_MBE_reproducibility_release"
license: "MIT"
keywords:
  - cosmology
  - bigravity
  - reproducibility
  - HAL
  - BPB/MBE
abstract: "Clean reproducibility capsule for the BPB/MBE HAL corpus release, containing final papers, source packages, claim-source locks, release indices and audit scripts."
'''
(REPO/'CITATION.cff').write_text(citation, encoding='utf-8')

# Zenodo draft
zenodo = {
    "title": "BPB/MBE HAL corpus reproducibility release capsule",
    "creators": [{"name": "Pelletier, F.", "affiliation": "Independent Research"}],
    "description": "Clean reproducibility capsule for the BPB/MBE HAL corpus release. Contains final papers, source packages, claim-source locks, release indices and audit scripts. Heavy private working-pipeline history is intentionally excluded.",
    "license": "MIT",
    "upload_type": "software",
    "keywords": ["cosmology", "bigravity", "reproducibility", "HAL", "BPB/MBE"],
    "version": "v2026.05.06"
}
(REPO/'.zenodo.json').write_text(json.dumps(zenodo, indent=2, ensure_ascii=False), encoding='utf-8')

# GitHub Actions optional audit.
workflow = '''name: release-audit
on:
  workflow_dispatch:
  push:
    paths:
      - 'release/**'
      - 'repro_scripts/**'
      - '00_*.csv'
      - '00_*.md'

jobs:
  audit:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Run release audit
        run: python repro_scripts/audit_release.py
      - name: Summarize claims
        run: python repro_scripts/claim_summary.py
'''
(REPO/'.github/workflows/release-audit.yml').write_text(workflow, encoding='utf-8')

# Scripts
script_audit = r'''#!/usr/bin/env python3
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
'''
(REPO/'repro_scripts/audit_release.py').write_text(script_audit, encoding='utf-8')
os.chmod(REPO/'repro_scripts/audit_release.py', 0o755)

script_claim = r'''#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLAIMS = ROOT / "release" / "release_index" / "00_CLAIMS_INDEX.csv"


def main() -> None:
    with CLAIMS.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    print(f"CLAIM_ROWS: {len(rows)}")
    # Try common column names robustly.
    keys = rows[0].keys() if rows else []
    paper_key = "paper_id" if "paper_id" in keys else ("paper" if "paper" in keys else None)
    level_key = "level" if "level" in keys else ("claim_level" if "claim_level" in keys else None)
    if paper_key:
        print("CLAIMS_BY_PAPER:")
        for paper, n in sorted(Counter(r[paper_key] for r in rows).items()):
            print(f"  {paper}: {n}")
    if level_key:
        print("CLAIMS_BY_LEVEL:")
        for level, n in sorted(Counter(r[level_key] for r in rows).items()):
            print(f"  {level}: {n}")

if __name__ == "__main__":
    main()
'''
(REPO/'repro_scripts/claim_summary.py').write_text(script_claim, encoding='utf-8')
os.chmod(REPO/'repro_scripts/claim_summary.py', 0o755)

# A helper to create a local GitHub/Zenodo release payload zip.
script_pack = r'''#!/usr/bin/env python3
from __future__ import annotations
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT.parent / f"{ROOT.name}.zip"
EXCLUDE_DIRS = {".git", "__pycache__"}

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
    for p in ROOT.rglob("*"):
        if any(part in EXCLUDE_DIRS for part in p.parts):
            continue
        if p.is_file():
            zf.write(p, p.relative_to(ROOT.parent))
print(OUT)
'''
(REPO/'repro_scripts/package_repo.py').write_text(script_pack, encoding='utf-8')
os.chmod(REPO/'repro_scripts/package_repo.py', 0o755)

# Repro level docs.
(REPO/'repro_levels/L0_checksum_audit.md').write_text('''# L0 checksum and release-index audit

Run:

```bash
python repro_scripts/audit_release.py
python repro_scripts/claim_summary.py
```

This verifies the final release payload against `release/release_index/00_PACKAGE_SHA256SUMS.csv`, checks that the manifest contains the expected paper versions, and summarizes the claim-source index.
''', encoding='utf-8')
(REPO/'repro_levels/L1_tables_figures.md').write_text('''# L1 tables and figures from frozen sources

The paper source ZIPs under `release/sources/` contain the source material used for the corpus manuscripts. This capsule freezes those source packages and does not unpack them by default.

A future L1 expansion can unpack each source ZIP and run its local figure/table build scripts, where included. The HAL-ready corpus does not require rerunning heavy production jobs for L1.
''', encoding='utf-8')
(REPO/'repro_levels/L2_light_diagnostics.md').write_text('''# L2 lightweight diagnostics

Some diagnostic gates are lightweight enough to rerun from frozen CSV/JSON artefacts. This public capsule currently prioritizes L0 release verification and claim indexing. L2 modules can be added incrementally without exposing the private working pipeline history.
''', encoding='utf-8')
(REPO/'repro_levels/L3_heavy_runs.md').write_text('''# L3 heavy production runs

Full production regeneration involves external data and heavy workflows: MCMC scans, CAMB/Planck-lite/CMB diagnostics, ACT likelihood work, SPARC raw-rotation-curve processing, and neutron-star EOS/tidal response runs.

These are documented by the claim-source locks and frozen artefacts. They are intentionally not bundled as a one-click public rerun in this clean capsule.
''', encoding='utf-8')

# Data instructions.
(REPO/'data_instructions/external_datasets.md').write_text('''# External data policy

This repository does not redistribute external datasets beyond the frozen release artefacts already packaged in the public corpus source bundles.

Relevant external families include Pantheon+/SNe, DESI BAO, KiDS, ACT DR6 lensing, SPARC mass models, EOS tables, and gravitational-wave/neutron-star reference data. Consult each original dataset's license and citation policy before attempting full L3 regeneration.
''', encoding='utf-8')
(REPO/'data_instructions/repo_scope.md').write_text('''# Scope of this clean reproduction repository

This is a release capsule. It excludes:

- private exploratory history;
- abandoned branches;
- debug scripts;
- very large intermediate products;
- local machine-specific paths;
- uncurated work-in-progress notebooks.

It includes exactly the material needed to audit the HAL corpus release and verify source-lock consistency.
''', encoding='utf-8')

# HAL metadata draft.
manifest_rows = []
with (release_dst/'release_index/00_MANIFEST_CORPUS.csv').open(newline='', encoding='utf-8') as f:
    manifest_rows = list(csv.DictReader(f))

def fmt_dep(row):
    return f'''## {row['paper_label']} — {row['title_short']}

- Version: `{row['version']}`
- Suggested HAL type: preprint / working paper / report
- Main file: `release/{row['pdf_path']}`
- Source attachment: `release/{row['source_zip_path']}`
- Audit status: `{row['audit_status']}`
- Suggested keywords: BPB/MBE; cosmology; bigravity; reproducibility; {row['role']}
- Comment to moderator: This file is part of the BPB/MBE HAL corpus release. The reproducibility capsule and claim-source index are provided as companion files.
'''
meta = '''# HAL deposit metadata draft

This is a working checklist for HAL deposits. Choose the final license deliberately in the HAL UI for each file. The repository-level code license is MIT, but manuscript licensing may be chosen separately.

Recommended corpus strategy:

1. Deposit `00_README_CORPUS.pdf` as the corpus guide.
2. Deposit Papers 0–VII as separate HAL records.
3. Attach source ZIPs to each paper record.
4. Put the GitHub/Zenodo DOI for this clean reproducibility capsule in the metadata or comments once available.

'''+"\n".join(fmt_dep(r) for r in manifest_rows)
(REPO/'docs/HAL_DEPOSIT_METADATA_DRAFT.md').write_text(meta, encoding='utf-8')

# Repo status.
status = f'''# Release status

Generated: 2026-05-06

Clean repo name: `BPB_MBE_reproducibility_release_v2026_05_06`

Final corpus id: `BPB_MBE_HAL_CORPUS_FINAL_20260506`

Expected paper versions:

| Paper | Version |
|---|---:|
| Paper 0 | v0.9 |
| Paper I | v0.4.1 |
| Paper II | v0.4.1 |
| Paper III | v0.4.1 |
| Paper IV | v0.4.1 |
| Paper V | v0.3.2 |
| Paper VI | v0.2.1 |
| Paper VII | v0.3 |

Run `python repro_scripts/audit_release.py` for the L0 audit.
'''
(REPO/'docs/RELEASE_STATUS.md').write_text(status, encoding='utf-8')

# Makefile
makefile = '''.PHONY: audit claims package

audit:
	python repro_scripts/audit_release.py

claims:
	python repro_scripts/claim_summary.py

package:
	python repro_scripts/package_repo.py
'''
(REPO/'Makefile').write_text(makefile, encoding='utf-8')

# .gitignore
(REPO/'.gitignore').write_text('''__pycache__/
*.pyc
.DS_Store
.venv/
venv/
.env
''', encoding='utf-8')

# Generate a top-level SHA sums for repo files.
sha_rows=[]
for p in sorted(REPO.rglob('*')):
    if p.is_file() and 'REPO_SHA256SUMS.csv' not in p.name:
        rel = p.relative_to(REPO).as_posix()
        sha_rows.append((rel, p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()))
with (REPO/'REPO_SHA256SUMS.csv').open('w', newline='', encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['path','bytes','sha256']); w.writerows(sha_rows)

# Run audit scripts.
import subprocess
subprocess.run(['python', str(REPO/'repro_scripts/audit_release.py')], check=True)
subprocess.run(['python', str(REPO/'repro_scripts/claim_summary.py')], check=True)

# Zip repo.
zip_path = BASE / f'{REPO.name}.zip'
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
    for p in sorted(REPO.rglob('*')):
        if p.is_file():
            zf.write(p, p.relative_to(REPO.parent))

sha = hashlib.sha256(zip_path.read_bytes()).hexdigest()
print('Created', zip_path)
print('SHA256', sha)
