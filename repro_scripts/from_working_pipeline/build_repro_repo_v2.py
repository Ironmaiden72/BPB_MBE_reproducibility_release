#!/usr/bin/env python3
from __future__ import annotations
import os, shutil, zipfile, hashlib, csv
from pathlib import Path

BASE = Path('/mnt/data')
SRC = BASE/'BPB_MBE_reproducibility_release_v2026_05_06'
OUT = BASE/'BPB_MBE_reproducibility_release_v2026_05_06_v2_code_data'
ZIP_OUT = BASE/'BPB_MBE_reproducibility_release_v2026_05_06_v2_code_data.zip'

if OUT.exists():
    shutil.rmtree(OUT)
shutil.copytree(SRC, OUT)

# Add clean reproduction material while keeping it isolated from the release manuscript layer.
(OUT/'data_locked').mkdir(exist_ok=True)
(OUT/'working_pipeline_snapshot').mkdir(exist_ok=True)
(OUT/'repro_scripts'/'from_working_pipeline').mkdir(parents=True, exist_ok=True)
(OUT/'docs').mkdir(exist_ok=True)

# 1) Unpack compact results and code archives into explicit folders.
def safe_extract_zip(zip_path: Path, dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        for info in z.infolist():
            name = info.filename
            if name.startswith('/') or '..' in Path(name).parts:
                raise RuntimeError(f'Unsafe zip member: {name}')
        z.extractall(dest)

archives_to_unpack = [
    ('results.zip', OUT/'data_locked'/'results_archive'),
    ('pipeline.zip', OUT/'working_pipeline_snapshot'/'pipeline'),
    ('scripts.zip', OUT/'working_pipeline_snapshot'/'scripts'),
    ('docs.zip', OUT/'working_pipeline_snapshot'/'docs'),
]
for zname, dest in archives_to_unpack:
    zpath = BASE/zname
    if zpath.exists():
        safe_extract_zip(zpath, dest)

# 2) Keep original raw archives too for checksum and provenance.
raw_archive_dir = OUT/'data_locked'/'raw_archives'
raw_archive_dir.mkdir(parents=True, exist_ok=True)
for name in [
    'results.zip', 'pipeline.zip', 'scripts.zip', 'papers.zip', 'docs.zip',
    'eos_APR.zip','eos_DD2.zip','eos_DDME2.zip','eos_NL3wrL55.zip','eos_SFHo.zip','eos_SLY4.zip','eos_TW.zip','bsk24_data (2).zip',
    'D2_closure_complete.zip','D4_closed.zip','D4_final.zip','D4_final_v2.zip','D4_sourcelock.zip',
    'cmb_paperIII.zip',
    'paperVI_EG_T1_rWdelta_sqrt_abg_outputs_verified.zip',
    'paperIII_cmb_p3a2_outputs_verified.zip','paperIII_cmb_p3a3_outputs_verified.zip','paperIII_cmb_p3a4_outputs_verified.zip','paperIII_cmb_p3a4b_final_audit_outputs_verified.zip',
    'paperIII_cmb_p3a2_generated_from_uploaded_zip.zip',
]:
    p = BASE/name
    if p.exists():
        shutil.copy2(p, raw_archive_dir/name)

# 3) Copy top-level reusable/audit scripts into a preserved snapshot.
for py in sorted(BASE.glob('*.py')):
    # exclude builders that create this packaging only? keep them too for provenance except image/PDF not an issue
    shutil.copy2(py, OUT/'repro_scripts'/'from_working_pipeline'/py.name)

# 4) Copy project metadata from uploaded clean repo root when available.
for name in ['README.md','MANIFEST.md','CHANGELOG.md','requirements.txt','LICENSE','D5_D6_final_lock.md','session_summary_20260504.md']:
    p = BASE/name
    if p.exists():
        target = OUT/'working_pipeline_snapshot'/('root_'+name if name in ['README.md','requirements.txt','LICENSE'] else name)
        shutil.copy2(p, target)

# 5) Add a v2 scope note.
scope = r'''# Reproducibility release v2 — code/data capsule

This repository is still a clean reproducibility capsule, not the private working laboratory.

The first capsule was intentionally L0/release-index oriented: PDFs, LaTeX sources, corpus indices,
claim locks and checksum audit. That is useful for HAL deposit validation, but it is too thin to feel
like a reproduction repository.

This v2 capsule therefore adds a preserved, isolated reproduction layer:

- `data_locked/results_archive/` — unpacked locked `results.zip` snapshot, including SPARC tables,
  chain-derived CSVs and figure artefacts used by the audited papers where present.
- `working_pipeline_snapshot/pipeline/` — unpacked pipeline scripts snapshot from the uploaded clean bundle.
- `working_pipeline_snapshot/scripts/` — helper shell reproduction scripts snapshot.
- `repro_scripts/from_working_pipeline/` — top-level gate/audit scripts used during the corpus lock.
- `data_locked/raw_archives/` — original uploaded zip archives copied unchanged for checksum/provenance.

Guardrail: this is not a promise that every historical exploratory run can be rerun end-to-end from this
capsule. The intended public contract is:

- L0: verify release checksums and claim/source index.
- L1: inspect and regenerate manuscript tables/figures from locked CSV/JSON where scripts are present.
- L2: rerun light diagnostic gates included in the snapshot.
- L3: heavy MCMC/CAMB/ACT/NS/full pipeline reruns are documented as future/full-pipeline reproduction,
  not required for the HAL corpus timestamp.

In short: v2 contains manuscripts + claim locks + locked data + selected code snapshots, while still excluding
private exploratory history, abandoned branches and non-release clutter.
'''
(OUT/'docs'/'REPRO_SCOPE_V2_CODE_DATA.md').write_text(scope, encoding='utf-8')

# 6) Extend README with v2 note.
readme = OUT/'README.md'
text = readme.read_text(encoding='utf-8') if readme.exists() else ''
add = '''\n\n## v2 code/data reproduction layer\n\nThis v2 capsule includes the original release-index layer plus isolated code/data snapshots:\n\n- `data_locked/results_archive/`\n- `data_locked/raw_archives/`\n- `working_pipeline_snapshot/pipeline/`\n- `working_pipeline_snapshot/scripts/`\n- `repro_scripts/from_working_pipeline/`\n\nSee `docs/REPRO_SCOPE_V2_CODE_DATA.md` for the exact reproduction contract and guardrails.\n'''
if '## v2 code/data reproduction layer' not in text:
    readme.write_text(text.rstrip()+add+'\n', encoding='utf-8')

# 7) Build checksum manifest for v2.
rows=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file():
        rel = p.relative_to(OUT).as_posix()
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        rows.append({'path': rel, 'bytes': p.stat().st_size, 'sha256': h})
with (OUT/'REPO_SHA256SUMS_V2.csv').open('w', newline='', encoding='utf-8') as f:
    w=csv.DictWriter(f, fieldnames=['path','bytes','sha256'])
    w.writeheader(); w.writerows(rows)

# 8) Zip folder.
if ZIP_OUT.exists():
    ZIP_OUT.unlink()
with zipfile.ZipFile(ZIP_OUT, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in sorted(OUT.rglob('*')):
        if p.is_file():
            z.write(p, OUT.name + '/' + p.relative_to(OUT).as_posix())

print('OUT', OUT)
print('ZIP', ZIP_OUT)
print('ZIP_SIZE', ZIP_OUT.stat().st_size)
print('SHA256', hashlib.sha256(ZIP_OUT.read_bytes()).hexdigest())
print('FILES', len(rows))
