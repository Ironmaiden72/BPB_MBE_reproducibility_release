#!/usr/bin/env python3
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
