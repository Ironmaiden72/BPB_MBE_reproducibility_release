#!/usr/bin/env python3
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
