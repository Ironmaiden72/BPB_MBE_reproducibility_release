# BPB/MBE reproducibility release capsule

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

## v2 code/data reproduction layer

This v2 capsule includes the original release-index layer plus isolated code/data snapshots:

- `data_locked/results_archive/`
- `data_locked/raw_archives/`
- `working_pipeline_snapshot/pipeline/`
- `working_pipeline_snapshot/scripts/`
- `repro_scripts/from_working_pipeline/`

See `docs/REPRO_SCOPE_V2_CODE_DATA.md` for the exact reproduction contract and guardrails.

## Citation and archive

Canonical archive DOI:

- Zenodo DOI: [10.5281/zenodo.20048680](https://doi.org/10.5281/zenodo.20048680)
- GitHub repository: https://github.com/Ironmaiden72/BPB_MBE_reproducibility_release
- GitHub release: https://github.com/Ironmaiden72/BPB_MBE_reproducibility_release/releases/tag/v2026.05.06-hal-corpus
- Author ORCID: [0009-0009-0862-2377](https://orcid.org/0009-0009-0862-2377)

Suggested citation:

> Pelletier, F. (2026). *BPB/MBE HAL Corpus Reproducibility Release v2026.05.06*. Zenodo. https://doi.org/10.5281/zenodo.20048680
