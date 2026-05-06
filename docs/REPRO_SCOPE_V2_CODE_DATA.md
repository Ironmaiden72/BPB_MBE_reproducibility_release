# Reproducibility release v2 — code/data capsule

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
