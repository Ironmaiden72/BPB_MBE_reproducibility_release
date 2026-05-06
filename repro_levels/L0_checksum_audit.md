# L0 checksum and release-index audit

Run:

```bash
python repro_scripts/audit_release.py
python repro_scripts/claim_summary.py
```

This verifies the final release payload against `release/release_index/00_PACKAGE_SHA256SUMS.csv`, checks that the manifest contains the expected paper versions, and summarizes the claim-source index.
