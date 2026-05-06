.PHONY: audit claims package

audit:
	python repro_scripts/audit_release.py

claims:
	python repro_scripts/claim_summary.py

package:
	python repro_scripts/package_repo.py
