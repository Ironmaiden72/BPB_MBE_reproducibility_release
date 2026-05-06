# BPB/MBE Framework — Reproducibility Repository

Companion repository to the BPB/MBE paper series submitted to MNRAS:

- **Paper I** — Framework, KiDS-1000 calibration, SPARC galactic test
- **Paper II** — Statistical validation of the galactic branch
- **Paper III** — Biface scalar field completion (LQC)
- **Paper IV** — Five-phase observational confrontation (BAO, CC, SNe, fσ8, KiDS, ACT, DESI, Euclid)
- **Paper V** — KiDS-only posterior (κbf fixed)

Every figure in Papers I, II, IV, and V can be reproduced from this
repository using the commands below. All intermediate data products
(CSV, NPZ, H5) are either shipped as frozen release artefacts or can
be regenerated from the provided scripts.

---

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash scripts/run_smoke_tests.sh      # 12/12 PASS expected
```

---

## Reproduce all figures — one command per paper

```bash
bash scripts/run_reproduce_paper1.sh   # Paper I  — 5 figures
bash scripts/run_reproduce_paper2.sh   # Paper II — 5 figures (+ CSV rebuild)
bash scripts/run_reproduce_paper4.sh   # Paper IV — bridge fig + MCMC figs*
bash scripts/run_reproduce_paperV.sh   # Paper V  — corner plot
```

\* Paper IV MCMC figures require `emcee h5py corner` and the frozen chain
`results/chains/bpb_allphases_chain.h5`. The galaxy-bridge figure
(`figIV_galaxy_bridge.pdf`) runs without the chain.

---

## Regenerate all derived CSV artefacts from raw inputs

```bash
bash scripts/run_build_chains.sh
```

This runs the full step pipeline in dependency order:

| Step | Script | Output | Time |
|------|--------|--------|------|
| 1 | `step1_compute_mouter.py` | `mouter_per_galaxy.csv` | ~1 min |
| 2 | `step2_sign_rule_loo.py` | `loo_predictions.csv` | ~1 min |
| 3 | `step3_null_test.py` | `null_distribution.csv` | ~10 min |
| 4 | `step4_delta_chi2.py` | `delta_chi2_per_galaxy.csv` | ~2 min |
| 4b | `step4b_type95_baseline.py` | `type95_per_galaxy.csv` | ~2 min |
| 5 | `step6_regime2_c200_scan.py` | `regime2_c200_scan.csv` | ~2 min |
| 6 | `step7_build_sparc_active_sample.py` | `sparc_active_sample.csv` | <1 min |
| 7 | `step8_build_chi2_summary.py` | `chi2_summary.csv` | <1 min |
| 8 | `step9_build_regime2_complete.py` | `regime2_complete.csv` | <1 min |

**Total**: ~20 min on a standard laptop.

---

## Frozen release artefacts

Some artefacts are shipped pre-computed because regeneration is expensive
or requires raw survey data not redistributable here:

| Artefact | Role | Regeneration |
|----------|------|--------------|
| `results/chains/bpb_allphases_chain.h5` | Paper IV MCMC chain (92 MB) | `pipeline/bpb_mcmc_v2.py` (~5 h) |
| `results/chains/paperV/kids_kbf_fixed_dynesty.npz` | Paper V posterior | `pipeline/bpb_kids_only_fixed_kbf_v4.py` (~2 h) |
| `results/chains/regime2_complete.csv` | Paper II Table 4 values | `step9_build_regime2_complete.py` (frozen table) |
| `results/chains/kids1000_xi_data.csv` | KiDS-1000 observations | External (Asgari et al. 2021) |
| `results/chains/bpb_power_spectrum.csv` | BPB growth/power grid | `pipeline/bpb_mbe_unified_relic.py` |

---

## Repository structure

```
BPB_MBE/
├── pipeline/               # All analysis and figure scripts
│   ├── step1_*.py  ..  step9_*.py   # CSV generation pipeline
│   ├── figI[1-5]_*.py      # Paper I figures
│   ├── fig[1-5]_*.py       # Paper II figures
│   ├── figIV_*.py          # Paper IV galaxy-bridge figure
│   ├── figPaperIV_all.py   # Paper IV MCMC figures (requires emcee)
│   ├── plot_corner_paperV.py
│   └── bpb_*.py            # Model and MCMC scripts
├── results/
│   ├── chains/             # CSV, H5, NPZ artefacts
│   └── figures/            # All output PDFs
├── scripts/
│   ├── run_smoke_tests.sh
│   ├── run_build_chains.sh
│   └── run_reproduce_paper[1|2|4|V].sh
├── MANIFEST.md             # Full artefact provenance table
├── requirements.txt
└── README.md
```

---

## Dependency graph — figure to data to script

```
KiDS-1000 data ──────────────────────────────► figI1, figI2, figI3

SPARC MRT + NFW reference fits
  ├─ step1 ──► mouter_per_galaxy.csv
  ├─ step2 ──► loo_predictions.csv
  ├─ step3 ──► null_distribution.csv
  ├─ step4 ──► delta_chi2_per_galaxy.csv
  ├─ step6 ──► regime2_c200_scan.csv
  ├─ step7 ──► sparc_active_sample.csv ──► figI4, fig1, fig2, fig3, figIV_bridge
  ├─ step9 ──► regime2_complete.csv   ──► fig5, figI4, figIV_bridge
  └─ step8 ──► chi2_summary.csv       ──► figI4, figIV_bridge

bpb_allphases_chain.h5 (frozen) ────────────► figPaperIV_all
kids_kbf_fixed_dynesty.npz (frozen) ────────► plot_corner_paperV
```

---

## AI disclosure

The preparation of this repository and the associated manuscripts involved
assistance from a large language model (Claude, Anthropic) for code review,
LaTeX editing, and manuscript drafting. All scientific results, analyses,
and conclusions are the sole responsibility of the author.

---

## Data sources

- **SPARC**: Lelli et al. (2016), AJ 152, 157 — http://astroweb.case.edu/SPARC/
- **NFW fits**: Li et al. (2020), ApJS 247, 31 — VizieR J/ApJS/247/31
- **KiDS-1000**: Asgari et al. (2021), A&A 645, A104
- **DESI DR1**: DESI Collaboration (2024), arXiv:2411.12021
- **ACT DR6**: Madhavacheril et al. (2024)
- **Euclid Year 1**: Euclid Collaboration (2024), arXiv:2405.13491

---

## Citation

If you use this repository, please cite Papers I–IV of the BPB/MBE series
(Pelletier 2026a–d, MNRAS submitted) and the data sources listed above.

Data and scripts archived on Zenodo: doi:10.5281/zenodo.19323855 (upon acceptance)
