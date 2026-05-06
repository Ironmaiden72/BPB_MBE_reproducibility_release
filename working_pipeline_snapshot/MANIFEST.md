# Reproducibility manifest

## Paper I

| Figure | Script | Inputs | Output |
|---|---|---|---|
| Fig I1 | `pipeline/figI1_kids_xi.py` | `results/chains/kids1000_xi_data.csv` | `results/figures/figI1_kids_xi.pdf` |
| Fig I2 | `pipeline/figI2_chi2_bridge.py` | `results/chains/kids1000_chi2_scan.csv` | `results/figures/figI2_chi2_bridge.pdf` |
| Fig I3 | `pipeline/figI3_power_growth.py` | `results/chains/bpb_power_spectrum.csv` | `results/figures/figI3_power_growth.pdf` |
| Fig I4 | `pipeline/figI4_galactic_summary.py` | `results/chains/sparc_active_sample.csv`, `results/chains/regime2_complete.csv` | `results/figures/figI4_galactic_summary.pdf` |
| Fig I5 | `pipeline/figI5_chi2_landscapes.py` | analytic display figure | `results/figures/figI5_chi2_landscapes.pdf` |

## Paper II

| Figure | Script | Inputs | Output |
|---|---|---|---|
| Fig 1a/1b | `pipeline/fig1_bpb_signal.py` | `results/chains/sparc_active_sample.csv` | `results/figures/fig1a_signal_hist.pdf`, `results/figures/fig1b_decision_space.pdf` |
| Fig 2a/2b | `pipeline/fig2_areq_scatter.py` | `results/chains/sparc_active_sample.csv` | `results/figures/fig2a_areq_physical.pdf`, `results/figures/fig2b_areq_oracle.pdf` |
| Fig 3a-3d | `pipeline/fig3_validation.py` | `results/chains/sparc_active_sample.csv`, `results/chains/regime2_complete.csv` | `results/figures/fig3a_null.pdf`, `fig3b_pergal.pdf`, `fig3c_regimes.pdf`, `fig3d_chi2summary.pdf` |
| Fig 4a/4b | `pipeline/fig4_chi2_landscapes.py` | analytic display figure | `results/figures/fig4a_IC2574.pdf`, `results/figures/fig4b_DDO170.pdf` |
| Fig 5 | `pipeline/fig5_regime2.py` | embedded table values | `results/figures/fig5_regime2.pdf` |

## Paper IV

| Figure / result | Script / file | Inputs | Output |
|---|---|---|---|
| All-phases chain | `results/chains/bpb_allphases_chain.h5` | frozen chain | frozen chain |
| Existing publication figures | shipped under `results/figures/` | frozen chain and shipped plotting scripts | existing PDF outputs |

## Paper V

| Figure | Script | Inputs | Output |
|---|---|---|---|
| Fig V1 (β_c posterior) | `pipeline/bpb_kids_only_fixed_kbf_v4.py --plot_only` | `results/chains/paperV/kids_kbf_fixed_dynesty.npz` | `results/figures/figV1_beta_c_posterior.pdf` |
| Fig V2 (corner) | `pipeline/plot_corner_paperV.py` | `results/chains/paperV/kids_kbf_fixed_dynesty.npz` | `results/figures/figV2_corner.pdf` |

### Frozen posterior

| File | Description | Samples | ln Z |
|---|---|---|---|
| `results/chains/paperV/kids_kbf_fixed_dynesty.npz` | dynesty posterior, κ_bf=2.706 fixed | 5709 | -196.65 ± 0.14 |

**Keys in NPZ:** `samples` (5709×14), `weights`, `logz`, `logzerr`, `logwt`, `logl`,
`param_names`, `theta_map`, `kappa_bf_fixed`, `converged`.

**Parameter order:** β_c, σ8, ω_m, A_IA, m1–m5, dz1–dz5.

### Full re-run from scratch (requires KiDS bridge file)

```bash
# Step 1 — generate 10k likelihood evaluations (~7 min, 12 cores)
python pipeline/bpb_kids_only_fixed_kbf_v4.py \
    --bridge outputs_shear/kids_real_bridge/kids_xipm_bridge.npz \
    --output_dir results/chains/paperV \
    --n_eval 10000 --n_cores 12

# Step 2 — train emulator + dynesty (~2h, nlive=500)
python pipeline/bpb_kids_only_fixed_kbf_v4.py \
    --bridge outputs_shear/kids_real_bridge/kids_xipm_bridge.npz \
    --output_dir results/chains/paperV \
    --skip_eval --nlive 500
```

| `type95_per_galaxy.csv` | `step4b_type95_baseline.py` | rebuildable | ~2 min |
