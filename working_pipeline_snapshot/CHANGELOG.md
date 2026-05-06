# Changelog

## v1.0.0 — 21 March 2026

### Initial release

**Pipeline:**
- `bpb_mbe_unified_relic.py` — core BPB/MBE physics (background + growth)
- `bpb_mcmc_v2.py` — main MCMC integrating Phases 1–5
- `bpb_phase2_kids.py` — KiDS-1000 S8 likelihood
- `bpb_phase3_act.py` — ACT DR6 + Planck CMB lensing (sigma8 constraint)
- `bpb_phase4_desi.py` — DESI DR1 f*sigma8 at 7 redshift bins
- `bpb_phase5_euclid.py` — Euclid Year 1 S8 (published) + fsig8 (forecast)

**Key results (all-phases MCMC, 5000 steps, 32 walkers):**
- beta_c = -0.15 +0.03/-0.04  (4.3σ, beta_c < 0)
- sigma8 = 0.82 ± 0.01
- H0 = 70.4 ± 1.6 km/s/Mpc
- kappa_bf = 2.61 +1.50/-1.52
- a_bg = 0.26 +0.13/-0.11  (z_t ~ 2.8)

**Papers:**
- Paper III (scalar field completion) — v.E, review-ready
- Papers I, II — submitted simultaneously to MNRAS
- Paper IV — in preparation

## v1.1.0 — 29 March 2026

### Paper V added

**New scripts:**
- `pipeline/bpb_kids_only_fixed_kbf_v4.py` — KiDS-1000 full-likelihood emulator
  (PCA+MLP, 10k evaluations, dynesty nested sampling, κ_bf fixed at Paper IV MAP)
- `pipeline/plot_corner_paperV.py` — Corner plot (β_c, σ8, ω_m, S8_eff)

**New frozen results:**
- `results/chains/paperV/kids_kbf_fixed_dynesty.npz` — posterior samples
  (5709 weighted samples, ln Z = -196.65 ± 0.14, RMSE = 3.25 lnL)
- `results/figures/figV1_beta_c_posterior.pdf` — marginalised β_c posterior
- `results/figures/figV2_corner.pdf` — corner (β_c, σ8, ω_m, S8_eff)

**Key result:**
- β_c = -0.428 +0.099/-0.052, P(β_c < 0) = 99.9%  (KiDS-1000 only, κ_bf fixed)

**Reproduction:**
```bash
bash scripts/run_reproduce_paperV.sh
```
