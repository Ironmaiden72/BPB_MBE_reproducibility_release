
## Provenance complète des artefacts

### Chaîne de génération (Papers I & II)

| Artefact | Script générateur | Statut | Temps |
|---|---|---|---|
| `mouter_per_galaxy.csv` | `step1_compute_mouter.py` | rebuildable | ~1 min |
| `loo_predictions.csv` | `step2_sign_rule_loo.py` | rebuildable | ~1 min |
| `null_distribution.csv` | `step3_null_test.py` | rebuildable | ~10 min |
| `delta_chi2_per_galaxy.csv` | `step4_delta_chi2.py` | rebuildable | ~2 min |
| `regime2_c200_scan.csv` | `step6_regime2_c200_scan.py` | rebuildable | ~2 min |
| `sparc_active_sample.csv` | `step7_build_sparc_active_sample.py` | rebuildable | <1 min |
| `chi2_summary.csv` | `step8_build_chi2_summary.py` | rebuildable | <1 min |
| `regime2_complete.csv` | `step9_build_regime2_complete.py` | rebuildable | <1 min |

### Artefacts figés (frozen release artifacts)

| Artefact | Rôle | Script producteur | Statut |
|---|---|---|---|
| `bpb_allphases_chain.h5` | Paper IV MCMC chain | `bpb_mcmc_v2.py` | frozen (~5h rerun) |
| `kids_kbf_fixed_dynesty.npz` | Paper V posterior | `bpb_kids_only_fixed_kbf_v4.py` | rebuildable (~2h) |
| `kids1000_xi_data.csv` | KiDS-1000 observations | external data | frozen |
| `bpb_power_spectrum.csv` | BPB growth/power | `bpb_mbe_unified_relic.py` | frozen |

### Pour tout régénérer

```bash
# Régénérer tous les CSV (hors chain Paper IV, ~20 min)
bash scripts/run_build_chains.sh

# Régénérer les figures Paper I
bash scripts/run_reproduce_paper1.sh

# Régénérer les figures Paper II
bash scripts/run_reproduce_paper2.sh

# Régénérer les figures Paper V (depuis NPZ gelé)
bash scripts/run_reproduce_paperV.sh
```
