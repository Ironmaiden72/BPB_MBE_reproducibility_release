# Paper I — Figure Scripts

## Structure
```
paperI/
  figI1_kids_xi.py          Fig 1: KiDS-1000 xi+/xi- (6 panels)
  figI2_chi2_bridge.py      Fig 2: chi2(gamma_p) landscape + KiDS->SPARC bridge
  figI3_power_growth.py     Fig 3: C_l^BPB/C_l^ref + D_BPB(z)/D_ref(z)
  figI4_galactic_summary.py Fig 4: Galactic branch summary (4 panels)
  figI5_chi2_landscapes.py  Fig 5: chi2(A) landscapes IC2574 + DDO170
  gen_data_paperI.py        Regenerate data CSVs from scratch
  data/
    kids1000_xi_data.csv    KiDS-1000 xi+/xi- data + BPB + LCDM models
    kids1000_chi2_scan.csv  chi2 vs gamma_p scan
    bpb_power_spectrum.csv  C_l ratio + D(z) growth ratio
    ../data/                Shared with Paper II:
      sparc_active_sample.csv
      regime2_complete.csv
```

## Usage
  python gen_data_paperI.py     # generate CSVs (already done)
  python figI1_kids_xi.py       # -> figI1_kids_xi.pdf
  ...etc

## Key hardcoded paper values (annotations)
  chi2_BPB=400.28, chi2_LCDM=401.30, Dchi2=-1.02
  gp_best=0.15, H0=68, Om=0.19, Om_eff=0.352
  C_l enhancement: 2-10% (diagonal bins)
  D(z) growth: ~5% at z=0.5

## Note on Figs 4 and 5
  Figs 4 and 5 of Paper I are summary versions of Paper II Figs 3 and 4.
  They use shared data from ../data/ (sparc_active_sample.csv etc.).
  Replace with real SPARC pipeline output for exact reproduction.
