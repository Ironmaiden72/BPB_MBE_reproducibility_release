# D5/D6 — Fermeture finale : amplitude SPARC + falsification kernel

## Statut

| Dérivation | Status |
|-----------|--------|
| D5 | `D5_AMPLITUDE_ROWDEPENDENT_KERNEL_RAW_RC_VALIDATED_LOCKED_FOR_PAPERII` |
| D6 | `D6_KERNEL_FALSIFICATION_LOCK_STRONG_PASS` |
| Closure | `PASS_FINAL_D6_LOCK` |

## Nombres verrouillés (Paper II)

### D5 — Kernel amplitude raw RC

| Quantité | Valeur |
|---------|--------|
| Primary raw recovered gap fraction | 0.773245496307141 |
| Primary raw total Δχ² | −4302.455767865691 |
| Primary raw delta vs archive | 302.77382413779014 |
| c200-only raw recovered gap fraction | 0.777719560160886 |

### D6 — Falsification par permutation

| Quantité | Valeur |
|---------|--------|
| N_PERM | 5000 |
| p-value | 0.0001999600079984003 = 1/(5000+1) |
| Interprétation | 0 permutations sur 5000 battent le kernel réel |

### Robustesse top-removal

| Retrait | Fraction récupérée |
|---------|-------------------|
| top-1 retiré | 0.832030 |
| top-3 retirés | 0.746563 |
| top-5 retirés | 0.580735 |
| top-10 retirés | 0.796819 |

## Wording Paper II (section D5/D6)

> D5 amplitude closure and raw rotation-curve validation.
>
> The global Weyl-matter normalization is not sufficient as a full SPARC
> amplitude closure. The missing factor is row-dependent, with a median
> archive-over-clean response scale of 0.9336 and a top-five residual
> localization fraction of 0.8079.
>
> The successful closure is a source-only row-dependent kernel. The primary
> leave-one-out model using halo mass and concentration source-side predictors
> recovers 0.778483 of the clean-to-archive gap (Δχ² = −4309.45). Repeated
> five-fold cross-validation gives a median recovered fraction of 0.7904.
>
> D5-A7 injects the fixed kernel into the raw point-level SPARC rotation-curve
> calculation. The primary raw model gives Δχ² = −4302.456, recovered gap
> fraction = 0.773245, with no invalid rows.
>
> D6 falsification: in 5000 source-feature permutations, no permuted model
> reaches the real recovered fraction of 0.773245, giving the empirical
> resolution-limited p-value p = 1/(5000+1) = 0.00019996.
>
> **Guardrail:** this result establishes that the row-dependent kernel
> recovers *most* of the clean-to-archive amplitude gap. Exact archive
> (Δχ² = −4605.23) or oracle (Δχ² = −4613.56) recovery is not claimed.

## Catalogue D — état complet post D5/D6

| ID | Statut | Verrou principal |
|----|--------|-----------------|
| D1 | ✅ FERMÉ | a₀ = cH₀(1−s₀) |
| D2 | ✅ FERMÉ | m_FP = 2.166 H₀, 16/16 |
| D3 | ⚠️ 3% | y_t²(1−a_bg)=a_bg³ + leakage HR→E_G |
| D4 | ✅ FERMÉ | α²_NS = 3/(32π²), ΔΛ/Λ = −2.72% |
| D5 | ✅ FERMÉ | kernel row-dep raw RC validé |
| D6 | ✅ FERMÉ | p = 1/5001, 0/5000 permutations |
| D7 | ⚠️ mécanisme | Δn_s = +0.010 via face-minus, CAMB requis |
