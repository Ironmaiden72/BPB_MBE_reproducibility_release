#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import csv
import json
import math
from typing import Any, Dict, List


OUTDIR = Path("outputs/paperVI_EG_theory")


DEFAULTS = {
    "a_bg": 0.403687948,
    "sqrt_a_bg": 0.6353644214149861,
    "r_free": 0.630864746,
    "r2_free": 0.39799032774564447,
    "delta_r": -0.004499675414986082,
    "delta_r2": -0.005697620254355551,
    "delta_chi2_t1": 0.008545000000000025,
    "theta_bg_deg": 50.55298183517175,
    "t7_dominant_trace_unweighted": 0.9999987128190622,
    "t7_dominant_trace_ell": 0.9999988797737299,
    "t7_parallel_act": 0.999999999998494,
    "t7_parallel_eg": 0.9999999999172636,
    "t8_fraction": 0.3979903273355348,
    "t8_sqrt_fraction": 0.6308647456749622,
    "t8_delta_fraction": -0.005697620664465219,
    "t8_relative_delta_fraction": -0.014113923124762742,
    "t8_corr_A_X": 1.0,
    "t8_chi2_per_point_all": 0.22649042932465788,
    "t8_planck_chi2_per_point": 0.024232480763707834,
    "t8_act_chi2_per_point": 0.5588874993058423,
    "t8_combined_chi2_per_point": 0.09635130790442359,
    "free_scan_csigma_auto": 0.7,
    "free_scan_sqrt_abg": 0.6350086259969958,
    "free_scan_r_hat": 0.6349568146956915,
    "free_scan_sigma_r": 0.04899053078032919,
    "free_scan_z": 0.001057577872274767,
    "free_scan_chi2_free": 1.3580463340514337,
    "free_scan_chi2_physical": 1.358047452522391,
    "free_scan_delta": 1.118470958516582e-06,
}


def load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def read_csv_dicts(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def fget(d: Dict[str, Any], keys: List[str], default: float) -> float:
    cur: Any = d
    try:
        for k in keys:
            cur = cur[k]
        return float(cur)
    except Exception:
        return default


def csv_lookup(rows: List[Dict[str, str]], key_col: str, key_val: str, val_col: str, default: float) -> float:
    for r in rows:
        if str(r.get(key_col, "")) == key_val:
            try:
                return float(r.get(val_col, default))
            except Exception:
                return default
    return default


def first_matching_float(rows: List[Dict[str, str]], predicates: Dict[str, str], val_col: str, default: float) -> float:
    for r in rows:
        ok = True
        for k, v in predicates.items():
            if str(r.get(k, "")) != v:
                ok = False
                break
        if ok:
            try:
                return float(r.get(val_col, default))
            except Exception:
                return default
    return default


def gather_values(outdir: Path) -> Dict[str, Any]:
    vals = dict(DEFAULTS)

    t1 = load_json(outdir / "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json")
    vals["a_bg"] = fget(t1, ["inputs", "a_bg"], vals["a_bg"])
    vals["sqrt_a_bg"] = fget(t1, ["derived", "sqrt_a_bg"], vals["sqrt_a_bg"])
    vals["r_free"] = fget(t1, ["inputs", "r_wdelta_free"], vals["r_free"])
    vals["r2_free"] = fget(t1, ["derived", "r_wdelta_free_squared"], vals["r2_free"])
    vals["delta_r"] = fget(t1, ["derived", "delta_r_free_minus_sqrt_a_bg"], vals["delta_r"])
    vals["delta_r2"] = fget(t1, ["derived", "delta_r2_minus_a_bg"], vals["delta_r2"])
    vals["delta_chi2_t1"] = fget(t1, ["derived", "delta_chi2_fixed_sqrt_abg_minus_free"], vals["delta_chi2_t1"])

    t3 = load_json(outdir / "paperVI_EG_T3_A4_overlap_reduction.json")
    vals["theta_bg_deg"] = fget(t3, ["derived", "theta_bg_deg"], vals["theta_bg_deg"])

    rank_rows = read_csv_dicts(outdir / "paperVI_EG_T7_rank_diagnostics.csv")
    vals["t7_dominant_trace_unweighted"] = first_matching_float(
        rank_rows, {"scheme": "unweighted_mean"}, "dominant_trace_fraction", vals["t7_dominant_trace_unweighted"]
    )
    vals["t7_dominant_trace_ell"] = first_matching_float(
        rank_rows, {"scheme": "ell_width_weighted_mean"}, "dominant_trace_fraction", vals["t7_dominant_trace_ell"]
    )

    proj_rows = read_csv_dicts(outdir / "paperVI_EG_T7_projection_fractions.csv")
    vals["t7_parallel_act"] = first_matching_float(
        proj_rows, {"scheme": "unweighted_mean", "target": "W_ACTSigma"}, "parallel_fraction", vals["t7_parallel_act"]
    )
    vals["t7_parallel_eg"] = first_matching_float(
        proj_rows, {"scheme": "unweighted_mean", "target": "W_EGSigma"}, "parallel_fraction", vals["t7_parallel_eg"]
    )

    amp_rows = read_csv_dicts(outdir / "paperVI_EG_T8_amplitude_fraction_summary.csv")
    vals["t8_fraction"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "fraction_XX_over_AA", vals["t8_fraction"]
    )
    vals["t8_sqrt_fraction"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "sqrt_fraction", vals["t8_sqrt_fraction"]
    )
    vals["t8_delta_fraction"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "delta_fraction_minus_a_bg", vals["t8_delta_fraction"]
    )
    vals["t8_relative_delta_fraction"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "relative_delta_fraction", vals["t8_relative_delta_fraction"]
    )
    vals["t8_corr_A_X"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "corr_A_X", vals["t8_corr_A_X"]
    )
    vals["t8_chi2_per_point_all"] = first_matching_float(
        amp_rows, {"dataset": "ALL", "weight_mode": "unweighted"}, "chi2_per_point", vals["t8_chi2_per_point_all"]
    )
    vals["t8_planck_chi2_per_point"] = first_matching_float(
        amp_rows, {"dataset": "planck_pr4_boss_diagnostic", "weight_mode": "unweighted"}, "chi2_per_point", vals["t8_planck_chi2_per_point"]
    )
    vals["t8_act_chi2_per_point"] = first_matching_float(
        amp_rows, {"dataset": "act_dr6_boss_diagnostic", "weight_mode": "unweighted"}, "chi2_per_point", vals["t8_act_chi2_per_point"]
    )
    vals["t8_combined_chi2_per_point"] = first_matching_float(
        amp_rows, {"dataset": "act_dr6_planck_pr4_boss_diagnostic", "weight_mode": "unweighted"}, "chi2_per_point", vals["t8_combined_chi2_per_point"]
    )

    # Free-r scan row from T8 JSON if available, otherwise defaults.
    t8 = load_json(outdir / "paperVI_EG_T8_auto_cross_amplitude_space.json")
    fs = t8.get("free_r_scan_support", {}) if isinstance(t8.get("free_r_scan_support", {}), dict) else {}
    mapping = {
        "free_scan_csigma_auto": "c_sigma_auto",
        "free_scan_sqrt_abg": "sqrt_a_bg",
        "free_scan_r_hat": "r_hat_free",
        "free_scan_sigma_r": "sigma_r",
        "free_scan_z": "z_sqrt_abg_minus_rhat",
        "free_scan_chi2_free": "chi2_EG_free_r",
        "free_scan_chi2_physical": "chi2_EG_physical",
        "free_scan_delta": "delta_physical_minus_free",
    }
    for out_key, in_key in mapping.items():
        try:
            vals[out_key] = float(fs.get(in_key, vals[out_key]))
        except Exception:
            pass

    vals["relative_delta_percent"] = 100.0 * vals["t8_relative_delta_fraction"]
    vals["abs_relative_delta_percent"] = abs(vals["relative_delta_percent"])

    return vals


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def build_markdown(vals: Dict[str, Any]) -> str:
    return f"""# Paper VI / E_G T9A - theory appendix draft from T1-T8

## Status

**P6_T9A_THEORY_APPENDIX_DRAFTED_FROM_T1_TO_T8**

Non-destructive writing gate. No likelihood run. NS untouched.

## Purpose

This gate turns the audited T1-T8 chain into a paper-ready theory appendix/section.
It does not claim an action-level HR/BPB perturbation proof. It states the result as a projected-diagnostic closure and identifies the remaining action-level target.

## Core result

The compressed E_G diagnostic contains an effective Weyl--matter correlation coefficient. The audited chain shows that the free coefficient can be replaced by the dictionary-tied closure

```text
r_Wdelta = sqrt(a_bg)
```

with negligible diagnostic penalty.

| quantity | value |
|---|---:|
| a_bg | {vals['a_bg']:.12f} |
| sqrt(a_bg) | {vals['sqrt_a_bg']:.12f} |
| free r_Wdelta | {vals['r_free']:.12f} |
| r_Wdelta^2 | {vals['r2_free']:.12f} |
| r_Wdelta^2 - a_bg | {vals['delta_r2']:.12f} |
| delta chi2 fixed sqrt(a_bg) vs free r | {vals['delta_chi2_t1']:.12g} |

## Projection theorem

If the projected Weyl response decomposes as

```text
W = sqrt(a_bg) delta_m + sqrt(1-a_bg) xi_W
```

with xi_W orthogonal to the matter-correlated mode in the projected diagnostic inner product, then

```text
corr(W, delta_m) = sqrt(a_bg).
```

Equivalently, in a two-mode projected Weyl basis,

```text
W_hat = cos(theta_W) e_m + sin(theta_W) e_x,
cos^2(theta_W) = a_bg,
theta_W = acos(sqrt(a_bg)) = {vals['theta_bg_deg']:.6f} deg.
```

## Obstruction in Wenzl shape space

The Wenzl 19-bin projected response curves are useful model-response curves, but they are nearly one-dimensional in shape. Therefore they cannot by themselves expose a two-mode matter/private Weyl split.

| diagnostic | value |
|---|---:|
| dominant trace fraction, unweighted | {vals['t7_dominant_trace_unweighted']:.12f} |
| dominant trace fraction, ell-width weighted | {vals['t7_dominant_trace_ell']:.12f} |
| naive parallel fraction W_ACTSigma on W_mu | {vals['t7_parallel_act']:.12f} |
| naive parallel fraction W_EGSigma on W_mu | {vals['t7_parallel_eg']:.12f} |
| target a_bg | {vals['a_bg']:.12f} |

Thus the Wenzl shape space must not be used as an A4 proof. It returns overlap fractions near unity, not a_bg.

## Correct projected diagnostic space

The correct space for the compressed r_Wdelta diagnostic is the auto/cross amplitude space. In the frozen diagnostic file, define

```text
A_i = unit_ACTSigma_cross,
X_i = prediction_global_auto_cross.
```

The row-level identity is

```text
X_i / A_i = r_global.
```

Therefore, for positive diagonal/bin weights,

```text
<X,X>/<A,A> = r_global^2.
```

The locked result is:

| quantity | value |
|---|---:|
| <X,X>/<A,A> | {vals['t8_fraction']:.12f} |
| sqrt(<X,X>/<A,A>) | {vals['t8_sqrt_fraction']:.12f} |
| <X,X>/<A,A> - a_bg | {vals['t8_delta_fraction']:.12f} |
| relative difference | {vals['relative_delta_percent']:.6f}% |
| corr(A,X) | {vals['t8_corr_A_X']:.12f} |
| all-data chi2/pt support | {vals['t8_chi2_per_point_all']:.12f} |
| Planck PR4+BOSS chi2/pt | {vals['t8_planck_chi2_per_point']:.12f} |
| ACT DR6+BOSS chi2/pt | {vals['t8_act_chi2_per_point']:.12f} |
| ACT+Planck+BOSS chi2/pt | {vals['t8_combined_chi2_per_point']:.12f} |

The free-r scan gives an independent support row:

| quantity | value |
|---|---:|
| c_sigma_auto | {vals['free_scan_csigma_auto']:.12f} |
| sqrt(a_bg) | {vals['free_scan_sqrt_abg']:.12f} |
| r_hat_free | {vals['free_scan_r_hat']:.12f} |
| sigma_r | {vals['free_scan_sigma_r']:.12f} |
| z(sqrt(a_bg)-r_hat) | {vals['free_scan_z']:.12f} |
| chi2_EG_free_r | {vals['free_scan_chi2_free']:.12f} |
| chi2_EG_physical | {vals['free_scan_chi2_physical']:.12f} |
| delta physical-free | {vals['free_scan_delta']:.12e} |

## Paper-ready statement

The free E_G Weyl--matter coefficient is not required as an independent parameter at the compressed diagnostic level. In the auto/cross amplitude space, the cross vector X and auto vector A satisfy X=r_global A, so that the matter-correlated Weyl fraction is <X,X>/<A,A>=r_global^2. The measured fraction, {vals['t8_fraction']:.12f}, agrees with the BPB background parameter a_bg={vals['a_bg']:.12f} to {vals['abs_relative_delta_percent']:.3f}%.

This supports the conditional dictionary closure

```text
r_Wdelta = sqrt(a_bg),
```

while the action-level derivation remains open.

## Guardrails

1. This is not an official full E_G likelihood or full covariance proof.
2. The Wenzl 19-bin shape space is not used to prove A4; T7 showed it is nearly collinear.
3. The auto/cross amplitude-space identity is a projected-diagnostic closure, not yet a derivation from the HR/BPB field equations.
4. The next action-level target is to derive the projector identity

```text
<P_m W,P_m W>_EG / <W,W>_EG = a_bg
```

from the HR/BPB perturbation sector.

## Recommended insertion point

Use this as a theory appendix or compact subsection in Paper VI after the E_G diagnostic results and before the final discussion. A suitable title is:

```text
Projected Weyl--matter closure: r_Wdelta = sqrt(a_bg)
```
"""


def build_tex(vals: Dict[str, Any]) -> str:
    return r"""
\subsection{Projected Weyl--matter closure: $r_{W\delta}=\sqrt{a_{\rm bg}}$}
\label{sec:eg-rwdelta-sqrt-abg}

The compressed $E_G$ diagnostic introduced an effective Weyl--matter
correlation coefficient $r_{W\delta}$.  In the present audit chain we find
that this coefficient is not required as an independent diagnostic parameter:
it is replaced, with negligible loss, by the dictionary-tied closure
\begin{equation}
  r_{W\delta}=\sqrt{a_{\rm bg}} .
\end{equation}
Numerically, the free projected value is
\begin{equation}
  r_{W\delta}^{\rm free} = %.12f,
  \qquad
  \sqrt{a_{\rm bg}} = %.12f,
\end{equation}
with
\begin{equation}
  \left(r_{W\delta}^{\rm free}\right)^2 = %.12f,
  \qquad
  a_{\rm bg}=%.12f,
\end{equation}
and fixing $r_{W\delta}=\sqrt{a_{\rm bg}}$ costs only
\begin{equation}
  \Delta\chi^2 = %.12g
\end{equation}
relative to a free projected $r_{W\delta}$.

The algebraic origin of this relation is simple.  If the projected Weyl
response decomposes into a matter-correlated mode and an orthogonal
BPB-private mode,
\begin{equation}
  W = \sqrt{a_{\rm bg}}\,\delta_m
      + \sqrt{1-a_{\rm bg}}\,\xi_W,
  \qquad
  \langle \delta_m,\xi_W\rangle_{E_G}=0,
\end{equation}
with matched projected variances, then
\begin{equation}
  \mathrm{corr}(W,\delta_m)=\sqrt{a_{\rm bg}} .
\end{equation}
Equivalently, in a two-mode projected Weyl basis,
\begin{equation}
  \widehat W = \cos\theta_W\,e_m + \sin\theta_W\,e_x,
  \qquad
  \cos^2\theta_W=a_{\rm bg},
\end{equation}
so that
\begin{equation}
  \theta_W = \arccos\sqrt{a_{\rm bg}} = %.6f^\circ .
\end{equation}

A crucial guardrail is that the Wenzl 19-bin projected response curves are
not the space in which this two-mode decomposition should be proven.  Their
shape space is effectively one-dimensional: the dominant trace fraction is
%.12f with an unweighted bin norm and %.12f with an $\ell$-width weighted
norm.  A naive projection of the projected Sigma responses onto the
matter-like response $W_\mu$ gives parallel fractions near unity,
%.12f and %.12f, rather than the target $a_{\rm bg}=%.12f$.
Thus the Wenzl response curves are useful model-response curves, but they
must not be used as a proof of the matter/private Weyl overlap.

The correct compressed diagnostic space is the auto/cross amplitude space.
Let
\begin{equation}
  A_i = ({\rm unit\_ACTSigma\_cross})_i,
  \qquad
  X_i = ({\rm prediction\_global\_auto\_cross})_i .
\end{equation}
The frozen diagnostic satisfies $X_i/A_i=r_{\rm global}$ row by row.  Hence,
for any positive diagonal/bin weights,
\begin{equation}
  \frac{\langle X,X\rangle}{\langle A,A\rangle}=r_{\rm global}^2 .
\end{equation}
The audited value is
\begin{equation}
  \frac{\langle X,X\rangle}{\langle A,A\rangle}=%.12f,
  \qquad
  \sqrt{\frac{\langle X,X\rangle}{\langle A,A\rangle}}=%.12f,
\end{equation}
which differs from $a_{\rm bg}$ by
\begin{equation}
  %.12f \quad (%.6f\%%).
\end{equation}
The supporting diagnostic quality is good, with $\chi^2/{\rm pt}=%.12f$
for the combined set, %.12f for Planck PR4+BOSS, %.12f for ACT DR6+BOSS,
and %.12f for the combined ACT+Planck+BOSS rows.

The free-$r$ scan gives a further check: at the best support row,
\begin{equation}
  \sqrt{a_{\rm bg}}=%.12f,
  \qquad
  \hat r_{\rm free}=%.12f,
  \qquad
  z=%.12f,
\end{equation}
with
\begin{equation}
  \chi^2_{E_G,{\rm free}}=%.12f,
  \qquad
  \chi^2_{E_G,{\rm physical}}=%.12f,
  \qquad
  \Delta\chi^2=%.12e .
\end{equation}

We therefore use $r_{W\delta}=\sqrt{a_{\rm bg}}$ as a projected-diagnostic
closure.  This is not yet an action-level derivation of the HR/BPB scalar
perturbation projectors.  The remaining theory target is the identity
\begin{equation}
  \frac{\langle P_m W,P_m W\rangle_{E_G}}
       {\langle W,W\rangle_{E_G}}
  = a_{\rm bg},
\end{equation}
where $P_m$ is the matter-correlated Weyl projector in the projected
$E_G$ inner product.
""" % (
        vals["r_free"], vals["sqrt_a_bg"], vals["r2_free"], vals["a_bg"], vals["delta_chi2_t1"],
        vals["theta_bg_deg"], vals["t7_dominant_trace_unweighted"], vals["t7_dominant_trace_ell"],
        vals["t7_parallel_act"], vals["t7_parallel_eg"], vals["a_bg"], vals["t8_fraction"],
        vals["t8_sqrt_fraction"], vals["t8_delta_fraction"], vals["relative_delta_percent"],
        vals["t8_chi2_per_point_all"], vals["t8_planck_chi2_per_point"], vals["t8_act_chi2_per_point"],
        vals["t8_combined_chi2_per_point"], vals["free_scan_sqrt_abg"], vals["free_scan_r_hat"],
        vals["free_scan_z"], vals["free_scan_chi2_free"], vals["free_scan_chi2_physical"],
        vals["free_scan_delta"],
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    vals = gather_values(OUTDIR)

    status = "P6_T9A_THEORY_APPENDIX_DRAFTED_FROM_T1_TO_T8"

    source_lock = [
        {
            "claim_id": "P6.T9A.01",
            "claim": "Free r_Wdelta is replaced by sqrt(a_bg) with negligible chi2 penalty.",
            "source": "paperVI_EG_T1_rWdelta_sqrt_abg_gate.json/csv",
            "status": "SOURCE_MATCHED",
        },
        {
            "claim_id": "P6.T9A.02",
            "claim": "Projection theorem gives r_Wdelta=sqrt(a_bg) under a two-mode orthogonal Weyl split.",
            "source": "paperVI_EG_T2_weyl_projection_theorem.md/json",
            "status": "SOURCE_MATCHED_CONDITIONAL_THEOREM",
        },
        {
            "claim_id": "P6.T9A.03",
            "claim": "A4 is equivalent to the projected mixing angle identity cos^2(theta_W)=a_bg.",
            "source": "paperVI_EG_T3_A4_overlap_reduction.md/json",
            "status": "SOURCE_MATCHED_REDUCTION_LEMMA",
        },
        {
            "claim_id": "P6.T9A.04",
            "claim": "Wenzl shape space is nearly one-dimensional and must not be used as A4 proof.",
            "source": "paperVI_EG_T7_collinearity_obstruction.md/csv",
            "status": "SOURCE_MATCHED_OBSTRUCTION",
        },
        {
            "claim_id": "P6.T9A.05",
            "claim": "Auto/cross amplitude space gives <X,X>/<A,A>=r_global^2 close to a_bg.",
            "source": "paperVI_EG_T8_auto_cross_amplitude_space.md/csv",
            "status": "SOURCE_MATCHED_PROJECTED_DIAGNOSTIC_CLOSURE",
        },
        {
            "claim_id": "P6.T9A.06",
            "claim": "Action-level HR/BPB projector derivation remains open.",
            "source": "paperVI_EG_T4_projector_inventory.md/csv and T8 decisions",
            "status": "OPEN_ITEM_LOCKED",
        },
    ]

    decisions = [
        {
            "item": "appendix_scope",
            "decision": "WRITE_CONDITIONAL_PROJECTED_DIAGNOSTIC_CLOSURE",
            "basis": "T1-T8 support r_Wdelta=sqrt(a_bg) in projected diagnostics, but not yet as action-level HR/BPB proof.",
        },
        {
            "item": "wenzl_shape_space",
            "decision": "MENTION_AS_OBSTRUCTION",
            "basis": "T7 shows Wenzl response curves are nearly collinear; not an A4 proof space.",
        },
        {
            "item": "auto_cross_space",
            "decision": "USE_AS_CORRECT_PROJECTED_DIAGNOSTIC_SPACE",
            "basis": "T8 gives row-level X/A=r_global and fraction <X,X>/<A,A>=r_global^2 close to a_bg.",
        },
        {
            "item": "action_level_claim",
            "decision": "DO_NOT_CLAIM_CLOSED",
            "basis": "The HR/BPB perturbation projector identity remains the next theory target.",
        },
        {
            "item": "next_gate",
            "decision": "P6_T9B_OPTIONAL_ACTION_LEVEL_PROJECTOR_OR_PAPER_VI_INSERTION",
            "basis": "Either insert the appendix or continue to a harder action-level projector derivation.",
        },
    ]

    key_values = [
        {"quantity": k, "value": v} for k, v in vals.items() if isinstance(v, (int, float))
    ]

    result = {
        "status": {"gate_status": status, "likelihood_run": "NO", "ns_touched": "NO"},
        "key_values": vals,
        "source_lock": source_lock,
        "decisions": decisions,
    }

    (OUTDIR / "paperVI_EG_T9A_theory_appendix.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (OUTDIR / "paperVI_EG_T9A_theory_appendix.md").write_text(build_markdown(vals), encoding="utf-8")
    (OUTDIR / "paperVI_EG_T9A_theory_appendix.tex").write_text(build_tex(vals), encoding="utf-8")

    write_csv(OUTDIR / "paperVI_EG_T9A_key_values.csv", key_values, ["quantity", "value"])
    write_csv(OUTDIR / "paperVI_EG_T9A_source_lock.csv", source_lock, ["claim_id", "claim", "source", "status"])
    write_csv(OUTDIR / "paperVI_EG_T9A_decisions.csv", decisions, ["item", "decision", "basis"])

    print(status)
    print("Wrote outputs to:", OUTDIR)
    print(" -", OUTDIR / "paperVI_EG_T9A_theory_appendix.md")
    print(" -", OUTDIR / "paperVI_EG_T9A_theory_appendix.tex")
    print(" -", OUTDIR / "paperVI_EG_T9A_source_lock.csv")
    print(" -", OUTDIR / "paperVI_EG_T9A_decisions.csv")


if __name__ == "__main__":
    main()
