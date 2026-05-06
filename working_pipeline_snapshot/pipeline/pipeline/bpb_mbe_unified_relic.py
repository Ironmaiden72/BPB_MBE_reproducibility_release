
from __future__ import annotations

import math
from typing import Any

import numpy as np

C_LIGHT = 299792.458  # km/s


def activation_tanh(a: np.ndarray | float, a_bg: float, Delta_bg: float) -> np.ndarray | float:
    if not (0.0 < a_bg < 1.0):
        raise ValueError("a_bg must lie in (0,1).")
    if Delta_bg <= 0.0:
        raise ValueError("Delta_bg must be > 0.")
    a_arr = np.asarray(a, dtype=float)
    x = (np.log(a_arr) - math.log(a_bg)) / Delta_bg
    out = 0.5 * (1.0 + np.tanh(x))
    if np.isscalar(a):
        return float(out)
    return out


def F_bg(a: np.ndarray | float, a_bg: float, Delta_bg: float) -> np.ndarray | float:
    s = activation_tanh(a, a_bg=a_bg, Delta_bg=Delta_bg)
    s0 = float(activation_tanh(1.0, a_bg=a_bg, Delta_bg=Delta_bg))
    if s0 <= 0.0:
        raise ValueError("S_bg(1) must be > 0.")
    out = np.asarray(s, dtype=float) / s0
    if np.isscalar(a):
        return float(out)
    return out


def interaction_integral_Jbf(
    *,
    a_bg: float,
    Delta_bg: float,
    a_int_min: float = 1.0e-3,
    n_int: int = 4096,
) -> float:
    if not (0.0 < a_int_min < 1.0):
        raise ValueError("a_int_min must lie in (0,1).")
    if n_int < 8:
        raise ValueError("n_int must be >= 8.")
    ln_a = np.linspace(math.log(a_int_min), 0.0, int(n_int), dtype=float)
    a = np.exp(ln_a)
    F = np.asarray(F_bg(a, a_bg=a_bg, Delta_bg=Delta_bg), dtype=float)
    J = np.trapz(F, ln_a)
    return float(J)


def relic_split_from_params(
    *,
    Om_eff0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
    n_int: int = 4096,
) -> dict[str, float]:
    if Om_eff0 < 0.0:
        raise ValueError("Om_eff0 must be >= 0.")
    if kappa_bf < 0.0:
        raise ValueError("kappa_bf must be >= 0.")
    Ode_eff0 = 1.0 - Or0 - Ok0 - Om_eff0
    if Ode_eff0 < 0.0:
        raise ValueError(
            f"Non-positive DE-like remainder: 1-Or0-Ok0-Om_eff0={Ode_eff0:.6g}."
        )

    J_bf = interaction_integral_Jbf(
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        a_int_min=a_int_min,
        n_int=n_int,
    )
    u = float(kappa_bf) * max(-float(beta_c), 0.0) * float(J_bf)
    s_surv = math.exp(-u)
    Om_surv0 = float(Om_eff0) * s_surv
    Om_induced0 = float(Om_eff0) - Om_surv0

    return {
        "J_bf": float(J_bf),
        "u_bf": float(u),
        "s_surv": float(s_surv),
        "Om_surv0": float(Om_surv0),
        "Omega_m_induced0": float(Om_induced0),
        "Omega_m_eff0": float(Om_eff0),
        "Omega_de_eff0": float(Ode_eff0),
        "dm_survival_fraction": float(s_surv) if Om_eff0 > 0 else 0.0,
        "dm_induced_fraction": float(1.0 - s_surv) if Om_eff0 > 0 else 0.0,
    }


def E2_unified_relic(
    z: np.ndarray,
    *,
    Om_eff0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    split = relic_split_from_params(
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    F = np.asarray(F_bg(a, a_bg=a_bg, Delta_bg=Delta_bg), dtype=float)
    Om_surv0 = split["Om_surv0"]
    Om_induced0 = split["Omega_m_induced0"]
    Ode_eff0 = split["Omega_de_eff0"]
    return (
        Or0 * a ** (-4.0)
        + Ok0 * a ** (-2.0)
        + Om_surv0 * a ** (-3.0)
        + F * (Om_induced0 * a ** (-3.0) + Ode_eff0)
    )


def E_unified_relic(
    z: np.ndarray,
    *,
    Om_eff0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
) -> np.ndarray:
    E2 = E2_unified_relic(
        z,
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    if np.any(~np.isfinite(E2)) or np.any(E2 <= 0.0):
        raise ValueError("Invalid E^2 in E_unified_relic.")
    return np.sqrt(E2)


def H_unified_relic_H0(
    z: np.ndarray,
    *,
    Om_eff0: float,
    H0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
) -> np.ndarray:
    return H0 * E_unified_relic(
        z,
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )


def DM_mpc_unified_relic(
    z: np.ndarray,
    *,
    Om_eff0: float,
    H0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
    n_int: int = 4000,
) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    zmax = float(np.max(z))
    zg = np.linspace(0.0, zmax, int(n_int), dtype=float)

    Ez = E_unified_relic(
        zg,
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    chi = np.zeros_like(zg)
    chi[1:] = np.cumsum(0.5 * (1.0 / Ez[1:] + 1.0 / Ez[:-1]) * np.diff(zg))
    chi = (C_LIGHT / H0) * chi
    chi_of_z = np.interp(z, zg, chi)

    if abs(Ok0) < 1.0e-12:
        return chi_of_z
    if Ok0 > 0.0:
        sq = np.sqrt(Ok0)
        return (C_LIGHT / H0) / sq * np.sinh(sq * H0 * chi_of_z / C_LIGHT)

    sq = np.sqrt(-Ok0)
    return (C_LIGHT / H0) / sq * np.sin(sq * H0 * chi_of_z / C_LIGHT)


def mu_theory_from_DM(DM_mpc: np.ndarray, z: np.ndarray) -> np.ndarray:
    DL_mpc = (1.0 + np.asarray(z, dtype=float)) * np.asarray(DM_mpc, dtype=float)
    return 5.0 * np.log10(DL_mpc) + 25.0


def omega_m_source_a(
    a: np.ndarray,
    *,
    Om_eff0: float,
    kappa_bf: float,
    beta_c: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    z = 1.0 / a - 1.0
    split = relic_split_from_params(
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    F = np.asarray(F_bg(a, a_bg=a_bg, Delta_bg=Delta_bg), dtype=float)
    Om_surv0 = split["Om_surv0"]
    Om_induced0 = split["Omega_m_induced0"]
    num = Om_surv0 * a ** (-3.0) + Om_induced0 * F * a ** (-3.0)
    E2 = E2_unified_relic(
        z,
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    return np.divide(num, E2, out=np.zeros_like(num, dtype=float), where=(E2 > 0.0))


def integrate_growth_mu_const_dynamic(
    *,
    Om_eff0: float,
    kappa_bf: float,
    beta_c: float,
    H0: float,
    a_bg: float,
    Delta_bg: float,
    Ok0: float = 0.0,
    Or0: float = 0.0,
    a_int_min: float = 1.0e-3,
    zmax_cov: float = 100.0,
    a_ini: float = 1.0e-3,
    n_a: int = 2500,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = 1.0 + float(beta_c)
    if mu <= 0.0:
        raise ValueError(f"mu=1+beta_c must be > 0, got {mu}")

    a_ini_eff = max(float(a_ini), 1.0 / (1.0 + float(zmax_cov)))
    a = np.geomspace(a_ini_eff, 1.0, int(n_a))
    z = 1.0 / a - 1.0

    H = H_unified_relic_H0(
        z,
        Om_eff0=Om_eff0,
        H0=H0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )
    if np.any(~np.isfinite(H)) or np.any(H <= 0.0):
        raise FloatingPointError("Invalid H(z) in integrate_growth_mu_const_dynamic.")

    ln_a = np.log(a)
    ln_H = np.log(H)
    dlnH_dlnA = np.gradient(ln_H, ln_a)
    Om_src = omega_m_source_a(
        a,
        Om_eff0=Om_eff0,
        kappa_bf=kappa_bf,
        beta_c=beta_c,
        a_bg=a_bg,
        Delta_bg=Delta_bg,
        Ok0=Ok0,
        Or0=Or0,
        a_int_min=a_int_min,
    )

    y1 = np.zeros_like(a)
    y2 = np.zeros_like(a)
    y1[0] = a[0]
    y2[0] = a[0]

    for i in range(len(a) - 1):
        h = ln_a[i + 1] - ln_a[i]

        def rhs(i0: int, D: float, Dp: float) -> tuple[float, float]:
            A = 2.0 + dlnH_dlnA[i0]
            B = 1.5 * Om_src[i0] * mu
            return Dp, (-A * Dp + B * D)

        D, Dp = y1[i], y2[i]
        k1_1, k1_2 = rhs(i, D, Dp)
        k2_1, k2_2 = rhs(i, D + 0.5 * h * k1_1, Dp + 0.5 * h * k1_2)
        k3_1, k3_2 = rhs(i, D + 0.5 * h * k2_1, Dp + 0.5 * h * k2_2)
        k4_1, k4_2 = rhs(i, D + h * k3_1, Dp + h * k3_2)

        y1[i + 1] = D + (h / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1)
        y2[i + 1] = Dp + (h / 6.0) * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2)

        if not np.isfinite(y1[i + 1]) or not np.isfinite(y2[i + 1]):
            raise FloatingPointError(
                f"Growth integration diverged at step {i}, a={a[i]:.5e}, z={z[i]:.5e}"
            )

    D = y1.copy()
    D0 = float(D[-1])
    if not np.isfinite(D0) or D0 == 0.0:
        raise FloatingPointError("Invalid D(0) normalization.")

    D /= D0
    y2_norm = y2 / D0
    f = y2_norm / np.clip(D, 1.0e-30, np.inf)

    if np.any(~np.isfinite(D)) or np.any(~np.isfinite(f)):
        raise FloatingPointError("Non-finite D(z) or f(z) after normalization.")

    return z[::-1], D[::-1], f[::-1]


def _walk_items(obj: Any):
    stack = [("", obj)]
    while stack:
        path, cur = stack.pop()
        if isinstance(cur, dict):
            for k, v in cur.items():
                p2 = f"{path}.{k}" if path else str(k)
                stack.append((p2, v))
        elif isinstance(cur, list):
            for i, v in enumerate(cur):
                p2 = f"{path}[{i}]"
                stack.append((p2, v))
        else:
            yield path, cur


def find_value_recursive(j: dict, key: str) -> float | None:
    if key in j:
        try:
            return float(j[key])
        except Exception:
            pass
    if "best" in j and isinstance(j["best"], dict) and key in j["best"]:
        try:
            return float(j["best"][key])
        except Exception:
            pass
    for path, val in _walk_items(j):
        if path.lower().endswith(f".{key.lower()}"):
            try:
                return float(val)
            except Exception:
                continue
    return None
