#!/usr/bin/env python3
"""Complete-Lyapunov collocation for a toy tumour–effector ODE.

The construction follows the meshless orbital-derivative collocation of
Argáez, Giesl and Hafstein (Wendland kernel, speed-normalised field,
failure set where the derivative will not stay negative). It is a
numerical approximation on a fixed window, not a computer-assisted proof
and not a treatment model.

Seed 20260921 is recorded for the batched return map. The collocation
itself is deterministic.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921

# Declared toy. Time is scaled by the tumour growth rate.
# T is burden in units of carrying capacity. E is effector density
# in units of a reference density. Not estimated from a patient or a cell line.
A_PAR = 1.0
B_PAR = 0.5
C_PAR = 0.65
D_PAR = 0.2

# Collocation window. Fixed before the production run; contains the
# periodic orbit and the coexistence focus, and excludes the axial saddles.
T_LO, T_HI = 0.03, 0.58
E_LO, E_HI = 0.25, 0.92
N_PRIMARY = 34
N_COARSE = 22
N_FINE = 42
SUPPORT_MULT = 9.0
# Among the three cuts below, the run selects the value whose rebuilt
# function has the smallest peak-to-peak variation of V on the sampled
# orbit, subject to covering at least 95 percent of that orbit within
# two grid spacings. The other two cuts are still computed and stored.
GAMMA_MENU = (-0.55, -0.40, -0.25)
STENCIL_M = 4
STENCIL_STEP = 0.4
DELTA2 = 1e-8
N_BINARY_ITERS = 5


def field(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    t = z[..., 0]
    e = z[..., 1]
    kill = A_PAR * t * e / (B_PAR + t)
    dT = t * (1.0 - t) - kill
    dE = C_PAR * kill - D_PAR * e
    return np.stack([dT, dE], axis=-1)


def equilibrium() -> np.ndarray:
    x = D_PAR * B_PAR / (C_PAR * A_PAR - D_PAR)
    y = (1.0 - x) * (B_PAR + x) / A_PAR
    return np.array([x, y], dtype=float)


def jacobian(z: np.ndarray) -> np.ndarray:
    t, e = float(z[0]), float(z[1])
    b = B_PAR
    den = b + t
    dfT = (1.0 - 2.0 * t) - A_PAR * e * b / den**2
    dfE = -A_PAR * t / den
    dgT = C_PAR * A_PAR * e * b / den**2
    dgE = C_PAR * A_PAR * t / den - D_PAR
    return np.array([[dfT, dfE], [dgT, dgE]], dtype=float)


def phi_derivatives(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """First and second derivatives of unscaled Wendland psi_{4,2}."""
    s = np.asarray(s, dtype=float)
    d1 = np.zeros_like(s)
    d2 = np.zeros_like(s)
    mask = (s > 0.0) & (s < 1.0)
    u = 1.0 - s[mask]
    d1[mask] = -(u**5 / 5.0 - 11.0 * u**6 / 30.0 + u**7 / 6.0)
    d2[mask] = u**4 - (11.0 / 5.0) * u**5 + (7.0 / 6.0) * u**6
    return d1, d2


def psi_pair(dist: np.ndarray, shape: float) -> tuple[np.ndarray, np.ndarray]:
    """psi_1 and psi_2 in the Giesl radial-collocation formulae, support 1/shape."""
    dist = np.asarray(dist, dtype=float)
    p1 = np.zeros_like(dist)
    p2 = np.zeros_like(dist)
    zero = dist <= 1e-14
    p1[zero] = -(shape**2) / 30.0
    p2[zero] = shape**4
    positive = ~zero
    radius = np.zeros_like(dist)
    radius[positive] = shape * dist[positive]
    inside = positive & (radius < 1.0)
    if np.any(inside):
        sm = radius[inside]
        d1, d2 = phi_derivatives(sm)
        p1[inside] = shape**2 * d1 / sm
        p2[inside] = shape**4 * (sm * d2 - d1) / sm**3
    return p1, p2


def normalise(f: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(f, axis=-1)
    return f / np.sqrt(DELTA2 + norm**2)[..., None]


def assemble(points: np.ndarray, f_hat: np.ndarray, shape: float) -> np.ndarray:
    diff = points[:, None, :] - points[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    p1, p2 = psi_pair(dist, shape)
    fi = f_hat[:, None, :]
    fj = f_hat[None, :, :]
    dot_i = np.sum(diff * fi, axis=-1)
    dot_j = np.sum((-diff) * fj, axis=-1)
    dot_f = np.sum(fi * fj, axis=-1)
    return p2 * dot_i * dot_j - p1 * dot_f


def values_and_derivative(
    query: np.ndarray,
    f_query: np.ndarray,
    points: np.ndarray,
    f_hat: np.ndarray,
    beta: np.ndarray,
    shape: float,
) -> tuple[np.ndarray, np.ndarray]:
    diff = query[:, None, :] - points[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    p1, p2 = psi_pair(dist, shape)
    fj = f_hat[None, :, :]
    fq = f_query[:, None, :]
    orbital = -p1 * np.sum(fq * fj, axis=-1) + p2 * np.sum(diff * fq, axis=-1) * np.sum(
        (-diff) * fj, axis=-1
    )
    coef = np.sum((-diff) * fj, axis=-1) * p1
    return coef @ beta, orbital @ beta


def cartesian_grid(n: int) -> tuple[np.ndarray, float]:
    ts = np.linspace(T_LO, T_HI, n)
    es = np.linspace(E_LO, E_HI, n)
    grid = np.stack(np.meshgrid(ts, es, indexing="xy"), axis=-1).reshape(-1, 2)
    spacing = float(max(ts[1] - ts[0], es[1] - es[0]))
    return grid, spacing


def stencil(points: np.ndarray, f_hat: np.ndarray, spacing: float) -> np.ndarray:
    direction = f_hat / np.linalg.norm(f_hat, axis=1)[:, None]
    step = STENCIL_STEP * spacing
    blocks = []
    for k in range(1, STENCIL_M + 1):
        blocks.append(points + k * step * direction)
        blocks.append(points - k * step * direction)
    return np.stack(blocks, axis=1)


def stencil_means(
    blocks: np.ndarray,
    points: np.ndarray,
    f_hat: np.ndarray,
    beta: np.ndarray,
    shape: float,
) -> np.ndarray:
    means = np.zeros(len(points))
    counts = np.zeros(len(points))
    for k in range(blocks.shape[1]):
        sample = blocks[:, k, :]
        f_sample = normalise(field(sample))
        _, orbital = values_and_derivative(sample, f_sample, points, f_hat, beta, shape)
        means += orbital
        counts += 1.0
    return means / counts


def sample_orbit(eq: np.ndarray) -> tuple[np.ndarray, float]:
    def ode(_t, z):
        return field(np.asarray(z, dtype=float))

    sol = solve_ivp(
        ode,
        (0.0, 420.0),
        eq * np.array([1.12, 0.88]),
        rtol=1e-8,
        atol=1e-8,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    times = np.linspace(300.0, 420.0, 9000)
    path = sol.sol(times)
    series = path[0]
    crossings = np.where((series[:-1] < eq[0]) & (series[1:] >= eq[0]))[0]
    if len(crossings) < 3:
        raise RuntimeError("periodic orbit was not detected by the section")
    periods = np.diff(times[crossings])
    period = float(np.median(periods[-4:]))
    i0, i1 = crossings[-3], crossings[-2]
    orbit = path[:, i0 : i1 + 1].T
    return orbit, period


def self_test() -> None:
    rng = np.random.default_rng(0)
    centres = rng.normal(size=(10, 2)) * 0.2
    beta = rng.normal(size=10)
    shape = 1.7

    def f_test(x):
        return np.stack([-x[..., 1], x[..., 0]], axis=-1)

    f_centres = f_test(centres)
    x = np.array([[0.15, -0.08]])
    f = f_test(x)
    value, orbital = values_and_derivative(x, f, centres, f_centres, beta, shape)
    eps = 1e-6
    forward, _ = values_and_derivative(x + eps * f, f_test(x + eps * f), centres, f_centres, beta, shape)
    backward, _ = values_and_derivative(x - eps * f, f_test(x - eps * f), centres, f_centres, beta, shape)
    numeric = (forward - backward) / (2.0 * eps)
    rel = abs(float(orbital[0] - numeric[0])) / max(1.0, abs(float(numeric[0])))
    if rel > 1e-6:
        raise RuntimeError(f"orbital derivative disagrees with a finite difference ({rel})")
    # Collocation residual on a tiny non-equilibrium set.
    pts = np.array([[0.2, 0.3], [0.25, 0.45], [0.4, 0.35], [0.33, 0.55], [0.5, 0.4]])
    fh = normalise(field(pts))
    matrix = assemble(pts, fh, 2.0)
    coeff = solve(matrix, -np.ones(len(pts)), assume_a="pos")
    _, got = values_and_derivative(pts, fh, pts, fh, coeff, 2.0)
    if np.max(np.abs(got + 1.0)) > 1e-8:
        raise RuntimeError("collocation residual at the nodes is not zero")


def geometry(fail: np.ndarray, points: np.ndarray, orbit: np.ndarray, eq: np.ndarray, spacing: float) -> dict:
    targets = np.vstack([orbit[::2], eq[None, :]])
    dist_union = cKDTree(targets).query(points)[0]
    dist_orbit = cKDTree(orbit[::2]).query(points)[0]
    dist_eq = np.linalg.norm(points - eq, axis=1)
    cover = cKDTree(points[fail]).query(orbit[::2])[0] if np.any(fail) else np.array([np.inf])
    stray = fail & (dist_union > 0.10)
    near_orbit = fail & (dist_orbit <= dist_eq)
    near_focus = fail & (dist_eq < dist_orbit)
    return {
        "fail_fraction": float(np.mean(fail)),
        "fail_count": int(np.sum(fail)),
        "distance_union_median": float(np.median(dist_union[fail])) if np.any(fail) else None,
        "distance_union_p90": float(np.quantile(dist_union[fail], 0.9)) if np.any(fail) else None,
        "orbit_cover_within_1_5_spacing": float(np.mean(cover < 1.5 * spacing)),
        "orbit_cover_within_2_spacing": float(np.mean(cover < 2.0 * spacing)),
        "stray_fraction_of_fail": float(np.mean(dist_union[fail] > 0.10)) if np.any(fail) else None,
        "stray_count": int(np.sum(stray)),
        "near_orbit_fail_count": int(np.sum(near_orbit)),
        "near_focus_fail_count": int(np.sum(near_focus)),
        "near_orbit_tube_median": float(np.median(dist_orbit[near_orbit])) if np.any(near_orbit) else None,
        "near_orbit_tube_p90": float(np.quantile(dist_orbit[near_orbit], 0.9)) if np.any(near_orbit) else None,
        "focus_blob_median": float(np.median(dist_eq[near_focus])) if np.any(near_focus) else None,
        "focus_blob_max": float(np.max(dist_eq[near_focus])) if np.any(near_focus) else None,
        "focus_to_nearest_fail": float(np.min(dist_eq[fail])) if np.any(fail) else None,
    }


def orbit_flatness(beta, points, f_hat, shape, orbit) -> dict:
    sample = orbit[:: max(1, len(orbit) // 120)]
    values, orbital = values_and_derivative(sample, normalise(field(sample)), points, f_hat, beta, shape)
    return {
        "V_peak_to_peak": float(np.max(values) - np.min(values)),
        "V_net": float(values[-1] - values[0]),
        "V_mean": float(np.mean(values)),
        "orbital_median": float(np.median(orbital)),
        "orbital_p10": float(np.quantile(orbital, 0.1)),
        "orbital_p90": float(np.quantile(orbital, 0.9)),
        "orbital_max": float(np.max(orbital)),
        "fraction_positive": float(np.mean(orbital > 0.0)),
    }


def collocate(n: int, gamma: float, orbit: np.ndarray, eq: np.ndarray, record_iters: bool = False) -> dict:
    points, spacing = cartesian_grid(n)
    shape = 1.0 / (SUPPORT_MULT * spacing)
    f_hat = normalise(field(points))
    matrix = assemble(points, f_hat, shape)
    eig_min = float(np.min(np.linalg.eigvalsh(matrix)))
    cond = float(np.linalg.cond(matrix))
    beta0 = solve(matrix, -np.ones(len(points)), assume_a="pos")
    blocks = stencil(points, f_hat, spacing)
    score0 = stencil_means(blocks, points, f_hat, beta0, shape)
    fail = score0 > gamma
    rhs = np.where(fail, 0.0, -1.0)
    beta = solve(matrix, rhs, assume_a="pos")
    score = stencil_means(blocks, points, f_hat, beta, shape)
    node_orbital = matrix @ beta
    out = {
        "n": n,
        "n_points": int(len(points)),
        "spacing": spacing,
        "shape": shape,
        "support_radius": float(1.0 / shape),
        "condition_number": cond,
        "min_eigenvalue": eig_min,
        "gamma": gamma,
        "geometry": geometry(fail, points, orbit, eq, spacing),
        "flatness": orbit_flatness(beta, points, f_hat, shape, orbit),
        "stencil_median_fail": float(np.median(score[fail])) if np.any(fail) else None,
        "stencil_median_pass": float(np.median(score[~fail])) if np.any(~fail) else None,
        "node_orbital_max_abs_error": float(
            np.max(np.abs(node_orbital - rhs))
        ),
        "any_point_fail_fraction_first_solve": None,
    }
    # Any-point rule on the first solve, for the comparison sentence only.
    any_fail = np.zeros(len(points), dtype=bool)
    for k in range(blocks.shape[1]):
        sample = blocks[:, k, :]
        _, orbital = values_and_derivative(sample, normalise(field(sample)), points, f_hat, beta0, shape)
        any_fail |= orbital > gamma
    out["any_point_fail_fraction_first_solve"] = float(np.mean(any_fail))
    out["mean_rule_fail_fraction_first_solve"] = float(np.mean(fail))

    if record_iters:
        history = []
        beta_it = beta0.copy()
        fail_it = np.zeros(len(points), dtype=bool)
        for it in range(N_BINARY_ITERS):
            sc = stencil_means(blocks, points, f_hat, beta_it, shape)
            fail_it = sc > gamma
            beta_it = solve(matrix, np.where(fail_it, 0.0, -1.0), assume_a="pos")
            geo = geometry(fail_it, points, orbit, eq, spacing)
            flat = orbit_flatness(beta_it, points, f_hat, shape, orbit)
            history.append(
                {
                    "iteration": it + 1,
                    "fail_fraction": geo["fail_fraction"],
                    "distance_union_p90": geo["distance_union_p90"],
                    "orbit_cover_within_2_spacing": geo["orbit_cover_within_2_spacing"],
                    "V_peak_to_peak": flat["V_peak_to_peak"],
                }
            )
        out["binary_iterations"] = history

    out["_arrays"] = {
        "points": points,
        "fail": fail,
        "beta": beta,
        "beta0": beta0,
        "f_hat": f_hat,
        "shape": shape,
        "spacing": spacing,
        "score0": score0,
        "score": score,
    }
    return out


def transient_drop(beta, points, f_hat, shape, period: float) -> dict:
    def ode(_t, z):
        return field(np.asarray(z, dtype=float))

    start = np.array([0.55, 0.88])
    horizon = 10.0 * period
    sol = solve_ivp(ode, (0.0, horizon), start, rtol=1e-7, atol=1e-7, dense_output=True)
    times = np.linspace(0.0, horizon, 1200)
    path = sol.sol(times).T
    values = np.empty(len(path))
    for i in range(0, len(path), 150):
        sl = slice(i, i + 150)
        values[sl], _ = values_and_derivative(
            path[sl], normalise(field(path[sl])), points, f_hat, beta, shape
        )
    # Drop over the first period, and from start to the last quarter.
    i_per = int(np.searchsorted(times, period))
    return {
        "start": start.tolist(),
        "V_start": float(values[0]),
        "V_after_one_period": float(values[i_per]),
        "drop_one_period": float(values[0] - values[i_per]),
        "V_end": float(values[-1]),
        "drop_total": float(values[0] - values[-1]),
        "times": times,
        "values": values,
        "path": path,
    }


def return_diagnostic(points: np.ndarray, fail: np.ndarray, orbit: np.ndarray, eq: np.ndarray, period: float) -> dict:
    """Finite-time return: large excursion, and a close approach to the start near one period.

    Staying near a weakly unstable focus is not counted as a return.
    """
    dt = period / 250.0
    steps = 500  # two periods
    z = points.copy()
    z0 = points.copy()
    max_exc = np.zeros(len(points))
    dist_at_period = np.full(len(points), np.inf)
    left = np.zeros(len(points), dtype=bool)
    margin_lo = np.array([T_LO - 0.02, E_LO - 0.02])
    margin_hi = np.array([T_HI + 0.02, E_HI + 0.02])
    for step in range(steps):
        k1 = field(z)
        k2 = field(z + 0.5 * dt * k1)
        k3 = field(z + 0.5 * dt * k2)
        k4 = field(z + dt * k3)
        z = z + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        time = (step + 1) * dt
        left |= np.any((z < margin_lo) | (z > margin_hi), axis=1)
        dist = np.linalg.norm(z - z0, axis=1)
        max_exc = np.maximum(max_exc, dist)
        if 0.85 * period <= time <= 1.15 * period:
            dist_at_period = np.minimum(dist_at_period, dist)
    returns = (max_exc >= 0.12) & (dist_at_period < 0.06)
    dist_orbit = cKDTree(orbit[::2]).query(points)[0]
    dist_eq = np.linalg.norm(points - eq, axis=1)
    tube = fail & (dist_orbit < 0.06) & (dist_orbit <= dist_eq)
    blob = fail & (dist_eq < 0.06) & (dist_eq < dist_orbit)
    return {
        "excursion_floor": 0.12,
        "return_radius": 0.06,
        "window_exit_fraction": float(np.mean(left)),
        "window_exit_fraction_fail": float(np.mean(left[fail])) if np.any(fail) else None,
        "window_exit_fraction_pass": float(np.mean(left[~fail])) if np.any(~fail) else None,
        "return_fraction_all": float(np.mean(returns)),
        "return_fraction_fail": float(np.mean(returns[fail])) if np.any(fail) else None,
        "return_fraction_pass": float(np.mean(returns[~fail])) if np.any(~fail) else None,
        "return_fraction_orbit_tube": float(np.mean(returns[tube])) if np.any(tube) else None,
        "return_fraction_focus_blob": float(np.mean(returns[blob])) if np.any(blob) else None,
        "orbit_tube_count": int(np.sum(tube)),
        "focus_blob_count": int(np.sum(blob)),
        "returns": returns,
    }


def phase_portrait(orbit, eq, period, path) -> None:
    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    ts = np.linspace(0.0, 1.05, 400)
    e_null_T = (1.0 - ts) * (B_PAR + ts) / A_PAR
    ax.plot(ts, e_null_T, color="#555555", lw=1.0, label="tumour nullcline")
    ax.axvline(eq[0], color="#555555", lw=1.0, ls="--", label="effector nullcline")
    # a few transients
    def ode(_t, z):
        return field(np.asarray(z, dtype=float))

    starts = [
        np.array([0.08, 0.80]),
        np.array([0.50, 0.30]),
        np.array([0.12, 0.32]),
        np.array([0.48, 0.86]),
    ]
    for i, z0 in enumerate(starts):
        sol = solve_ivp(ode, (0.0, 2.2 * period), z0, rtol=1e-6, atol=1e-6, dense_output=True)
        tt = np.linspace(0.0, 2.2 * period, 400)
        zz = sol.sol(tt)
        ax.plot(zz[0], zz[1], color="#4C78A8", lw=0.9, alpha=0.9, label="transient" if i == 0 else None)
    ax.plot(orbit[:, 0], orbit[:, 1], color="#E45756", lw=1.8, label="periodic orbit")
    ax.scatter([eq[0]], [eq[1]], c="black", s=28, zorder=5, label="coexistence focus")
    ax.set_xlim(0.0, 0.70)
    ax.set_ylim(0.15, 1.05)
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "phase_portrait.png", dpi=160)
    plt.close(fig)


def derivative_panels(pack, orbit, eq) -> None:
    arr = pack["_arrays"]
    points, shape, f_hat = arr["points"], arr["shape"], arr["f_hat"]
    ng = 160
    gt, ge = np.meshgrid(
        np.linspace(T_LO, T_HI, ng),
        np.linspace(E_LO, E_HI, ng),
        indexing="xy",
    )
    grid = np.stack([gt.ravel(), ge.ravel()], axis=1)
    panels = []
    for beta in (arr["beta0"], arr["beta"]):
        orbital = np.empty(len(grid))
        for i in range(0, len(grid), 400):
            sl = slice(i, i + 400)
            _, orbital[sl] = values_and_derivative(
                grid[sl], normalise(field(grid[sl])), points, f_hat, beta, shape
            )
        panels.append(orbital.reshape(ng, ng))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3), sharey=True)
    titles = ("first solve, target −1", "rebuild, target 0 on the failing set")
    for ax, data, title in zip(axes, panels, titles):
        mesh = ax.pcolormesh(gt, ge, data, cmap="coolwarm", vmin=-1.4, vmax=0.6, shading="auto")
        ax.contour(gt, ge, data, levels=[-0.4], colors="black", linewidths=0.7)
        ax.plot(orbit[:, 0], orbit[:, 1], color="black", lw=1.15)
        ax.scatter([eq[0]], [eq[1]], c="black", s=16, zorder=5)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("tumour burden T")
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("effector density E")
    fig.colorbar(mesh, ax=axes, fraction=0.046, pad=0.04, label="orbital derivative along the normalised field")
    fig.savefig(FIG / "orbital_derivative.png", dpi=160)
    plt.close(fig)


def partition_plot(pack, orbit, eq) -> None:
    arr = pack["_arrays"]
    points, fail = arr["points"], arr["fail"]
    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    ax.scatter(
        points[~fail, 0],
        points[~fail, 1],
        s=9,
        c="#4C78A8",
        linewidths=0,
        label="decrease nodes (stencil mean ≤ −0.4)",
    )
    ax.scatter(
        points[fail, 0],
        points[fail, 1],
        s=11,
        c="#E45756",
        linewidths=0,
        label="near-level nodes (stencil mean > −0.4)",
    )
    ax.plot(orbit[:, 0], orbit[:, 1], color="black", lw=1.3, label="periodic orbit")
    ax.scatter([eq[0]], [eq[1]], c="black", s=28, zorder=5, label="coexistence focus")
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "partition.png", dpi=160)
    plt.close(fig)


def level_plot(pack, orbit, eq) -> None:
    arr = pack["_arrays"]
    points, shape, f_hat, beta = arr["points"], arr["shape"], arr["f_hat"], arr["beta"]
    ng = 160
    gt, ge = np.meshgrid(
        np.linspace(T_LO, T_HI, ng),
        np.linspace(E_LO, E_HI, ng),
        indexing="xy",
    )
    grid = np.stack([gt.ravel(), ge.ravel()], axis=1)
    values = np.empty(len(grid))
    for i in range(0, len(grid), 400):
        sl = slice(i, i + 400)
        values[sl], _ = values_and_derivative(
            grid[sl], normalise(field(grid[sl])), points, f_hat, beta, shape
        )
    values = values.reshape(ng, ng)
    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    fill = ax.contourf(gt, ge, values, levels=18, cmap="viridis")
    ax.contour(gt, ge, values, levels=12, colors="black", linewidths=0.35, alpha=0.55)
    ax.plot(orbit[:, 0], orbit[:, 1], color="white", lw=1.5)
    ax.scatter([eq[0]], [eq[1]], c="white", s=22, zorder=5)
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(fill, ax=ax, fraction=0.046, pad=0.04, label="rebuilt collocation function V")
    fig.tight_layout()
    fig.savefig(FIG / "level_sets.png", dpi=160)
    plt.close(fig)


def decrease_plot(transient: dict, flat: dict, orbit, pack) -> None:
    arr = pack["_arrays"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    ax = axes[0]
    band_lo = flat["V_mean"] - 0.5 * flat["V_peak_to_peak"]
    band_hi = flat["V_mean"] + 0.5 * flat["V_peak_to_peak"]
    ax.axhspan(band_lo, band_hi, color="#E45756", alpha=0.25, label="range of V on the periodic orbit")
    ax.plot(transient["times"], transient["values"], color="#4C78A8", lw=1.4, label="approach from (0.55, 0.88)")
    ax.set_xlabel("time")
    ax.set_ylabel("V")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    score = arr["score"]
    fail = arr["fail"]
    bins = np.linspace(-1.8, 0.8, 36)
    ax.hist(score[~fail], bins=bins, color="#4C78A8", alpha=0.85, label="decrease nodes")
    ax.hist(score[fail], bins=bins, color="#E45756", alpha=0.75, label="near-level nodes")
    ax.axvline(pack["gamma"], color="black", lw=0.8, ls="--")
    ax.set_xlabel("stencil-mean orbital derivative after the rebuild")
    ax.set_ylabel("collocation nodes")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "decrease_vs_level.png", dpi=160)
    plt.close(fig)
    # silence unused in case flat is only for the caller
    _ = flat


def return_plot(points, fail, returns, orbit, eq) -> None:
    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    both = fail & returns
    fail_only = fail & ~returns
    return_only = ~fail & returns
    neither = ~fail & ~returns
    ax.scatter(points[neither, 0], points[neither, 1], s=8, c="#d0d0d0", linewidths=0, label="neither")
    ax.scatter(points[return_only, 0], points[return_only, 1], s=10, c="#54A24B", linewidths=0, label="returns, not a failing node")
    ax.scatter(points[fail_only, 0], points[fail_only, 1], s=11, c="#E45756", linewidths=0, label="failing node, does not return")
    ax.scatter(points[both, 0], points[both, 1], s=11, c="#F2C14E", linewidths=0, label="failing node and returns")
    ax.plot(orbit[:, 0], orbit[:, 1], color="black", lw=1.1)
    ax.scatter([eq[0]], [eq[1]], c="black", s=22, zorder=5)
    ax.set_xlabel("tumour burden T")
    ax.set_ylabel("effector density E")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "return_vs_failure.png", dpi=160)
    plt.close(fig)


def public_pack(pack: dict) -> dict:
    return {k: v for k, v in pack.items() if k != "_arrays"}


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    self_test()
    eq = equilibrium()
    jac = jacobian(eq)
    eig = np.linalg.eigvals(jac)
    orbit, period = sample_orbit(eq)
    linear_period = float(2.0 * np.pi / np.abs(np.imag(eig[0])))

    menu = [collocate(N_PRIMARY, g, orbit, eq, record_iters=False) for g in GAMMA_MENU]

    def selectable(pack: dict) -> bool:
        return pack["geometry"]["orbit_cover_within_2_spacing"] >= 0.95

    eligible = [pack for pack in menu if selectable(pack)]
    if not eligible:
        raise RuntimeError("no cut covered the periodic orbit")
    primary_pick = min(eligible, key=lambda pack: pack["flatness"]["V_peak_to_peak"])
    primary = collocate(N_PRIMARY, primary_pick["gamma"], orbit, eq, record_iters=True)
    coarse = collocate(N_COARSE, primary["gamma"], orbit, eq, record_iters=False)
    fine = collocate(N_FINE, primary["gamma"], orbit, eq, record_iters=False)
    variants = [pack for pack in menu if pack["gamma"] != primary["gamma"]]

    arr = primary["_arrays"]
    transient = transient_drop(arr["beta"], arr["points"], arr["f_hat"], arr["shape"], period)
    returned = return_diagnostic(arr["points"], arr["fail"], orbit, eq, period)

    phase_portrait(orbit, eq, period, transient["path"])
    derivative_panels(primary, orbit, eq)
    partition_plot(primary, orbit, eq)
    level_plot(primary, orbit, eq)
    decrease_plot(transient, primary["flatness"], orbit, primary)
    return_plot(arr["points"], arr["fail"], returned["returns"], orbit, eq)

    # Ratio that Chapter 4 quotes: cycle variation against the one-period transient drop.
    ratio = primary["flatness"]["V_peak_to_peak"] / transient["drop_one_period"]

    def strip_variant(pack):
        data = public_pack(pack)
        return data

    results = {
        "seed": SEED,
        "parameters": {"a": A_PAR, "b": B_PAR, "c": C_PAR, "d": D_PAR},
        "window": {"T": [T_LO, T_HI], "E": [E_LO, E_HI]},
        "equilibrium": {"T": float(eq[0]), "E": float(eq[1])},
        "jacobian": jac.tolist(),
        "eigenvalues_real": [float(np.real(w)) for w in eig],
        "eigenvalues_imag": [float(np.imag(w)) for w in eig],
        "linear_period": linear_period,
        "nonlinear_period": period,
        "orbit_box": {
            "T_min": float(orbit[:, 0].min()),
            "T_max": float(orbit[:, 0].max()),
            "E_min": float(orbit[:, 1].min()),
            "E_max": float(orbit[:, 1].max()),
        },
        "method": {
            "wendland": "psi_4_2",
            "delta2": DELTA2,
            "support_multiple_of_spacing": SUPPORT_MULT,
            "gamma_menu": list(GAMMA_MENU),
            "gamma_selected": primary["gamma"],
            "gamma_rule": "smallest peak-to-peak of V on the orbit among cuts covering at least 95 percent of the orbit within two spacings",
            "stencil_each_side": STENCIL_M,
            "stencil_step_in_spacings": STENCIL_STEP,
            "decision": "mean orbital derivative on the flow-aligned stencil",
            "rebuild": "one solve with target 0 on the exceedance set and -1 off it",
        },
        "primary": strip_variant(primary),
        "coarse": strip_variant(coarse),
        "fine": strip_variant(fine),
        "gamma_variants": [strip_variant(v) for v in variants],
        "transient": {
            "start": transient["start"],
            "V_start": transient["V_start"],
            "V_after_one_period": transient["V_after_one_period"],
            "drop_one_period": transient["drop_one_period"],
            "V_end": transient["V_end"],
            "drop_total": transient["drop_total"],
        },
        "flatness_over_one_period_drop": ratio,
        "return": {k: v for k, v in returned.items() if k not in {"returns", "min_after"}},
    }
    (ROOT / "results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({
        "equilibrium": results["equilibrium"],
        "eigenvalues_real": results["eigenvalues_real"],
        "eigenvalues_imag": results["eigenvalues_imag"],
        "periods": [linear_period, period],
        "orbit_box": results["orbit_box"],
        "primary_geometry": results["primary"]["geometry"],
        "primary_flatness": results["primary"]["flatness"],
        "stencil_medians": [
            results["primary"]["stencil_median_fail"],
            results["primary"]["stencil_median_pass"],
        ],
        "any_vs_mean": [
            results["primary"]["any_point_fail_fraction_first_solve"],
            results["primary"]["mean_rule_fail_fraction_first_solve"],
        ],
        "transient": results["transient"],
        "ratio": ratio,
        "return": results["return"],
        "coarse_geometry": results["coarse"]["geometry"],
        "fine_geometry": results["fine"]["geometry"],
        "variants": [
            {"gamma": v["gamma"], "geometry": v["geometry"], "flatness": v["flatness"]}
            for v in results["gamma_variants"]
        ],
        "iters": results["primary"].get("binary_iterations"),
        "cond": [
            results["primary"]["condition_number"],
            results["coarse"]["condition_number"],
            results["fine"]["condition_number"],
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
