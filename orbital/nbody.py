"""Planar Newtonian N-body integrator for orbital-memory.

Units: G = 1. State layout: [x1,y1,...,xN,yN, vx1,vy1,...,vxN,vyN].
Bodies may be massless (true restricted test particles): acceleration never
divides by a body's own mass, so m = 0 is a valid, well-behaved test particle
that feels the field but perturbs nothing.
"""

import numpy as np
from scipy.integrate import solve_ivp

G = 1.0


def rhs(t, y, masses, drag=None):
    n = len(masses)
    pos = y[: 2 * n].reshape(n, 2)
    vel = y[2 * n :].reshape(n, 2)
    acc = np.zeros_like(pos)
    for i in range(n):
        dr = pos - pos[i]
        r2 = np.einsum("ij,ij->i", dr, dr)
        r2[i] = np.inf
        acc[i] = np.sum((G * masses / (r2 * np.sqrt(r2)))[:, None] * dr, axis=0)
    if drag is not None:
        # weak drag toward corotation (Omega = z_hat, n=1): damps a body's
        # velocity in the rotating frame, so it relaxes toward a Lagrange point.
        # Physically the gas-drag / tidal channel that captures Trojans.
        k, idx = drag
        for i in idx:
            v_rot = vel[i] - np.array([-pos[i, 1], pos[i, 0]])  # v - Omega x r
            acc[i] -= k * v_rot
    return np.concatenate([vel.ravel(), acc.ravel()])


def total_energy(y, masses):
    n = len(masses)
    pos = y[: 2 * n].reshape(n, 2)
    vel = y[2 * n :].reshape(n, 2)
    kinetic = 0.5 * np.sum(masses * np.einsum("ij,ij->i", vel, vel))
    potential = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            potential -= G * masses[i] * masses[j] / np.linalg.norm(pos[i] - pos[j])
    return kinetic + potential


def integrate(bodies, t_end, n_samples=2000, rtol=1e-11, atol=1e-12, t0=0.0,
              drag=None):
    """Integrate a list of bodies (dicts with 'm','pos','vel') to t_end.

    drag=(k, [indices]) applies weak corotation drag to those bodies (turns the
    conservative cell into an attracting one — capture, refresh, relaxation).

    Returns times (n_samples,), trajectories (N, n_samples, 2), velocities,
    the final full state vector, and relative energy drift.
    """
    masses = np.array([b["m"] for b in bodies], dtype=float)
    y0 = np.concatenate(
        [np.array([b["pos"] for b in bodies], float).ravel(),
         np.array([b["vel"] for b in bodies], float).ravel()]
    )
    times = np.linspace(t0, t_end, n_samples)
    sol = solve_ivp(rhs, (t0, t_end), y0, args=(masses, drag), method="DOP853",
                    t_eval=times, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    n = len(masses)
    traj = sol.y[: 2 * n].T.reshape(-1, n, 2).transpose(1, 0, 2)
    vels = sol.y[2 * n :].T.reshape(-1, n, 2).transpose(1, 0, 2)
    e0 = total_energy(y0, masses)
    e1 = total_energy(sol.y[:, -1], masses)
    drift = abs((e1 - e0) / e0) if e0 != 0 else abs(e1 - e0)
    return {"t": times, "traj": traj, "vel": vels, "yf": sol.y[:, -1],
            "masses": masses, "energy_drift": drift}


def longitude(pos, center=(0.0, 0.0)):
    """True longitude (radians) of pos about center, in (-pi, pi]."""
    d = np.asarray(pos) - np.asarray(center)
    return np.arctan2(d[..., 1], d[..., 0])


def wrap_pi(a):
    return (a + np.pi) % (2 * np.pi) - np.pi
