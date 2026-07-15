"""Planar OR spatial Newtonian N-body integrator for orbital-memory.

Units: G = 1. Dimension d (2 or 3) is inferred from the body positions, so the
same code runs the flat flip-flop and the inclined 3D Trojan. State layout:
[pos_0..pos_{N-1}, vel_0..vel_{N-1}] flattened, each vector length d.

Bodies may be massless (true restricted test particles): acceleration never
divides by a body's own mass, so m = 0 is a valid, well-behaved test particle
that feels the field but perturbs nothing.
"""

import numpy as np
from scipy.integrate import solve_ivp

G = 1.0


def _corotation_velocity(pos):
    """Omega x r for Omega = n z_hat, n = 1 (rigid rotation of the frame)."""
    v = np.zeros_like(pos)
    v[0] = -pos[1]
    v[1] = pos[0]
    return v  # z-component stays 0


def rhs(t, y, masses, dim, drag=None):
    n = len(masses)
    pos = y[: n * dim].reshape(n, dim)
    vel = y[n * dim :].reshape(n, dim)
    acc = np.zeros_like(pos)
    for i in range(n):
        dr = pos - pos[i]
        r2 = np.einsum("ij,ij->i", dr, dr)
        r2[i] = np.inf
        acc[i] = np.sum((G * masses / (r2 * np.sqrt(r2)))[:, None] * dr, axis=0)
    if drag is not None:
        k, idx = drag
        for i in idx:
            acc[i] -= k * (vel[i] - _corotation_velocity(pos[i]))
    return np.concatenate([vel.ravel(), acc.ravel()])


def total_energy(y, masses, dim):
    n = len(masses)
    pos = y[: n * dim].reshape(n, dim)
    vel = y[n * dim :].reshape(n, dim)
    kinetic = 0.5 * np.sum(masses * np.einsum("ij,ij->i", vel, vel))
    potential = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            potential -= G * masses[i] * masses[j] / np.linalg.norm(pos[i] - pos[j])
    return kinetic + potential


def integrate(bodies, t_end, n_samples=2000, rtol=1e-11, atol=1e-12, t0=0.0,
              drag=None):
    """Integrate a list of bodies (dicts with 'm','pos','vel') to t_end.

    Positions/velocities may be 2- or 3-vectors (all bodies must agree).
    drag=(k, [indices]) applies weak corotation drag to those bodies — an
    experimental non-gravitational term kept ONLY to demonstrate that
    dissipation destroys this memory (see TestDissipation).

    Returns times, trajectories (N, n_samples, d), velocities, the final full
    state vector, the dimension, and relative energy drift.
    NOTE: energy_drift weights bodies by mass, so it measures the MASSIVE
    bodies only — a massless test particle contributes nothing. The particle's
    own invariant is the Jacobi constant (theory.jacobi_of). Result keys here
    (traj/vel) intentionally differ from rotating.integrate's (xy/v); the two
    engines' results are not interchangeable.
    """
    masses = np.array([b["m"] for b in bodies], dtype=float)
    dim = len(bodies[0]["pos"])
    y0 = np.concatenate(
        [np.array([b["pos"] for b in bodies], float).ravel(),
         np.array([b["vel"] for b in bodies], float).ravel()]
    )
    times = np.linspace(t0, t_end, n_samples)
    sol = solve_ivp(rhs, (t0, t_end), y0, args=(masses, dim, drag),
                    method="DOP853", t_eval=times, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    n = len(masses)
    traj = sol.y[: n * dim].T.reshape(-1, n, dim).transpose(1, 0, 2)
    vels = sol.y[n * dim :].T.reshape(-1, n, dim).transpose(1, 0, 2)
    e0 = total_energy(y0, masses, dim)
    e1 = total_energy(sol.y[:, -1], masses, dim)
    drift = abs((e1 - e0) / e0) if e0 != 0 else abs(e1 - e0)
    return {"t": times, "traj": traj, "vel": vels, "yf": sol.y[:, -1],
            "masses": masses, "dim": dim, "energy_drift": drift}


def wrap_pi(a):
    return (a + np.pi) % (2 * np.pi) - np.pi
