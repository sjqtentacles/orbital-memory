"""A gravitational memory bit: which Lagrange island a test particle librates in.

Circular restricted three-body problem (G=1, total mass 1, separation 1, mean
motion n=1, period 2*pi). A heavy star and a planet orbit their barycenter;
a massless test particle sits in a tadpole orbit around the leading (L4) or
trailing (L5) triangular Lagrange point. Those two islands are separated by a
separatrix (the horseshoe region around L3), so a small perturbation can't move
the particle between them — the island it librates in is a nonvolatile bit.

Encoding: librating around L4 (+60 deg ahead of the planet) = '1';
librating around L5 (-60 deg, trailing) = '0'; a horseshoe that wraps past L3,
or an escape, = an erased/invalid cell.
"""

import numpy as np

from . import nbody

MU = 0.003          # planet mass fraction (< 0.0385 -> L4/L5 linearly stable)
PERIOD = 2 * np.pi  # planet orbital period


def _rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


def primaries(mu=MU):
    star = {"m": 1 - mu, "pos": [-mu, 0.0], "vel": [0.0, -mu]}
    planet = {"m": mu, "pos": [1 - mu, 0.0], "vel": [0.0, 1 - mu]}
    return star, planet


def lagrange_state(mu, point):
    """Position and corotating velocity of the exact L4/L5 point (inertial
    frame, t=0). Corotation velocity is z_hat x r with n = 1."""
    sy = 1.0 if point == "L4" else -1.0
    pos = np.array([0.5 - mu, sy * np.sqrt(3) / 2])
    vel = np.array([-pos[1], pos[0]])  # z_hat x pos, n=1
    return pos, vel


def make_cell(state="L4", mu=MU, libration_deg=6.0):
    """[star, planet, test particle] with the particle seeded into a tadpole
    of roughly `libration_deg` amplitude around the chosen island."""
    star, planet = primaries(mu)
    pos, vel = lagrange_state(mu, state)
    R = _rot(np.radians(libration_deg))  # nudge along the corotation circle
    pos, vel = R @ pos, R @ vel
    particle = {"m": 0.0, "pos": pos.tolist(), "vel": vel.tolist()}
    return [star, planet, particle]


def make_cell_3d(state="L4", mu=MU, libration_deg=6.0, inclination_z=0.16):
    """3D inclined Trojan: the in-plane tadpole of make_cell, plus the particle
    lifted to height z = inclination_z with zero vertical velocity, so it bobs
    through the orbital plane once per orbit (vertical epicyclic motion) while
    slowly librating around L4/L5. Star and planet stay in the z=0 plane."""
    star, planet = primaries(mu)
    for b in (star, planet):
        b["pos"] = b["pos"] + [0.0]
        b["vel"] = b["vel"] + [0.0]
    pos, vel = lagrange_state(mu, state)
    R = _rot(np.radians(libration_deg))
    pos, vel = R @ pos, R @ vel
    particle = {"m": 0.0,
                "pos": [pos[0], pos[1], inclination_z],
                "vel": [vel[0], vel[1], 0.0]}
    return [star, planet, particle]


def resonant_angle(res, particle=2, planet=1):
    """phi(t) = longitude(particle) - longitude(planet), degrees in (-180,180].
    Barycenter is the origin by construction."""
    lam_p = np.arctan2(res["traj"][particle, :, 1], res["traj"][particle, :, 0])
    lam_pl = np.arctan2(res["traj"][planet, :, 1], res["traj"][planet, :, 0])
    return np.degrees(nbody.wrap_pi(lam_p - lam_pl))


def _circular_mean_deg(phi):
    a = np.radians(phi)
    return np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean()))


def classify(phi):
    """Read the bit from a resonant-angle time series.

    Returns (label, center_deg, amplitude_deg). A tadpole stays on one side of
    the star-planet line and never crosses the L3 point at 180 deg; a horseshoe
    sweeps across it.
    """
    crosses_L3 = np.any(np.abs(np.diff(np.unwrap(np.radians(phi)))) > np.pi / 2)
    center = _circular_mean_deg(phi)
    amp = float(np.max(np.abs(nbody.wrap_pi(np.radians(phi) - np.radians(center))))
                * 180 / np.pi)
    span = phi.max() - phi.min()
    if crosses_L3 or span > 150:
        return "erased", center, amp
    return ("L4" if center > 0 else "L5"), center, amp


def hold(state="L4", periods=40, mu=MU, libration_deg=6.0, n_per_period=60):
    """Integrate a freshly written cell and return the run + resonant angle."""
    bodies = make_cell(state, mu, libration_deg)
    t_end = periods * PERIOD
    res = nbody.integrate(bodies, t_end, n_samples=int(periods * n_per_period))
    res["phi"] = resonant_angle(res)
    return res


def kick(bodies_state, dv, particle=2):
    """Return a new body list with the particle's velocity kicked by dv=(dvx,dvy).
    `bodies_state` is a list of dicts (e.g. from split_run's checkpoint)."""
    out = [dict(b) for b in bodies_state]
    v = list(out[particle]["vel"])
    out[particle] = dict(out[particle],
                         vel=[v[0] + dv[0], v[1] + dv[1]])
    return out


def state_to_bodies(res, mu=MU):
    """Reconstruct a body list from the final state of a run (for kicks)."""
    n = len(res["masses"])
    d = res["dim"]
    yf = res["yf"]
    pos = yf[: n * d].reshape(n, d)
    vel = yf[n * d:].reshape(n, d)
    return [{"m": float(res["masses"][i]), "pos": pos[i].tolist(),
             "vel": vel[i].tolist()} for i in range(n)]
