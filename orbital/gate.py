"""GATE: a flyby acts on stored memory — the two projects meet.

slingshot-computing does logic with flybys; orbital-memory stores bits in
Lagrange islands. This module points the first at the second: a massive
'bullet' on a fast hyperbolic pass is AIMED (1-D shooting on its launch
direction, root-found on the full 4-body simulation — the slingshot
compiler's trick reused) to shave past the stored moonlet.

    bullet present  -> the moonlet is scattered off its island: bit ERASED
    bullet absent   -> bit unaffected
    bullet grazing  -> (MISS_SAFE, 25x farther) bit unaffected — locality,
                       the datum a dual-rail write head would build on

Aiming lesson (measured): what erases is the impulse's TANGENTIAL share —
the guiding center, which stores the bit, only moves with tangential kicks.
A radial ring-crossing at miss 0.004 delivers a mostly-radial impulse that
pumps a huge epicycle yet leaves the bit readable (amp ~77 deg, still L4);
oblique approaches are gentler still (lower relative speed = adiabatic tug).
The pass must be DEEP: at miss ~0.002 the tangential share finally beats
memory.ERASE_KICK and the moonlet is flung from the resonance.

A conditional erase is logic acting on memory: the bullet's presence is the
input line, the stored bit is the target register.

Two physical facts run the implementation:
  * The bullet must be MASSIVE (a massless body exerts no force), and a
    massive intruder drags the system barycenter — so the readout uses the
    COM-corrected resonant angle (longitudes about the star-planet mass
    center), not the origin-anchored one.
  * Aiming needs the full problem: the bullet is deflected by star and
    planet en route, so its launch direction is root-found on simulated
    closest approach (signed, so brentq can cross zero), exactly as
    slingshot-computing calibrates its gates.
"""

import numpy as np
from scipy.optimize import brentq

from . import memory, nbody

BULLET_M = 2e-4          # lighter than the secondary (3e-3), heavy enough to kick
BULLET_R0 = 4.0          # launch radius
BULLET_V = 1.5           # launch speed: v^2/2 - 1/r0 > 0 -> hyperbolic, exits
MISS_ERASE = 0.002       # aimed closest approach for an erasing pass
MISS_SAFE = 0.05         # a graze at this distance must NOT erase
T_AIM = 6.0              # probe horizon for aiming (encounter at t ~ 2)
T_VERDICT = 45 * memory.PERIOD   # post-shot horizon for the readout


def resonant_angle_com(res, particle=2, primaries=(0, 1)):
    """Resonant angle measured about the INSTANTANEOUS star-planet mass
    center. The massive bullet drags the system barycenter off the origin;
    longitudes about the origin would smear the readout."""
    m = res["masses"]
    ms, mp = m[primaries[0]], m[primaries[1]]
    com = ((ms * res["traj"][primaries[0]] + mp * res["traj"][primaries[1]])
           / (ms + mp))
    d_p = res["traj"][particle] - com
    d_pl = res["traj"][primaries[1]] - com
    lam_p = np.arctan2(d_p[:, 1], d_p[:, 0])
    lam_pl = np.arctan2(d_pl[:, 1], d_pl[:, 0])
    return np.degrees(nbody.wrap_pi(lam_p - lam_pl))


def predict_target(cell, t_enc):
    """Moonlet's inertial position at the encounter time, from a bullet-free
    run of the cell — the straight-line aiming guess."""
    res = nbody.integrate([dict(b) for b in cell], t_enc, n_samples=200)
    return res["traj"][2, -1]


def bullet_body(psi, launch_at):
    """Bullet launched from `launch_at` with speed BULLET_V in direction psi
    (absolute angle, radians)."""
    return {"m": BULLET_M, "pos": list(launch_at),
            "vel": [BULLET_V * np.cos(psi), BULLET_V * np.sin(psi)]}


def _launch_point(cell):
    """Fixed launch position: radially outward of the predicted encounter
    point, so the bullet crosses the co-orbital ring nearly perpendicular
    (fast crossing, minimal dwell near other lanes)."""
    t_enc = (BULLET_R0 - 1.0) / BULLET_V
    target = predict_target(cell, t_enc)
    return BULLET_R0 * target / np.linalg.norm(target), target


def signed_miss(cell, psi, launch_at=None):
    """Closest bullet-moonlet approach within T_AIM, signed by which side
    the bullet passes (z of r_rel x v_rel at closest approach) — a smooth,
    sign-changing function of psi that brentq can root."""
    if launch_at is None:
        launch_at, _ = _launch_point(cell)
    bodies = [dict(b) for b in cell] + [bullet_body(psi, launch_at)]
    try:
        res = nbody.integrate(bodies, T_AIM, n_samples=600, rtol=1e-9,
                              atol=1e-10)
    except RuntimeError:
        return 0.0  # probe collided: treat as a dead-center hit
    rel = res["traj"][3] - res["traj"][2]
    d2 = np.einsum("ij,ij->i", rel, rel)
    i = int(np.argmin(d2))
    # for a near-linear crossing d^2(t) is a parabola: its vertex gives the
    # true minimum even when it falls between samples (sampling otherwise
    # floors the measured miss at ~v_rel * dt / 2)
    if 0 < i < len(d2) - 1:
        y0, y1, y2 = d2[i - 1], d2[i], d2[i + 1]
        denom = y0 - 2 * y1 + y2
        if denom > 0:
            s = 0.5 * (y0 - y2) / denom
            d2_min = y1 - 0.25 * (y0 - y2) * s
        else:
            d2_min = y1
    else:
        d2_min = d2[i]
    vrel = res["vel"][3, i] - res["vel"][2, i]
    side = np.sign(rel[i, 0] * vrel[1] - rel[i, 1] * vrel[0])
    return float(np.sqrt(max(d2_min, 0.0)) * (side if side != 0 else 1.0))


def aim(cell, miss=MISS_ERASE):
    """Find the launch direction whose signed miss equals +miss: scan a
    narrow fan around the straight-line guess for a bracket, then brentq.
    Returns (psi, achieved_unsigned_miss)."""
    launch_at, target = _launch_point(cell)
    guess = float(np.arctan2(target[1] - launch_at[1],
                             target[0] - launch_at[0]))

    def f(psi):
        return signed_miss(cell, psi, launch_at) - miss

    span = 0.06
    grid = np.linspace(guess - span, guess + span, 13)
    vals = [f(p) for p in grid]
    psi = None
    for a, b, fa, fb in zip(grid, grid[1:], vals, vals[1:]):
        if fa == 0.0:
            psi = float(a)
            break
        if fa * fb < 0:
            psi = float(brentq(f, a, b, xtol=1e-7))
            break
    if psi is None:
        raise RuntimeError(f"aim: no bracket around guess; misses={vals}")
    achieved = abs(signed_miss(cell, psi, launch_at))
    return psi, achieved


def fire(cell, psi, t_end=T_VERDICT, n_samples=2000, rtol=1e-9, atol=1e-10):
    """The verdict run: cell + aimed bullet, integrated past the pass and
    long enough to classify the moonlet's fate. Body order: star, planet,
    moonlet, bullet."""
    launch_at, _ = _launch_point(cell)
    bodies = [dict(b) for b in cell] + [bullet_body(psi, launch_at)]
    return nbody.integrate(bodies, t_end, n_samples=n_samples, rtol=rtol,
                           atol=atol)
