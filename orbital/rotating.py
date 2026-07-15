"""The restricted problem in its natural coordinates: the co-rotating frame.

Here the primaries sit still on the x-axis, L4/L5 are literal fixed points,
the resonant angle is atan2(y, x), and the memory's phase-space anatomy
(tadpole islands, horseshoe band, circulation) is visible to the naked eye.
Integrating 4 ODEs instead of 12 makes parameter sweeps ~2 orders of
magnitude cheaper than the full inertial N-body — and, crucially, the mass
ratio may be TIME-DEPENDENT, which is the write mechanism: growing the
secondary adiabatically captures a circulating body into a chosen island,
the same way a growing Jupiter captured its Trojans.

Equations of motion (n = 1, total mass 1, separation 1, barycenter origin):
    x'' - 2 y' = dOmega/dx,   y'' + 2 x' = dOmega/dy
with the primaries pinned at (-mu, 0) and (1-mu, 0). For time-dependent
mu(t) the mass is transferred star -> planet (total fixed at 1), so the
circular two-body kinematics stay exact throughout the ramp.
"""

import numpy as np
from scipy.integrate import solve_ivp

from .memory import MU


def _accel(x, y, vx, vy, mu):
    r1 = np.hypot(x + mu, y)
    r2 = np.hypot(x - (1 - mu), y)
    ax = (x + 2 * vy
          - (1 - mu) * (x + mu) / r1 ** 3
          - mu * (x - (1 - mu)) / r2 ** 3)
    ay = (y - 2 * vx
          - (1 - mu) * y / r1 ** 3
          - mu * y / r2 ** 3)
    return ax, ay


def _rhs(t, s, mu_of_t):
    x, y, vx, vy = s
    ax, ay = _accel(x, y, vx, vy, mu_of_t(t))
    return [vx, vy, ax, ay]


def smooth_ramp(mu0, mu1, t_ramp, t0=0.0):
    """A C^1 monotone mass-growth schedule: mu0 before t0, smoothstep over
    [t0, t0 + t_ramp], mu1 after. The write pulse."""
    def mu_of_t(t):
        s = np.clip((t - t0) / t_ramp, 0.0, 1.0)
        return mu0 + (mu1 - mu0) * (3 * s * s - 2 * s ** 3)
    return mu_of_t


def integrate(state0, t_end, n_samples=2000, mu=MU, rtol=1e-11, atol=1e-12,
              t0=0.0):
    """Integrate a rotating-frame state [x, y, vx, vy].

    mu may be a scalar (constant) or a callable mu(t). Returns t, xy, v,
    phi (the resonant angle, degrees), mu(t) samples, and the Jacobi
    'constant' time series (exactly conserved only for constant mu — during
    a write pulse it MOVES, and that motion is the write energy).
    """
    mu_of_t = mu if callable(mu) else (lambda t, _m=mu: _m)
    times = np.linspace(t0, t_end, n_samples)
    sol = solve_ivp(_rhs, (t0, t_end), np.asarray(state0, float),
                    args=(mu_of_t,), method="DOP853", t_eval=times,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    x, y, vx, vy = sol.y
    mus = np.array([mu_of_t(t) for t in times])
    r1 = np.hypot(x + mus, y)
    r2 = np.hypot(x - (1 - mus), y)
    omega = 0.5 * (x ** 2 + y ** 2) + (1 - mus) / r1 + mus / r2
    jacobi = 2 * omega - (vx ** 2 + vy ** 2)
    return {"t": times, "xy": np.column_stack([x, y]),
            "v": np.column_stack([vx, vy]),
            "phi": np.degrees(np.arctan2(y, x)),
            "mu": mus, "jacobi": jacobi,
            "yf": np.array([x[-1], y[-1], vx[-1], vy[-1]])}


# ---------------------------------------------------------------------------
# state constructors
# ---------------------------------------------------------------------------

def from_inertial(pos, vel, t):
    """Convert an inertial planar state to the rotating frame at time t."""
    c, s = np.cos(-t), np.sin(-t)
    x = c * pos[0] - s * pos[1]
    y = s * pos[0] + c * pos[1]
    vx = c * vel[0] - s * vel[1] + y
    vy = s * vel[0] + c * vel[1] - x
    return np.array([x, y, vx, vy])


def to_rotating_frame(res, body):
    """One body's trajectory from an inertial nbody.integrate() run, rotated
    into the co-rotating frame (n = 1). Returns (n_samples, 2) positions; a
    3D run's z column is dropped (the frame rotates about z). The single
    shared implementation of the transform every demo, test and figure needs."""
    t = res["t"]
    xy = res["traj"][body]
    c, s = np.cos(-t), np.sin(-t)
    return np.column_stack([c * xy[:, 0] - s * xy[:, 1],
                            s * xy[:, 0] + c * xy[:, 1]])


def to_inertial_bodies(state, t, mu=MU):
    """Rotating-frame particle state at time t -> full inertial 3-body list
    (star, planet, particle), e.g. to hand a written bit to nbody.integrate."""
    x, y, vx, vy = state
    c, s = np.cos(t), np.sin(t)

    def rot(px, py):
        return [c * px - s * py, s * px + c * py]

    star_p = rot(-mu, 0.0)
    planet_p = rot(1 - mu, 0.0)
    part_p = rot(x, y)
    # inertial velocity = R(t) (v_rot + omega x r_rot), omega = z_hat
    part_v = rot(vx - y, vy + x)
    star_v = [-star_p[1], star_p[0]]
    planet_v = [-planet_p[1], planet_p[0]]
    return [{"m": 1 - mu, "pos": star_p, "vel": star_v},
            {"m": mu, "pos": planet_p, "vel": planet_v},
            {"m": 0.0, "pos": part_p, "vel": part_v}]


def lagrange_tadpole(point="L4", libration_deg=6.0, mu=MU):
    """Rotating-frame state seeded on the corotation circle, `libration_deg`
    ahead of the chosen triangular point. (Angles are measured from exactly
    ±60°, so this differs from memory.make_cell's rotated-L4-vector seed by
    the tiny atan2(√3/2, ½−μ) − 60° ≈ 0.17° offset — both are valid tadpole
    seeds; they are not bit-for-bit identical.)"""
    sy = 1.0 if point == "L4" else -1.0
    ang = np.radians(60.0 * sy + libration_deg)
    r = np.hypot(0.5 - mu, np.sqrt(3) / 2)
    return np.array([r * np.cos(ang), r * np.sin(ang), 0.0, 0.0])


def circular_coorbital(phi0_deg, da=0.03):
    """A blank medium: a body on an inertial circular orbit of radius 1 + da
    at longitude phi0 relative to the secondary. In the rotating frame it
    slowly circulates (drift rate ~ -1.5 da per unit time), waiting to be
    written."""
    a = 1.0 + da
    phi = np.radians(phi0_deg)
    pos = np.array([a * np.cos(phi), a * np.sin(phi)])
    v_circ = 1.0 / np.sqrt(a)           # inertial circular speed, G M_tot = 1
    that = np.array([-np.sin(phi), np.cos(phi)])
    v_rot = (v_circ - a) * that          # subtract corotation
    return np.array([pos[0], pos[1], v_rot[0], v_rot[1]])
