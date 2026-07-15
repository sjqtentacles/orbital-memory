"""Closed-form checks on the memory cell: the rotating-frame effective
potential, the Lagrange points, and the Jacobi constant.

The Jacobi constant C_J is the one conserved quantity of the circular
restricted three-body problem — the rotating-frame analogue of energy. It is
the rigorous backbone of the memory:

  * it is conserved along any orbit (a strong integrator check),
  * a held bit sits at C_J ≈ C_L4 = 3 − μ(1−μ) (the exact triangular value),
  * writing energy into the cell (a kick) lowers C_J; cross the separatrix and
    the tadpole becomes a horseshoe — the bit erases. So C_J is both a readout
    and the quantity the noise margin is measured in.
"""

import numpy as np
from scipy.optimize import brentq

from .memory import MU


def effective_potential(x, y, z=0.0, mu=MU):
    """Rotating-frame (Roche) potential Ω, with n = 1, barycenter at origin.
    Primary (1−μ) at (−μ, 0, 0), secondary (μ) at (1−μ, 0, 0). The
    centrifugal term is planar (the frame rotates about z); gravity is 3D."""
    r1 = np.sqrt((x + mu) ** 2 + y ** 2 + z ** 2)
    r2 = np.sqrt((x - (1 - mu)) ** 2 + y ** 2 + z ** 2)
    return 0.5 * (x ** 2 + y ** 2) + (1 - mu) / r1 + mu / r2


def _collinear_root(mu, bracket):
    def dOmega_dx(x):
        r1 = abs(x + mu)
        r2 = abs(x - (1 - mu))
        return x - (1 - mu) * (x + mu) / r1 ** 3 - mu * (x - (1 - mu)) / r2 ** 3
    return brentq(dOmega_dx, *bracket, xtol=1e-13)


def lagrange_points(mu=MU):
    """The three co-orbital-relevant equilibria: L3 (collinear, opposite the
    secondary) and the triangular L4 (leading) / L5 (trailing)."""
    return {
        "L3": (_collinear_root(mu, (-1.6, -0.9)), 0.0),
        "L4": (0.5 - mu, +np.sqrt(3) / 2),
        "L5": (0.5 - mu, -np.sqrt(3) / 2),
    }


def C_L4(mu=MU):
    """Exact Jacobi constant at the triangular points: C_L4 = C_L5 = 3 − μ(1−μ)."""
    return 3.0 - mu * (1 - mu)


def jacobi_constant(pos, vel, t, mu=MU):
    """C_J from an INERTIAL state, planar or spatial (broadcast over samples).

    Rotate into the frame co-rotating at n = 1 about z, then
    C_J = 2Ω(x_r, y_r, z) − |v_rot|². For 3D states the potential carries z
    and the speed carries vz (the frame rotation leaves both unchanged), so
    this is the full three-dimensional Jacobi integral.
    """
    pos = np.asarray(pos); vel = np.asarray(vel); t = np.asarray(t)
    x, y = pos[..., 0], pos[..., 1]
    vx, vy = vel[..., 0], vel[..., 1]
    c, s = np.cos(-t), np.sin(-t)
    xr, yr = c * x - s * y, s * x + c * y
    # v_rot = R(-t) v_inertial − Ω × r_rot, with Ω = z_hat (n = 1)
    vxr = c * vx - s * vy + yr
    vyr = s * vx + c * vy - xr
    if pos.shape[-1] == 3:
        z = pos[..., 2]
        v2 = vxr ** 2 + vyr ** 2 + vel[..., 2] ** 2
    else:
        z = 0.0
        v2 = vxr ** 2 + vyr ** 2
    return 2 * effective_potential(xr, yr, z, mu) - v2


def jacobi_of(res, body=2, mu=MU):
    """Jacobi-constant time series for one body of an integrate() result."""
    return jacobi_constant(res["traj"][body], res["vel"][body], res["t"], mu)
