"""The averaged 1-DOF co-orbital Hamiltonian — the analytic backbone.

Averaging the circular restricted problem over the fast orbital period reduces
the tadpole to a single degree of freedom: the guiding center moving in the
resonant angle phi = lambda - lambda_planet, in an effective potential well
(standard co-orbital averaging — Nesvorny et al. 2002; Murray & Dermott sec 3.13).

The whole model rests on one shape:

    f(phi) = cos(phi) - 1 / (2 sin(phi/2))          (direct + indirect terms)

with guiding-center equation of motion  phi'' = 3 mu f'(phi)  (n = 1 units).
Its equilibria are exactly L4/L5 (phi = +-60 deg, stable) and L3 (phi = 180 deg,
the separatrix saddle); f''(60 deg) = -9/4 gives the small-amplitude frequency
omega^2 = (27/4) mu (matching memory.libration_period), and the L4->L3 barrier
is exactly 3 mu. Everything below is quadrature and root-finding on f.

This module turns the memory's measured constants into closed forms:
  * libration_period(amp)   — the nonlinear, amplitude-dependent period
  * separatrix_amplitude()  — the widest tadpole before it becomes a horseshoe
  * erase_margin(amp)       — the tangential dv that erases a bit (derives ERASE_KICK)

The full n-body sim is then the *verifier* of these forms, not the design tool.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from . import memory

PHI_L4 = 60.0       # deg — leading triangular point (bit 1)
PHI_L3 = 180.0      # deg — the tadpole/horseshoe separatrix saddle


def f(phi_rad):
    """The averaged co-orbital shape (direct 1/distance + indirect cos term)."""
    return np.cos(phi_rad) - 1.0 / (2.0 * np.sin(phi_rad / 2.0))


def df(phi_rad):
    """d f / d phi."""
    return -np.sin(phi_rad) + 0.25 * np.cos(phi_rad / 2.0) / np.sin(phi_rad / 2.0) ** 2


def potential(phi_deg, mu=memory.MU):
    """The averaged potential well V(phi): 0 at L4 (60 deg), rising to the
    barrier 3*mu at L3 (180 deg). The tadpole librates inside this well."""
    p = np.radians(phi_deg)
    return 3.0 * mu * (f(np.radians(PHI_L4)) - f(p))


def barrier_height(mu=memory.MU):
    """The L4 -> L3 potential barrier: exactly 3*mu."""
    return 3.0 * mu


def _separatrix_conjunction_deg():
    """Conjunction-side turning point of the widest tadpole (E = barrier),
    i.e. the phi < 60 deg where f(phi) = f(180 deg). mu-independent."""
    f180 = f(np.radians(PHI_L3))
    return brentq(lambda d: f(np.radians(d)) - f180, 1.0, PHI_L4 - 1e-9)


def separatrix_amplitude():
    """The maximum tadpole libration amplitude (classify-style: half the
    peak-to-peak excursion) before the orbit crosses L3 into a horseshoe.
    mu-independent (mu only scales the energy, not the angular extent) — ~78 deg."""
    phi_c = _separatrix_conjunction_deg()
    return 0.5 * (PHI_L3 - phi_c)


def _turning_points_deg(amp_deg, mu):
    """Given a release from rest at phi = 60 + amp_deg (the L3 side, matching
    memory.make_cell's seed), the two libration turning points in degrees."""
    phi_hi = PHI_L4 + amp_deg
    E = potential(phi_hi, mu)
    # conjunction-side turning point: V(phi_lo) = E, phi_lo in (phi_c, 60)
    phi_c = _separatrix_conjunction_deg()
    phi_lo = brentq(lambda d: potential(d, mu) - E, phi_c + 1e-6, PHI_L4 - 1e-9)
    return phi_lo, phi_hi


def libration_period(amp_deg, mu=memory.MU):
    """Nonlinear libration period for a tadpole released from rest at
    phi = 60 + amp_deg (matching memory.make_cell(libration_deg=amp_deg)), by
    quadrature T = 2 * integral d phi / sqrt(2 (E - V(phi))) between turning
    points. In nondimensional time (planet period 2*pi). Reduces to
    memory.libration_period() = 2*pi/sqrt(27/4 mu) as amp -> 0."""
    if amp_deg < 1e-6:
        return memory.libration_period(mu)
    phi_lo, phi_hi = _turning_points_deg(amp_deg, mu)
    E = potential(phi_hi, mu)
    lo, hi = np.radians(phi_lo), np.radians(phi_hi)

    def integrand(p):
        val = 2.0 * (E - 3.0 * mu * (f(np.radians(PHI_L4)) - f(p)))
        return 1.0 / np.sqrt(val) if val > 0 else 0.0

    # endpoints are integrable square-root singularities -> use 'weight' via
    # substitution is overkill; quad handles them with the points hint.
    half, _ = quad(integrand, lo, hi, limit=200, points=None)
    return 2.0 * half


def erase_margin(amp_deg=2.0, mu=memory.MU):
    """The tangential burn dv (as a fraction of orbital speed) that lifts a bit
    at libration amplitude amp_deg over the L3 barrier — the analytic noise
    margin, to be compared with the measured memory.ERASE_KICK.

    Mapping: a tangential kick dv changes the semimajor axis (da = 2 dv, vis-viva
    at unit radius/speed), hence the guiding-center rate phi' by -3 dv (n=1). It
    adds kinetic energy 1/2 (3 dv)^2 = 4.5 dv^2 to the 1-DOF system. Erasing needs
    that to carry the bit from its current energy over the barrier:
        4.5 dv^2 >= barrier - E(amp)
    NOTE: this is the L3-crossing margin. The averaged model breaks down near
    conjunction (phi -> 0, close planetary encounter), where the real sim can
    scatter a bit out somewhat more easily, so the measured ERASE_KICK runs a
    bit below this leading-order estimate."""
    E = potential(PHI_L4 + amp_deg, mu)
    reach = max(barrier_height(mu) - E, 0.0)
    return float(np.sqrt(reach / 4.5))
