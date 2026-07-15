"""WRITE: set the bit by growing the secondary — the horseshoe pinch.

The blank medium is a HORSESHOE co-orbital: bound to the 1:1 resonance but
holding no bit (it sweeps through both islands; reads 'erased'). The write
pulse is a slow, smooth growth of the secondary's mass (mass transferred from
the primary, total fixed — the circular kinematics stay exact). Growing mu
widens the tadpole islands as ~sqrt(mu); when the separatrix sweeps past the
horseshoe's radial offset, the horseshoe PINCHES and the moonlet is captured
into whichever island it is transiting at that moment. The pulse TIMING —
not its shape — selects the bit.

Why this and not something simpler:
  * an impulsive kick cannot write — one conservative kick can only move the
    state along the energy surface, and (verified) super-threshold kicks
    scatter the body out of the resonance entirely: kicks ERASE;
  * dissipation cannot write — drag destabilizes L4/L5 (verified earlier);
  * a slow parameter change CAN — adiabatic capture is robust precisely
    because it doesn't need aim, only timing (Neishtadt, Henrard). This is
    the mechanism by which a growing Jupiter captured its Trojans.

Canonical schedule (mu: 3e-4 -> memory.MU over 100 orbits, blank horseshoe at
da = 0.030): pulse delays of ~8-12 orbits write '1' (L4), ~28-36 write '0'
(L5), calibrated by demos/write_demo.py's timing scan.
"""

import numpy as np

from . import memory, nbody, rotating

MU0 = 3e-4            # blank-medium mass ratio (horseshoe regime at DA)
DA = 0.030            # blank horseshoe's radial offset
PHI0 = 180.0          # blank starts opposite the secondary (far from pinch)
T_RAMP = 100 * memory.PERIOD   # write-pulse duration (adiabatic: >> libration)
T_HOLD = 50 * memory.PERIOD    # settle-and-read window after the pulse
WRITE_DELAY = {"1": 12 * memory.PERIOD,   # fire while transiting the L4 side
               "0": 30 * memory.PERIOD}   # fire while transiting the L5 side
# Deep-capture delays calibrated by the timing scan in demos/write_demo.py:
# captures clear the L3 separatrix by ~46-49 deg and conjunction by ~25 deg.
# Neighbouring delays sit on capture fringes (slow horseshoes) — the guard
# bands between write windows are real and part of the datasheet.


def blank():
    """The unwritten medium: a horseshoe moonlet at mass ratio MU0."""
    return rotating.circular_coorbital(phi0_deg=PHI0, da=DA)


def write_bit(bit, n_samples=1500, rtol=1e-9, atol=1e-10):
    """Write '1' (L4) or '0' (L5) into a blank cell by pulse timing alone.

    Same blank state, same pulse shape — only the firing time differs.
    Returns the full run; read the result with read()."""
    delay = WRITE_DELAY[str(bit)]
    ramp = rotating.smooth_ramp(MU0, memory.MU, t_ramp=T_RAMP, t0=delay)
    t_end = delay + T_RAMP + T_HOLD
    run = rotating.integrate(blank(), t_end, n_samples=n_samples, mu=ramp,
                             rtol=rtol, atol=atol)
    run["bit"] = str(bit)
    return run


def read(run, window=450):
    """Read the stored bit from the tail of a run: 'L4' -> 1, 'L5' -> 0.
    The window must exceed one wide-tadpole libration cycle (~40 orbits),
    else a slow horseshoe can masquerade as a stored bit."""
    label, center, amp = memory.classify(run["phi"][-window:])
    return {"L4": "1", "L5": "0"}.get(label, "erased"), center, amp


def written_cell_bodies(run):
    """Hand the freshly written cell to the full inertial N-body integrator
    (star, planet, moonlet) for an honest hold test at constant mu."""
    return rotating.to_inertial_bodies(run["yf"], run["t"][-1], mu=memory.MU)
