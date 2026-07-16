"""How long does a bit last? — quantified retention, honestly.

A stored bit is NOT permanent. Deep, low-amplitude tadpoles are Nekhoroshev-
stable — exponentially long-lived, effectively censored survivors of any
feasible integration. But as the libration amplitude approaches the L3
separatrix (orbital/hamiltonian.separatrix_amplitude, ~78 deg), the escape time
drops steeply, and a real perturber (Saturn) erodes bits by pumping their
amplitude toward that edge. Real Jupiter Trojans are metastable on 10 kyr–100
Myr timescales (Greenstreet et al., ApJL 2024), not eternal.

This module measures the feasible part — the escape band near the separatrix —
by direct integration, and is explicit about the deep regime being a censored
survivor whose far-longer lifetime is a theory extrapolation, not a measurement.
Pure scipy; no symplectic integrator, so we never claim a Myr number we did not
integrate.
"""

import numpy as np

from . import memory, nbody

CENSORED = None      # escape_time returns this when the bit survives the horizon


def escape_time(cell, max_orbits=600, chunk_orbits=25, n_per_orbit=40,
                rtol=1e-9, atol=1e-10):
    """Integrate a cell in chunks until the bit reads 'erased' (escaped the
    island) or the horizon is reached. Returns the escape time in orbits, or
    CENSORED (None) if it survives — a lower bound, not eternity. The clock is
    continued across chunks (t0) so the readout frame stays aligned."""
    bodies = [dict(b) for b in cell]
    done = 0
    while done < max_orbits:
        nxt = min(done + chunk_orbits, max_orbits)
        res = nbody.integrate(bodies, nxt * memory.PERIOD,
                              n_samples=int((nxt - done) * n_per_orbit),
                              t0=done * memory.PERIOD, rtol=rtol, atol=atol)
        if memory.classify(memory.resonant_angle(res))[0] == "erased":
            return float(nxt)
        bodies = memory.state_to_bodies(res)
        done = nxt
    return CENSORED


def escape_scan(amplitudes_deg, mu=memory.MU, max_orbits=600, **kw):
    """Escape time vs seeded libration amplitude — the retention datasheet.
    Returns a list of {'amp', 'escape_orbits', 'censored'} dicts. Deep bits come
    back censored (survivors); near-separatrix bits return a finite escape time."""
    out = []
    for amp in amplitudes_deg:
        e = escape_time(memory.make_cell("L4", mu=mu, libration_deg=amp),
                        max_orbits=max_orbits, **kw)
        out.append({"amp": float(amp), "escape_orbits": e,
                    "censored": e is CENSORED})
    return out


def is_censored_survivor(amp_deg, mu=memory.MU, max_orbits=400, **kw):
    """True if a bit at this amplitude survives the horizon (a deep, long-lived
    bit whose true lifetime is beyond what we integrate)."""
    return escape_time(memory.make_cell("L4", mu=mu, libration_deg=amp_deg),
                       max_orbits=max_orbits, **kw) is CENSORED
