"""WRITE: set the bit by ORBIT INSERTION — a real, actuatable placement.

The bit is which island the body librates in, so writing is *delivering the
body into the chosen island*. That is exactly how you would station a
spacecraft at L4/L5, and the conservative analog of how nature parks a body
in a Trojan orbit: arrive on a co-orbital transfer and perform one insertion
burn that drops you onto a tadpole.

Concretely: the target is a tadpole turning point around L4 (bit '1') or L5
(bit '0') — a point on the corotation circle, `amplitude_deg` from the exact
triangular point, where the rotating-frame velocity is nearly zero. The body
arrives there on a Hohmann-like transfer whose apoapsis just touches the
co-orbital ring from `approach_da` inside; the insertion burn is the prograde
delta-v that circularizes it at that longitude. After the burn the body is in
a tadpole of the requested amplitude — a stored bit. We report the burn in
metres per second (orbital.units): a bit costs propellant to place, and that
is the honest cost, not a grown planet.

No mass is grown or shrunk anywhere. Rewriting and erasing are done by real
gravitational flybys (orbital/gate.py), not by changing the planet.
"""

import numpy as np

from . import memory, nbody

# Default insertion geometry.
AMP0 = 6.0          # deg — default written libration amplitude (a deep bit)
APPROACH_DA = 0.15  # transfer apoapsis offset: large enough that the arrival is
                    # NOT a bit (erased) — the insertion burn is what captures it


def _target_and_arrival(bit, amplitude_deg, approach_da, mu):
    """The target tadpole turning point (position + corotation velocity) and the
    slower velocity the body has when it COASTS into that point on a co-orbital
    transfer whose apoapsis touches the ring from `approach_da` inside."""
    point = "L4" if str(bit) == "1" else "L5"
    target = memory.make_cell(point, mu=mu, libration_deg=amplitude_deg)
    P = np.array(target[2]["pos"])
    v_tadpole = np.array(target[2]["vel"])        # corotation velocity at P
    r = float(np.linalg.norm(P))
    that = np.array([-P[1], P[0]]) / r            # prograde tangential unit
    a_t = r - approach_da / 2.0                    # transfer semi-major axis
    v_apo = float(np.sqrt(max(2.0 / r - 1.0 / a_t, 0.0)))   # vis-viva at apoapsis
    v_arrival = v_apo * that
    return target, P, v_tadpole, v_arrival


def arrival_state(bit, amplitude_deg=AMP0, approach_da=APPROACH_DA, mu=memory.MU):
    """The body at the target point BEFORE the insertion burn: coasting in on
    the transfer at apoapsis speed. On its own this does NOT hold the bit — it
    is on a transfer, slower than corotation, so it drifts off the island. The
    burn is what writes; integrate this to see the difference."""
    target, P, _, v_arrival = _target_and_arrival(bit, amplitude_deg,
                                                   approach_da, mu)
    moon = {"m": 0.0, "pos": P.tolist(), "vel": v_arrival.tolist()}
    return [dict(target[0]), dict(target[1]), moon]


def insert(bit, amplitude_deg=AMP0, approach_da=APPROACH_DA, mu=memory.MU):
    """Write `bit` by ACTUALLY inserting the body into its island.

    The body coasts to the target point on a co-orbital transfer (arrival_state)
    and the insertion burn dv_vec = v_tadpole − v_arrival is APPLIED to its
    velocity, dropping it onto the tadpole. Returns (cell, dv) where the moonlet
    velocity is literally v_arrival + dv_vec and `dv` is the magnitude of that
    applied burn (nondimensional; convert with units.System.mps). The burn is a
    real velocity change, not an accounting number — arrival_state without it
    does not read as the bit."""
    target, P, v_tadpole, v_arrival = _target_and_arrival(bit, amplitude_deg,
                                                          approach_da, mu)
    dv_vec = v_tadpole - v_arrival                  # the insertion burn
    moon = {"m": 0.0, "pos": P.tolist(),
            "vel": (v_arrival + dv_vec).tolist()}   # burn applied
    cell = [dict(target[0]), dict(target[1]), moon]
    return cell, float(np.linalg.norm(dv_vec))


def write_bit(bit, amplitude_deg=AMP0, periods=45, n_samples=None,
              mu=memory.MU):
    """Insert `bit` and integrate the resulting cell in the full inertial
    N-body engine. Returns the run with 'phi' (resonant angle), 'bit', 'dv'
    (nondimensional insertion burn), and 'cell' attached. Read it with read()."""
    cell, dv = insert(bit, amplitude_deg=amplitude_deg, mu=mu)
    if n_samples is None:
        n_samples = int(periods * 60)
    run = nbody.integrate(cell, periods * memory.PERIOD, n_samples=n_samples)
    run["phi"] = memory.resonant_angle(run)
    run["bit"] = str(bit)
    run["dv"] = dv
    run["cell"] = cell
    return run


def blank_cell(mu=memory.MU, da=0.08, phi0_deg=90.0):
    """The unwritten medium: a co-orbital body on a ring `da` outside the
    planet's orbit, wide enough to circulate/horseshoe rather than librate — it
    holds no bit and reads 'erased'. This is the arrival state that insertion
    acts on, and the state an erasing flyby leaves behind."""
    star, planet = memory.primaries(mu)
    a = 1.0 + da
    phi = np.radians(phi0_deg)
    pos = [a * np.cos(phi), a * np.sin(phi)]
    v_circ = np.sqrt(1.0 / a)
    vel = [-v_circ * np.sin(phi), v_circ * np.cos(phi)]
    return [star, planet, {"m": 0.0, "pos": pos, "vel": vel}]


def read(run, window_orbits=40.0):
    """Read the stored bit from the tail of a run.

    Returns (bit, center_deg, amp_deg) with bit '1' (L4), '0' (L5), or the
    string 'erased'. The window is physical time (orbits), converted to samples
    from the run's own clock, so the verdict is independent of sampling density.
    It must exceed one wide-tadpole libration cycle, else a slow horseshoe can
    masquerade as a stored bit."""
    t = run["t"]
    t_cut = t[-1] - window_orbits * memory.PERIOD
    idx = int(np.searchsorted(t, t_cut))
    label, center, amp = memory.classify(run["phi"][idx:])
    return {"L4": "1", "L5": "0"}.get(label, "erased"), center, amp
