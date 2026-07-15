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
APPROACH_DA = 0.05  # co-orbital transfer apoapsis offset (sets the burn size)


def insert(bit, amplitude_deg=AMP0, approach_da=APPROACH_DA, mu=memory.MU):
    """Write `bit` by inserting the body into its island.

    Returns (cell, dv): `cell` is [star, planet, moonlet] with the moonlet
    placed on a tadpole of amplitude ~amplitude_deg around L4 ('1') or L5
    ('0'); `dv` is the NONDIMENSIONAL insertion-burn magnitude (circularizing
    a co-orbital transfer that reaches the ring from `approach_da` inside).
    Convert to m/s with units.System.mps."""
    point = "L4" if str(bit) == "1" else "L5"
    cell = memory.make_cell(point, mu=mu, libration_deg=amplitude_deg)
    r = float(np.linalg.norm(cell[2]["pos"]))     # target radius (~1)
    # circularization burn: v_circ(r) minus the transfer's apoapsis speed, for a
    # transfer with apoapsis r and periapsis r - approach_da (G*M_tot = 1).
    v_circ = np.sqrt(1.0 / r)
    a_t = r - approach_da / 2.0
    v_apo = np.sqrt(max(2.0 / r - 1.0 / a_t, 0.0))
    dv = float(v_circ - v_apo)
    return cell, dv


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
