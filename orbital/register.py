"""A multi-bit register: several co-orbital cells at spaced radii.

A single planet has only one L4/L5, so N independent bits need N distinct hosts.
This is real hardware: Saturn's moons Telesto/Calypso ride Tethys, Helene/
Polydeuces ride Dione — different secondaries at different radii, each with its
own co-orbital Trojan. Here a central star carries N light secondaries on spaced
circular orbits; each secondary's Trojan sits at its L4 (bit 1) or L5 (bit 0).

Body order: [star, planet_0, trojan_0, planet_1, trojan_1, ...]. The whole thing
runs in the honest inertial `nbody` engine (a nibble is ~9 bodies; the swarm's
single-frame vectorized field does NOT apply — the hosts move relative to each
other, so there is no shared rotating frame). Register capacity is bounded by
mutual (Hill) stability of the secondaries: a nibble (4 bits) is reliable at
spacing >= 1.8; push further only as spacing/mass allow.
"""

import numpy as np

from . import memory, nbody

SPACING = 1.8               # geometric radius ratio between adjacent hosts
MU_EACH = memory.MU         # secondary mass fraction per host


def radii(n_bits, spacing=SPACING):
    """Host orbital radii: r_i = spacing**i (innermost host at radius 1)."""
    return [spacing ** i for i in range(n_bits)]


def make_register(bits, spacing=SPACING, mu_each=MU_EACH):
    """Build the register body list for `bits` (e.g. "1010"): star + one
    (secondary, Trojan) pair per bit. Total momentum is zeroed via the star and
    the barycenter placed at the origin, so the origin-anchored resonant-angle
    readout stays clean."""
    bits = str(bits)
    rs = radii(len(bits), spacing)
    bodies = [{"m": 1.0, "pos": [0.0, 0.0], "vel": [0.0, 0.0]}]
    for r, b in zip(rs, bits):
        v = np.sqrt(1.0 / r)                       # circular speed about the star
        bodies.append({"m": mu_each, "pos": [r, 0.0], "vel": [0.0, v]})
        ang = np.radians(60.0 if b == "1" else -60.0)
        bodies.append({"m": 0.0, "pos": [r * np.cos(ang), r * np.sin(ang)],
                       "vel": [-v * np.sin(ang), v * np.cos(ang)]})
    m = np.array([b["m"] for b in bodies])
    P = sum(b["m"] * np.array(b["vel"]) for b in bodies)
    bodies[0]["vel"] = (-P / m[0]).tolist()        # zero total momentum
    com = sum(b["m"] * np.array(b["pos"]) for b in bodies) / m.sum()
    for b in bodies:
        b["pos"] = (np.array(b["pos"]) - com).tolist()
    return bodies


def _cell_libration_period(i, spacing, mu_each):
    orbit = 2 * np.pi * radii(i + 1, spacing)[i] ** 1.5
    return orbit / np.sqrt(27.0 / 4.0 * mu_each)


def hold_register(bits, spacing=SPACING, mu_each=MU_EACH, n_librations=8,
                  n_per_orbit=60):
    """Integrate a written register long enough for the slowest (outer) cell to
    librate `n_librations` times. Returns the run with register metadata
    attached (`n_bits`, `spacing`, `mu_each`) for read_register."""
    bits = str(bits)
    rs = radii(len(bits), spacing)
    t_end = n_librations * _cell_libration_period(len(bits) - 1, spacing, mu_each)
    n_samples = max(4000, int(t_end / (2 * np.pi) * n_per_orbit))
    run = nbody.integrate(make_register(bits, spacing, mu_each), t_end,
                          n_samples=n_samples)
    run["n_bits"] = len(bits)
    run["spacing"] = spacing
    run["mu_each"] = mu_each
    return run


def read_register(run, window_librations=4):
    """Read the register back as a bit-string. Each cell is classified over the
    last `window_librations` of ITS OWN libration period (outer cells librate
    slower), relative to its own host. A cell that has escaped reads '?'."""
    n, spacing, mu = run["n_bits"], run["spacing"], run["mu_each"]
    t = run["t"]
    out = ""
    for i in range(n):
        lib = _cell_libration_period(i, spacing, mu)
        idx = int(np.searchsorted(t, max(t[-1] - window_librations * lib, 0.0)))
        phi = memory.resonant_angle(run, particle=2 + 2 * i, planet=1 + 2 * i)
        label = memory.classify(phi[idx:])[0]
        out += {"L4": "1", "L5": "0"}.get(label, "?")
    return out


def cell_amplitude(run, i, window_librations=4):
    """Libration amplitude (deg) of cell i over its own read window — the probe
    for crosstalk (a coupled cell librates wider than an isolated one)."""
    spacing, mu, t = run["spacing"], run["mu_each"], run["t"]
    lib = _cell_libration_period(i, spacing, mu)
    idx = int(np.searchsorted(t, max(t[-1] - window_librations * lib, 0.0)))
    phi = memory.resonant_angle(run, particle=2 + 2 * i, planet=1 + 2 * i)
    return memory.classify(phi[idx:])[2]


def crosstalk(bits, spacing=SPACING, mu_each=MU_EACH):
    """Max libration amplitude induced across the cells by their neighbours.
    Each cell is seeded exactly at its triangular point (amplitude ~0 in
    isolation), so any libration it acquires is pure coupling. Falls with
    spacing — the crosstalk-vs-spacing datasheet."""
    run = hold_register(bits, spacing, mu_each)
    return max(cell_amplitude(run, i) for i in range(len(str(bits))))
