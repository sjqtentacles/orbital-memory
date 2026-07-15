"""A gravitational memory bit: which Lagrange island a body librates in.

Circular restricted three-body problem (G=1, total mass 1, separation 1, mean
motion n=1, period 2*pi). A heavy primary and a secondary orbit their
barycenter; a massless test body sits in a tadpole orbit around the leading
(L4) or trailing (L5) triangular Lagrange point. Those two islands are
separated by a separatrix (the horseshoe region around L3), so a small
perturbation can't move the body between them — the island it librates in is a
nonvolatile bit.

This is a moon-scale mechanism, not just a toy: Saturn's moons Telesto and
Calypso sit in the L4/L5 points of Tethys, and Helene/Polydeuces in those of
Dione. Those co-orbital moons are exactly this cell — real hardware that has
held its "bit" for the age of the solar system. The dynamics depend only on
the mass ratio `mu`, so the same code covers star+planet, planet+moon, and
moon+co-orbital-moonlet; only the libration timescale (~1/sqrt(mu)) changes.

Encoding: librating around L4 (+60 deg ahead of the secondary) = '1';
around L5 (-60 deg, trailing) = '0'; a horseshoe past L3 or an escape = erased.
"""

import numpy as np

from . import nbody

MU = 0.003          # secondary mass fraction (< 0.0385 -> L4/L5 linearly stable)
PERIOD = 2 * np.pi  # orbital period of the primary pair

# The measured noise margin: the smallest tangential kick (as a fraction of
# orbital speed) that ejects a deep tadpole across the separatrix. Kicks
# below this leave the bit intact; at or above it, the bit erases.
ERASE_KICK = 0.035

# A few real co-orbital mass ratios (secondary / total), for reference and
# for parameterizing the cell at moon scale. Libration slows as ~1/sqrt(mu).
MU_SUN_JUPITER = 9.54e-4       # the classic Jupiter Trojans
MU_SUN_EARTH = 3.00e-6         # 2010 TK7, Earth's (large-amplitude) Trojan
MU_SATURN_TETHYS = 1.09e-6     # Telesto & Calypso ride Tethys's L4/L5
MU_SATURN_DIONE = 1.85e-6      # Helene & Polydeuces ride Dione's L4/L5


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


def primaries_eccentric(mu=MU, e=0.0489, a=1.0):
    """Star + planet on a Kepler orbit of relative semimajor axis `a` and
    eccentricity `e`, started at periapsis (barycenter at the origin,
    G*M_tot = 1). Real Jupiter has e ~ 0.0489. Reduces to primaries() at e=0."""
    r_peri = a * (1 - e)
    v_peri = np.sqrt((1 + e) / (1 - e) / a)     # vis-viva at periapsis, GM=1
    r_rel = np.array([r_peri, 0.0])
    v_rel = np.array([0.0, v_peri])
    star = {"m": 1 - mu, "pos": (-mu * r_rel).tolist(),
            "vel": (-mu * v_rel).tolist()}
    planet = {"m": mu, "pos": ((1 - mu) * r_rel).tolist(),
              "vel": ((1 - mu) * v_rel).tolist()}
    return star, planet


def make_cell_eccentric(state="L4", mu=MU, e=0.0489, libration_deg=6.0, a=1.0):
    """A tadpole bit around L4/L5 of an ECCENTRIC (elliptic restricted) system.

    The equilateral point is a central configuration: L4/L5 relative to the
    barycenter is a fixed linear image of the star->planet vector,
    M = (1/2 - mu) I + (sqrt3/2) J (J a +/-90-deg rotation for L4/L5). Seeding
    the particle at M r_rel with velocity M v_rel puts it exactly on the
    pulsating-rotating equilateral solution; a small `libration_deg` nudge sets
    it librating. The whole tadpole then breathes at the planet's orbital
    frequency."""
    star, planet = primaries_eccentric(mu, e, a)
    r_rel = np.array(planet["pos"]) - np.array(star["pos"])
    v_rel = np.array(planet["vel"]) - np.array(star["vel"])
    sy = 1.0 if state == "L4" else -1.0
    J = sy * np.array([[0.0, -1.0], [1.0, 0.0]])       # +90 (L4) / -90 (L5)
    M = (0.5 - mu) * np.eye(2) + (np.sqrt(3) / 2) * J
    pos, vel = M @ r_rel, M @ v_rel
    R = _rot(np.radians(libration_deg))                # nudge along the orbit
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


# Real Saturn, in Sun-Jupiter units (total Sun+Jupiter mass = 1, a_Jupiter = 1).
MASS_SATURN = 2.855e-4     # M_Saturn / (M_Sun + M_Jupiter)
A_SATURN = 1.8412          # a_Saturn / a_Jupiter (9.582 AU / 5.204 AU)


def sun_jupiter_saturn_cell(state="L4", libration_deg=6.0,
                            mu=MU_SUN_JUPITER, with_saturn=True,
                            saturn_phase_deg=180.0):
    """A Jupiter Trojan bit, optionally with real Saturn as a 4th body.

    Bodies: [Sun, Jupiter, Trojan(massless), Saturn]. Saturn rides a circular
    orbit at the real spacing A_SATURN with the real mass ratio MASS_SATURN,
    started `saturn_phase_deg` from Jupiter. This is the real mechanism that
    erodes the Trojan swarm; with_saturn=False is the unperturbed control."""
    cell = make_cell(state, mu=mu, libration_deg=libration_deg)
    if not with_saturn:
        return cell
    ang = np.radians(saturn_phase_deg)
    r = A_SATURN
    pos = np.array([r * np.cos(ang), r * np.sin(ang)])
    v = np.sqrt(1.0 / r)                        # circular about the interior mass
    vel = np.array([-v * np.sin(ang), v * np.cos(ang)])
    saturn = {"m": MASS_SATURN, "pos": pos.tolist(), "vel": vel.tolist()}
    return cell + [saturn]


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
    """Read the bit from a resonant-angle time series (degrees, wrapped).

    Returns (label, center_deg, amplitude_deg), decided by which separatrix
    lines the trajectory crosses — the physical definition, valid for
    arbitrarily wide tadpoles:

      * a TADPOLE never reaches conjunction (phi = 0) nor L3 (phi = ±180):
        no sign changes at all -> the bit, L4 (+) or L5 (−) by mean side;
      * a HORSESHOE turns around near the secondary but sweeps through L3:
        sign changes only at ±180 -> 'erased';
      * CIRCULATION sweeps through both -> 'erased'.
    """
    phi = np.asarray(phi, dtype=float)
    if phi.size == 0:
        raise ValueError("classify: empty resonant-angle series")
    center = _circular_mean_deg(phi)
    amp = float(np.max(np.abs(nbody.wrap_pi(np.radians(phi) - np.radians(center))))
                * 180 / np.pi)
    if phi.size < 10:
        return "erased", center, amp  # too short to be evidence of libration
    flips = np.where(np.diff(np.sign(phi)) != 0)[0]
    if len(flips) == 0:
        return ("L4" if center > 0 else "L5"), center, amp
    return "erased", center, amp  # horseshoe or circulation: no stored bit


def hold(state="L4", periods=40, mu=MU, libration_deg=6.0, n_per_period=60):
    """Integrate a freshly written cell and return the run + resonant angle."""
    bodies = make_cell(state, mu, libration_deg)
    t_end = periods * PERIOD
    res = nbody.integrate(bodies, t_end, n_samples=int(periods * n_per_period))
    res["phi"] = resonant_angle(res)
    return res


def kick(bodies_state, dv, particle=2):
    """Return a new body list with the particle's velocity kicked by dv.

    Dimension-safe: dv may be shorter than the velocity (a 2-vector dv applied
    to a 3D cell leaves vz untouched); extra dv components are ignored.
    """
    out = [dict(b) for b in bodies_state]
    v = list(out[particle]["vel"])
    v = [v[i] + (dv[i] if i < len(dv) else 0.0) for i in range(len(v))]
    out[particle] = dict(out[particle], vel=v)
    return out


def kicked_cell(frac, state="L4", libration_deg=2.0, settle_periods=5,
                mu=MU, direction=None):
    """The standard noise-margin experiment, packaged: seed a deep tadpole,
    let it settle `settle_periods`, then kick the moonlet tangentially by
    `frac` of its orbital speed. Returns the kicked body list, ready for
    nbody.integrate. `direction` overrides the unit kick direction."""
    cell = make_cell(state, mu, libration_deg)
    pre = nbody.integrate(cell, settle_periods * PERIOD, n_samples=200)
    bodies = state_to_bodies(pre)
    v = np.array(bodies[2]["vel"])
    speed = np.linalg.norm(v)
    uhat = np.asarray(direction, float) if direction is not None else v / speed
    return kick(bodies, (uhat * speed * frac).tolist())


def libration_period(mu=MU):
    """Analytic tadpole (long-period) libration period around L4/L5 for small
    mu: T = 2*pi / sqrt(27/4 * mu), in the same time units (planet period 2*pi).
    Used to validate the simulation against linearized Trojan theory."""
    return 2 * np.pi / np.sqrt(27.0 / 4.0 * mu)


def measured_libration_period(t, phi):
    """Measure the tadpole libration period from a resonant-angle series phi(t):
    twice the mean spacing between successive crossings of the libration center
    (sub-sample interpolated). This is the simulation's answer to compare
    against libration_period() and against observed Trojan periods once
    converted to years. Requires at least two center crossings."""
    t = np.asarray(t, dtype=float)
    phi = np.asarray(phi, dtype=float)
    s = phi - _circular_mean_deg(phi)
    idx = np.where(np.diff(np.sign(s)) != 0)[0]
    if len(idx) < 2:
        raise ValueError("measured_libration_period: need >= 2 center crossings")
    tc = [t[i] + s[i] / (s[i] - s[i + 1]) * (t[i + 1] - t[i]) for i in idx]
    return 2.0 * float(np.mean(np.diff(tc)))


def state_to_bodies(res, mu=MU):
    """Reconstruct a body list from the final state of a run (for kicks)."""
    n = len(res["masses"])
    d = res["dim"]
    yf = res["yf"]
    pos = yf[: n * d].reshape(n, d)
    vel = yf[n * d:].reshape(n, d)
    return [{"m": float(res["masses"][i]), "pos": pos[i].tolist(),
             "vel": vel[i].tolist()} for i in range(n)]
