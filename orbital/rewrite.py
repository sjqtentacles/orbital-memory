"""ERASE — and the rewrite wall: this memory is write-once, erase-once.

Erase is the write pulse run backwards: shrink mu from the cell value back
to the blank value with the same smoothstep. The tadpole island narrows as
~sqrt(mu); when its width falls below the stored orbit's extent, the
separatrix RELEASES the moonlet back into a bounded, bit-less horseshoe.
That much works, and is tested.

THE REWRITE WALL (a measured negative result, the honest headline here):

  * Release inflates the medium. Crossing the separatrix jumps the adiabatic
    invariant: the blank that was captured at radial offset da = 0.030 comes
    back at da ~ 0.071 (measured; pinned by a test).
  * Recapture has a ceiling. A growth pulse can only re-pinch a horseshoe
    whose offset is below the point where the tadpole band outruns the Hill
    radius: W/r_H = 2.35 mu^(1/6) < 1 throughout the reachable mu range, so
    offsets much beyond ~0.045 are unrecapturable by ANY mu pulse — the
    would-be pinch happens inside Hill-scattering territory. All scanned
    second-write delays read 'erased' (pinned by a test).
  * Every conditioning channel is closed. Station-keeping burns that cool a
    tadpole PUMP a horseshoe (burning toward da = 0 pushes it into the
    sticky chaotic layer around the L3 separatrix); continuous drag
    destabilizes L4/L5 outright; and a slower, "gentler" erase is WORSE —
    at 3x the erase duration the moonlet lingers so long in the separatrix
    layer that it is ejected from the co-orbital region entirely (measured:
    restored da 10+, i.e. gone; pinned by a slow test).

So erasing the bit trashes the medium beyond this design's ability to
rewrite it — an erasure-cost statement with a distinctly Landauer flavor.
A rewritable cell needs a release channel that beats the inflation (e.g. a
shepherding body during release, or capture physics beyond the pure mu
pulse). That is the roadmap item, and it is open.

For contrast: a COOLED bit does not release at all (its area is below the
small-mu island's) — cooling write-protects (orbital/cool.py).
"""

import numpy as np

from . import memory, rotating, write

T_ERASE = 100 * memory.PERIOD    # erase-pulse duration; slower is WORSE (ejects)
T_STORE = 30 * memory.PERIOD     # hold at full mu between write and erase
T_SETTLE = 20 * memory.PERIOD    # coast at blank mu after the release
RELEASE_DA = (0.05, 0.09)        # measured restored-horseshoe offset band
RECAPTURE_CEILING = 0.045        # ~max offset a mu pulse can re-pinch


def erase(state, t0, n_samples=1500, rtol=1e-9, atol=1e-10):
    """Shrink mu MU -> MU0 over T_ERASE, then coast T_SETTLE at blank mu.
    Returns the rotating run; its tail classifies 'erased' and stays
    horseshoe-bounded, at an inflated offset (see RELEASE_DA)."""
    ramp = rotating.smooth_ramp(memory.MU, write.MU0, t_ramp=T_ERASE, t0=t0)
    return rotating.integrate(state, t0 + T_ERASE + T_SETTLE,
                              n_samples=n_samples, mu=ramp,
                              rtol=rtol, atol=atol, t0=t0)


def horseshoe_offset(run, window_orbits=15.0):
    """The medium's radial offset max|r - 1| over the run's tail — the probe
    for release inflation (and for ejection: values >> 0.1 mean 'gone')."""
    t = run["t"]
    idx = int(np.searchsorted(t, t[-1] - window_orbits * memory.PERIOD))
    r = np.hypot(run["xy"][idx:, 0], run["xy"][idx:, 1])
    return float(np.max(np.abs(r - 1.0)))


def l3_crossings(run):
    """Times where the resonant angle sweeps through ±180 (the L3 line) —
    the horseshoe's phase marker."""
    phi = run["phi"]
    flips = np.where(np.abs(np.diff(phi)) > 180.0)[0]
    return run["t"][flips]


def cycle_markers(first, t_erase_duration=T_ERASE):
    """Absolute times of the frozen write -> store -> erase -> settle schedule."""
    t_w1 = write.WRITE_DELAY[str(first)]
    t_erase = t_w1 + write.T_RAMP + T_STORE
    t_settled = t_erase + t_erase_duration + T_SETTLE
    return {"write1": t_w1, "erase": t_erase, "settled": t_settled}


def erased_prefix(first="1", n_samples=2500, rtol=1e-9, atol=1e-10,
                  t_erase_duration=T_ERASE):
    """One continuous run of blank -> write `first` -> store -> erase ->
    settle, via mu_schedule. The state every recapture experiment starts
    from (the release phase is schedule-coupled, so scan and experiment
    must share this prefix bit-for-bit)."""
    m = cycle_markers(first, t_erase_duration)
    sched = rotating.mu_schedule(write.MU0, [
        (m["write1"], write.T_RAMP, memory.MU),
        (m["erase"], t_erase_duration, write.MU0),
    ])
    run = rotating.integrate(write.blank(), m["settled"],
                             n_samples=n_samples, mu=sched,
                             rtol=rtol, atol=atol)
    run["markers"] = m
    return run


def scan_rewrite_delays(delays, first="1", n_samples=2000, rtol=1e-9,
                        atol=1e-10):
    """The recapture experiment: from the settled post-erase state, fire a
    standard growth pulse at each candidate delay and read the result.
    At this design point every delay reads 'erased' — the apparatus that
    measured the rewrite wall, kept as the test harness that pins it."""
    prefix = erased_prefix(first, rtol=rtol, atol=atol)
    e_state, e_t = prefix["yf"], float(prefix["t"][-1])
    out = []
    for delay in delays:
        ramp = rotating.smooth_ramp(write.MU0, memory.MU,
                                    t_ramp=write.T_RAMP, t0=e_t + delay)
        run = rotating.integrate(e_state, e_t + delay + write.T_RAMP
                                 + write.T_HOLD, n_samples=n_samples,
                                 mu=ramp, rtol=rtol, atol=atol, t0=e_t)
        bit, center, amp = write.read(run)
        out.append({"delay": delay, "bit": bit, "center": center,
                    "amp": amp})
    return out
