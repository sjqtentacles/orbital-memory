"""COOL: shrink a wide bit into a deep one with station-keeping burns.

This is honest spacecraft station-keeping: the stored body is a spacecraft
with a small propulsion budget, and cooling spends a few metres per second of
delta-v to damp a wide libration into a deep, robust one. A wide bit librates
±70°; cooling brings it below 30°. Getting it right took three wrong schemes,
each instructive enough to record:

  * retro-kicks at max rotating-frame speed — expelled the moonlet through
    the L1 neck within kicks (C_J does not confine radially; an inner
    star-grazing orbit shares the tadpole's C_J exactly);
  * "raise C_J toward the L4 ceiling" — BACKWARDS. The measured geometry:
    deep bits sit at C_J just ABOVE C_L4 (the band's floor), the
    tadpole/horseshoe separatrix sits at C_L3 above it, and wide bits near
    that. Raising C_J widens the bit toward the separatrix — impulsive drag,
    the same poison that makes continuous drag destabilize L4/L5;
  * blind prograde kicks — right C_J direction, wrong phase: they walk C_J
    straight through the deep band and out the bottom (C < C_L4 = ejected).

What works is textbook pendulum damping in the co-orbital's own variables.
The tadpole is a slow pendulum whose coordinate is the resonant angle phi
and whose MOMENTUM is the radial offset da = r - 1. A tangential burn
changes da by ~2*dv, so: fire where |da| is large (mid-swing, the momentum
extreme), directed to shrink it — prograde when inside corotation, retro
when outside — sized ~DAMP*|da|/2 and capped well under the erase
threshold. Three or four such burns take an amplitude-66° bit below 30°.
This is how real co-orbital station-keeping works.

A lesson the failures leave behind: C_J STRATIFIES BUT DOES NOT CLASSIFY.
The full Jacobi constant mixes the slow pendulum's energy with the fast
epicyclic (eccentricity) energy, so a deep-but-slightly-eccentric tadpole
can sit BELOW C_L4 while a scattered, erased orbit sits at nearly the same
C_J. Cooling's correctness certificate is operational (amplitude down,
value preserved, survives the honest engine) — the jacobi_gap it reports
records the direction of travel, nothing more.

A deep bit is also a hardened bit: a cooled tadpole sits farther from every
separatrix, so it takes a larger perturbation (a closer flyby) to erase it.
Cold storage is robust storage.
"""

import numpy as np

from . import memory, rotating, theory

DAMP = 0.35                    # fraction of the radial offset removed per burn
CAP = 0.008                    # hard cap on |dv| (memory.ERASE_KICK = 0.035)
MAX_ROUNDS = 20                # coast-and-burn rounds before giving up
WINDOW = 8 * memory.PERIOD     # coast per round (>= one fast epicycle sweep)
SETTLE = 25 * memory.PERIOD    # post-cooling coast for a clean read
TARGET_AMP = 35.0              # deg — the contract bound for a "deep" bit
PHI_SAFE = (55.0, 135.0)       # fire only far from the planet and from L3
DA_MIN = 0.01                  # below this offset, don't bother burning


def jacobi_rotating(state, mu=memory.MU):
    x, y, vx, vy = state
    return float(2 * theory.effective_potential(x, y, mu=mu)
                 - (vx ** 2 + vy ** 2))


def cool(state, t0=0.0, mu=memory.MU, target_amp=TARGET_AMP, damp=DAMP,
         cap=CAP, max_rounds=MAX_ROUNDS, rtol=1e-9, atol=1e-10):
    """Cool a rotating-frame state (a wide tadpole) at constant mu.

    Coast up to WINDOW, read the libration; if still wider than target_amp,
    burn tangentially at the largest safe radial offset and repeat. Returns a
    stitched run dict: t, xy, phi, kick_times, kick_sizes (all <= cap),
    n_kicks, jacobi_gap (|C_J - C_L4| before and after), and yf/t after a
    SETTLE-orbit read coast. Deterministic."""
    state = np.asarray(state, float).copy()
    t_cur = float(t0)
    gap0 = abs(jacobi_rotating(state, mu) - theory.C_L4(mu))
    ts, xys, phis = [], [], []
    kick_times, kick_sizes = [], []

    for _ in range(max_rounds):
        seg = rotating.integrate(state, t_cur + WINDOW, n_samples=320, mu=mu,
                                 rtol=rtol, atol=atol, t0=t_cur)
        label, _, amp = memory.classify(seg["phi"])
        if amp < target_amp and label in ("L4", "L5"):
            ts.append(seg["t"]); xys.append(seg["xy"]); phis.append(seg["phi"])
            state = seg["yf"].copy(); t_cur = float(seg["t"][-1])
            break
        r = np.hypot(seg["xy"][:, 0], seg["xy"][:, 1])
        da = r - 1.0
        aphi = np.abs(seg["phi"])
        safe = ((aphi > PHI_SAFE[0]) & (aphi < PHI_SAFE[1])
                & (np.abs(da) > DA_MIN))
        if not safe.any():
            ts.append(seg["t"]); xys.append(seg["xy"]); phis.append(seg["phi"])
            state = seg["yf"].copy(); t_cur = float(seg["t"][-1])
            continue
        idx = np.where(safe)[0]
        i = int(idx[np.argmax(np.abs(da[idx]))])
        t_kick = float(seg["t"][i])
        ts.append(seg["t"][:i + 1]); xys.append(seg["xy"][:i + 1])
        phis.append(seg["phi"][:i + 1])
        if t_kick > t_cur + 1e-9:
            seg2 = rotating.integrate(state, t_kick, n_samples=8, mu=mu,
                                      rtol=rtol, atol=atol, t0=t_cur)
            state = seg2["yf"].copy()
        # tangential burn sized to shrink the radial offset (da += 2 dv)
        x, y = state[0], state[1]
        rr = float(np.hypot(x, y))
        that = np.array([-y, x]) / rr
        dv = float(np.clip(-damp * (rr - 1.0) / 2.0, -cap, cap))
        state[2] += dv * that[0]
        state[3] += dv * that[1]
        kick_times.append(t_kick)
        kick_sizes.append(abs(dv))
        t_cur = t_kick

    tail = rotating.integrate(state, t_cur + SETTLE, n_samples=1200, mu=mu,
                              rtol=rtol, atol=atol, t0=t_cur)
    ts.append(tail["t"]); xys.append(tail["xy"]); phis.append(tail["phi"])
    gap1 = abs(jacobi_rotating(tail["yf"], mu) - theory.C_L4(mu))

    return {"t": np.concatenate(ts), "xy": np.vstack(xys),
            "phi": np.concatenate(phis),
            "kick_times": np.array(kick_times),
            "kick_sizes": np.array(kick_sizes),
            "n_kicks": len(kick_times),
            "jacobi_gap": (gap0, gap1),
            "yf": tail["yf"]}
