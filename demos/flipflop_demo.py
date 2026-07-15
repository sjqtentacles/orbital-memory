"""Demo 01: a gravitational memory bit that holds — and its noise margin.

A test particle librates in a tadpole around the leading (L4) or trailing (L5)
Trojan point of a star+planet system. Which island it's in is a nonvolatile
bit: the two islands are separated by a separatrix, so any perturbation smaller
than a threshold leaves the bit intact, and a larger one erases it.

The protection is TOPOLOGICAL, not dissipative — L4/L5 are Coriolis-stabilized
potential maxima, so adding drag would destabilize them. The bit is held by the
phase-space geometry (an invariant island), the KAM/Nekhoroshev way.

Writes out/flipflop.json (for the viewer) and out/flipflop.png.
Usage: python -m demos.flipflop_demo
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from orbital import memory, nbody

OUT = pathlib.Path(__file__).resolve().parent.parent / "out"


def rotating_frame(res, body):
    """Body trajectory rotated into the frame co-rotating with the planet
    (n = 1), so the star, planet and Lagrange points sit still."""
    t = res["t"]
    xy = res["traj"][body]
    c, s = np.cos(-t), np.sin(-t)
    x = c * xy[:, 0] - s * xy[:, 1]
    y = s * xy[:, 0] + c * xy[:, 1]
    return np.column_stack([x, y])


def noise_margin(periods=50, samples=2500):
    """Kick an L4 bit with growing tangential impulses; find the fraction of
    orbital speed at which the tadpole is destroyed (bit erased)."""
    base = memory.make_cell("L4", libration_deg=2.0)
    pre = nbody.integrate(base, 5 * memory.PERIOD, n_samples=200)
    bodies = memory.state_to_bodies(pre)
    v = np.array(bodies[2]["vel"])
    uhat = v / np.linalg.norm(v)
    speed = np.linalg.norm(v)
    rows = []
    for frac in np.arange(0.005, 0.055, 0.005):
        dv = (uhat * speed * frac).tolist()
        res = nbody.integrate(memory.kick(bodies, dv), periods * memory.PERIOD,
                              n_samples=samples)
        phi = memory.resonant_angle(res)
        label, center, amp = memory.classify(phi)
        rows.append((frac, label, amp))
    return rows, speed


def main():
    OUT.mkdir(exist_ok=True)
    payload = {}

    print("=== HOLD: two nonvolatile states ===")
    holds = {}
    for st in ("L4", "L5"):
        res = memory.hold(st, periods=80, libration_deg=6.0)
        phi = res["phi"]
        label, center, amp = memory.classify(phi)
        holds[st] = res
        print(f"  {st}: librates around {center:+5.1f} deg (amp {amp:.1f}), "
              f"reads {label!r} for 80 orbits, drift {res['energy_drift']:.1e}")
        payload[st] = {
            "t": res["t"].tolist(),
            "phi": np.round(phi, 3).tolist(),
            "rot_particle": np.round(rotating_frame(res, 2), 4).tolist(),
            "center": round(center, 2), "amp": round(amp, 2), "read": label,
        }
    assert memory.classify(holds["L4"]["phi"])[0] == "L4"
    assert memory.classify(holds["L5"]["phi"])[0] == "L5"
    print("  -> the bit holds. (Jupiter's real Trojans have held theirs for 4.5 Gyr.)")

    print("\n=== NOISE MARGIN: how hard a kick erases the bit ===")
    rows, speed = noise_margin()
    threshold = next((f for f, lab, _ in rows if lab == "erased"), None)
    for frac, label, amp in rows:
        bar = "kept " if label != "erased" else "ERASED"
        print(f"  kick {frac*100:4.1f}% of orbital speed -> {bar} "
              f"(libration amp {amp:5.1f} deg, reads {label})")
    print(f"  -> separatrix threshold ~{threshold*100:.1f}% of orbital speed: "
          f"a real noise margin, set by phase-space geometry, not dissipation.")
    payload["noise_margin"] = {
        "speed": speed, "threshold_frac": threshold,
        "rows": [{"frac": f, "label": l, "amp": a} for f, l, a in rows],
    }

    # geometry for the viewer
    star, planet = memory.primaries()
    l4, _ = memory.lagrange_state(memory.MU, "L4")
    l5, _ = memory.lagrange_state(memory.MU, "L5")
    payload["geometry"] = {
        "mu": memory.MU, "period": memory.PERIOD,
        "star": star["pos"], "planet": planet["pos"],
        "L4": l4.tolist(), "L5": l5.tolist(),
    }
    (OUT / "flipflop.json").write_text(json.dumps(payload))
    print(f"\nwrote {OUT / 'flipflop.json'}")

    # --- figure: rotating-frame tadpoles + phi(t) traces ---
    fig = plt.figure(figsize=(12, 5.2), facecolor="#0a0e17")
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.22)

    ax = fig.add_subplot(gs[0]); ax.set_facecolor("#0a0e17")
    ax.set_aspect("equal"); ax.set_title("rotating frame — the bit is *where* it orbits",
                                         color="w", fontsize=12)
    ax.plot(0, 0, "+", color="#445", ms=10)  # barycenter
    ax.plot(*star["pos"], "o", color="#ffd166", ms=16, label="star")
    ax.plot(*planet["pos"], "o", color="#9fb0dd", ms=9, label="planet")
    for pt, lab, col in ((l4, "L4 = 1", "#54d1ff"), (l5, "L5 = 0", "#ff8fa3")):
        ax.plot(*pt, "x", color=col, ms=9)
        ax.annotate(lab, pt, textcoords="offset points", xytext=(8, 6),
                    color=col, fontsize=11)
    for st, col in (("L4", "#54d1ff"), ("L5", "#ff8fa3")):
        rp = rotating_frame(holds[st], 2)
        ax.plot(rp[:, 0], rp[:, 1], color=col, lw=1.3, alpha=0.9)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.2, 1.2)
    ax.tick_params(colors="#556"); [s.set_color("#334") for s in ax.spines.values()]

    ax2 = fig.add_subplot(gs[1]); ax2.set_facecolor("#0a0e17")
    ax2.set_title("resonant angle over 80 orbits — it just holds",
                  color="w", fontsize=12)
    for st, col in (("L4", "#54d1ff"), ("L5", "#ff8fa3")):
        r = holds[st]
        ax2.plot(r["t"] / memory.PERIOD, r["phi"], color=col, lw=1.1,
                 label=f"{st}  (reads {'1' if st=='L4' else '0'})")
    ax2.axhline(60, color="#54d1ff", ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(-60, color="#ff8fa3", ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(180, color="#556", ls="--", lw=0.8); ax2.axhline(-180, color="#556", ls="--", lw=0.8)
    ax2.set_xlabel("orbits", color="w"); ax2.set_ylabel("λ(particle) − λ(planet)  [deg]", color="w")
    ax2.set_ylim(-200, 200)
    ax2.tick_params(colors="#889"); [s.set_color("#334") for s in ax2.spines.values()]
    ax2.legend(facecolor="#141a2a", labelcolor="w", edgecolor="#334", fontsize=10)

    fig.suptitle("Orbital Memory · a gravitational flip-flop", color="w", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "flipflop.png", dpi=140, facecolor=fig.get_facecolor())
    print(f"wrote {OUT / 'flipflop.png'}")


if __name__ == "__main__":
    main()
