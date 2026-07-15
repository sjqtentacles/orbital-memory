"""Demo 01: a gravitational memory bit that holds — and its noise margin.

A test particle librates in a tadpole around the leading (L4) or trailing (L5)
Trojan point of a star+planet system. Which island it's in is a nonvolatile
bit: the two islands are separated by a separatrix, so any perturbation smaller
than a threshold leaves the bit intact, and a larger one erases it.

The protection is TOPOLOGICAL, not dissipative — L4/L5 are Coriolis-stabilized
equilibria, so adding drag would destabilize them. The bit is held by the
phase-space geometry (an invariant island), the KAM/Nekhoroshev way.

Writes out/flipflop.png.  Usage: python -m demos.flipflop_demo
"""

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, GROUND, L4C, L5C, PANEL, PLANET, STAR
from orbital import memory, nbody, rotating, theory

OUT = pathlib.Path(__file__).resolve().parent.parent / "out"


def noise_margin(periods=50, samples=2500):
    """Kick an L4 bit with growing tangential impulses; find the fraction of
    orbital speed at which the tadpole is destroyed (bit erased)."""
    rows = []
    for frac in np.arange(0.005, 0.055, 0.005):
        res = nbody.integrate(memory.kicked_cell(frac),
                              periods * memory.PERIOD, n_samples=samples)
        label, center, amp = memory.classify(memory.resonant_angle(res))
        rows.append((frac, label, amp))
    return rows


def main():
    OUT.mkdir(exist_ok=True)

    print("=== HOLD: two nonvolatile states ===")
    holds = {}
    for st in ("L4", "L5"):
        res = memory.hold(st, periods=80, libration_deg=6.0)
        phi = res["phi"]
        label, center, amp = memory.classify(phi)
        holds[st] = res
        cj = theory.jacobi_of(res)
        print(f"  {st}: librates around {center:+5.1f} deg (amp {amp:.1f}), "
              f"reads {label!r} for 80 orbits")
        print(f"      moonlet Jacobi constant C_J = {cj.mean():.6f} "
              f"(conserved to {cj.max() - cj.min():.0e}; analytic C_L4 = "
              f"{theory.C_L4():.6f}); primaries' energy drift "
              f"{res['energy_drift']:.0e}")
    assert memory.classify(holds["L4"]["phi"])[0] == "L4"
    assert memory.classify(holds["L5"]["phi"])[0] == "L5"
    print("  -> the bit holds. (Jupiter's real Trojans have held theirs for 4.5 Gyr.)")

    print("\n=== NOISE MARGIN: how hard a kick erases the bit ===")
    rows = noise_margin()
    threshold = next((f for f, lab, _ in rows if lab == "erased"), None)
    for frac, label, amp in rows:
        bar = "kept " if label != "erased" else "ERASED"
        print(f"  kick {frac*100:4.1f}% of orbital speed -> {bar} "
              f"(libration amp {amp:5.1f} deg, reads {label})")
    print(f"  -> separatrix threshold ~{threshold*100:.1f}% of orbital speed "
          f"(promoted to memory.ERASE_KICK = {memory.ERASE_KICK}): a real "
          f"noise margin, set by phase-space geometry, not dissipation.")

    star, planet = memory.primaries()
    l4, _ = memory.lagrange_state(memory.MU, "L4")
    l5, _ = memory.lagrange_state(memory.MU, "L5")

    # --- figure: rotating-frame tadpoles + phi(t) traces ---
    fig = plt.figure(figsize=(12, 5.2), facecolor=GROUND)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.22)

    ax = fig.add_subplot(gs[0]); ax.set_facecolor(GROUND)
    ax.set_aspect("equal")
    ax.set_title("rotating frame — the bit is *where* it orbits",
                 color="w", fontsize=12)
    ax.plot(0, 0, "+", color="#39445f", ms=10)  # barycenter
    ax.plot(*star["pos"], "o", color=STAR, ms=16, label="star")
    ax.plot(*planet["pos"], "o", color=PLANET, ms=9, label="planet")
    for pt, lab, col in ((l4, "L4 = 1", L4C), (l5, "L5 = 0", L5C)):
        ax.plot(*pt, "x", color=col, ms=9)
        ax.annotate(lab, pt, textcoords="offset points", xytext=(8, 6),
                    color=col, fontsize=11)
    for st, col in (("L4", L4C), ("L5", L5C)):
        rp = rotating.to_rotating_frame(holds[st], 2)
        ax.plot(rp[:, 0], rp[:, 1], color=col, lw=1.3, alpha=0.9)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.2, 1.2)
    ax.tick_params(colors="#556")
    for s in ax.spines.values():
        s.set_color("#334")

    ax2 = fig.add_subplot(gs[1]); ax2.set_facecolor(GROUND)
    ax2.set_title("resonant angle over 80 orbits — it just holds",
                  color="w", fontsize=12)
    for st, col in (("L4", L4C), ("L5", L5C)):
        r = holds[st]
        ax2.plot(r["t"] / memory.PERIOD, r["phi"], color=col, lw=1.1,
                 label=f"{st}  (reads {'1' if st=='L4' else '0'})")
    ax2.axhline(60, color=L4C, ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(-60, color=L5C, ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(180, color="#556", ls="--", lw=0.8)
    ax2.axhline(-180, color="#556", ls="--", lw=0.8)
    ax2.set_xlabel("orbits", color="w")
    ax2.set_ylabel("λ(particle) − λ(planet)  [deg]", color="w")
    ax2.set_ylim(-200, 200)
    ax2.tick_params(colors="#889")
    for s in ax2.spines.values():
        s.set_color("#334")
    ax2.legend(facecolor=PANEL, labelcolor="w", edgecolor="#334", fontsize=10)

    fig.suptitle("Orbital Memory · a gravitational flip-flop", color="w",
                 fontsize=14)
    from demos.style import DOCS
    DOCS.mkdir(exist_ok=True)
    for target in (OUT / "flipflop.png", DOCS / "flipflop.png"):
        fig.savefig(target, dpi=140, facecolor=fig.get_facecolor(),
                    bbox_inches="tight")
        print(f"wrote {target}")


if __name__ == "__main__":
    main()
