"""PHASE PORTRAIT: the co-orbital pendulum's whole phase space, one figure.

The stored bit is a slow pendulum whose coordinate is the resonant angle phi
and whose momentum is the radial offset da = r - 1. Nested closed loops around
L4 (+60, bit 1) and L5 (-60, bit 0) are TADPOLES; the big loops that wrap all
the way around through L3 (180) are HORSESHOES; beyond them the orbit
CIRCULATES. The tadpole islands are the two memory states; the separatrix
between them is why a small perturbation can't flip the bit.

Drawn at a small mass ratio where all three families coexist cleanly (the same
regime as the anatomy figure).

Writes docs/phase_portrait.png.  Usage: python -m demos.phase_portrait
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import (CIRC, DIM, DOCS, GRID, GROUND, HORSE, L4C, L5C,
                         TEXT_W)
from orbital import memory, rotating

MU = 3e-4                       # small ratio: tadpoles + horseshoes coexist
LIB = 2 * np.pi / np.sqrt(6.75 * MU)     # linear libration period (nondim)


def _phi_da(run):
    ang = np.degrees(np.arctan2(run["xy"][:, 1], run["xy"][:, 0]))
    ang = ((ang + 180) % 360) - 180
    r = np.hypot(run["xy"][:, 0], run["xy"][:, 1])
    return ang, r - 1.0


def _tadpole(point, libration_deg):
    st = rotating.lagrange_tadpole(point, libration_deg=libration_deg, mu=MU)
    return _phi_da(rotating.integrate(st, 1.15 * LIB, n_samples=2500, mu=MU))


def main():
    DOCS.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 5.2), facecolor=GROUND)
    ax.set_facecolor(GROUND)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.tick_params(colors=DIM)

    # horseshoe / circulation zones, labelled (their orbits are chaotic to
    # trace at this mu; the memory lives in the tadpole islands below)
    ax.axhspan(0.028, 0.075, color=HORSE, alpha=0.07)
    ax.axhspan(-0.075, -0.028, color=HORSE, alpha=0.07)
    ax.axhspan(0.075, 0.11, color=CIRC, alpha=0.06)
    ax.axhspan(-0.11, -0.075, color=CIRC, alpha=0.06)
    ax.text(0, 0.093, "circulation  ·  blank", color=CIRC, fontsize=9,
            ha="center", va="center")
    ax.text(0, 0.051, "horseshoe  ·  erased  (sweeps through L3)", color=HORSE,
            fontsize=9, ha="center", va="center")

    # tadpoles: nested loops around each island — the two memory states
    for point, col in (("L4", L4C), ("L5", L5C)):
        for ld in np.linspace(5, 38, 6):
            ang, da = _tadpole(point, ld)
            ax.plot(ang, da, color=col, lw=0.9, alpha=0.9)

    ax.axvline(0, color=DIM, lw=0.6, alpha=0.4)
    for x, lab, c in ((60, "L4 · bit 1", L4C), (-60, "L5 · bit 0", L5C),
                      (180, "L3", DIM), (-180, "L3", DIM)):
        ax.scatter([x], [0], color=c, s=26, zorder=6)
        ax.annotate(lab, (x, 0), textcoords="offset points", xytext=(0, 8),
                    color=c, fontsize=8.5, ha="center")

    ax.set_xlim(-185, 185)
    ax.set_ylim(-0.11, 0.11)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_xlabel("resonant angle  φ  (deg)", color=TEXT_W)
    ax.set_ylabel("radial offset  da = r − 1", color=TEXT_W)
    ax.set_title("phase portrait of the co-orbital pendulum — "
                 "tadpoles (bits), horseshoes, circulation",
                 color=TEXT_W, fontsize=12.5, pad=12)
    fig.tight_layout()
    out = DOCS / "phase_portrait.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
