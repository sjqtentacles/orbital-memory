"""STABILITY MAP: where a co-orbital bit survives, in (mass ratio, amplitude).

For a grid of mass ratios and seeded libration amplitudes, integrate a tadpole
and ask whether it is still a stored bit afterwards. The triangular points are
linearly stable only for mu < 0.0385 (the Gascheau/Routh limit) — beyond it
L4/L5 come apart and no bit survives. Real co-orbitals (Jupiter's Trojans,
Earth's, Saturn's moons) all sit far to the stable side; they are marked.

Writes docs/stability.png.  Usage: python -m demos.stability_map
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, DOCS, ERASE, GOOD, GRID, GROUND, TEXT_W
from orbital import memory, nbody

ROUTH = 0.0385
MUS = np.geomspace(8e-7, 0.12, 20)
AMPS = np.linspace(6, 66, 7)
REAL = [("Jupiter Trojans", memory.MU_SUN_JUPITER),
        ("Earth Trojan", memory.MU_SUN_EARTH),
        ("Saturn co-orbitals", memory.MU_SATURN_TETHYS)]


def _survives(mu, amp_deg):
    cell = memory.make_cell("L4", mu=mu, libration_deg=amp_deg)
    res = nbody.integrate(cell, 25 * memory.PERIOD, n_samples=1200)
    return memory.classify(memory.resonant_angle(res))[0] == "L4"


def main():
    DOCS.mkdir(exist_ok=True)
    grid = np.array([[_survives(mu, amp) for mu in MUS] for amp in AMPS])

    fig, ax = plt.subplots(figsize=(8.6, 5.0), facecolor=GROUND)
    ax.set_facecolor(GROUND)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.tick_params(colors=DIM)

    from matplotlib.colors import ListedColormap
    ax.pcolormesh(MUS, AMPS, grid, cmap=ListedColormap([ERASE, GOOD]),
                  shading="nearest", alpha=0.75)

    ax.axvline(ROUTH, color=TEXT_W, lw=1.4, ls=(0, (5, 3)))
    ax.text(ROUTH * 1.04, AMPS[-1] - 4, "Routh limit  μ = 0.0385\n"
            "→ L4/L5 unstable", color=TEXT_W, fontsize=9, va="top")

    for name, mu in REAL:
        ax.axvline(mu, color=DIM, lw=1.0)
        ax.text(mu, AMPS[0] - 3.5, name, color=DIM, fontsize=8.2, rotation=90,
                va="top", ha="center")

    ax.set_xscale("log")
    ax.set_xlabel("mass ratio  μ  =  m₂ / (m₁ + m₂)", color=TEXT_W)
    ax.set_ylabel("seeded libration amplitude  (deg)", color=TEXT_W)
    ax.set_title("where a co-orbital bit survives — green stable, red lost",
                 color=TEXT_W, fontsize=12.5, pad=12)
    fig.tight_layout()
    out = DOCS / "stability.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    n_ok = int(grid.sum())
    print(f"wrote {out}  ({n_ok}/{grid.size} grid cells stable)")


if __name__ == "__main__":
    main()
