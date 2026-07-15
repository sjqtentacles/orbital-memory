"""VALIDATION, as imagery: the simulated bit against the real sky.

Runs the cell at the real Sun-Jupiter mass ratio, converts the clock to years
with orbital.units, and shows the resonant angle librating with a ~148-year
period — the period actually observed for Jupiter's Trojan asteroids. Nothing
is fit to that number; it falls out of Newtonian gravity at the real ratio.

Writes docs/validation.png.  Usage: python -m demos.validation
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, DOCS, GOOD, GRID, GROUND, L4C, PLANET, STAR, TEXT_W
from orbital import memory, nbody, rotating, units

SYS = units.SUN_JUPITER
OBSERVED = (147.0, 160.0)   # observed Jupiter-Trojan long-period libration band


def main():
    DOCS.mkdir(exist_ok=True)
    mu = memory.MU_SUN_JUPITER
    cell = memory.make_cell("L4", mu=mu, libration_deg=12.0)
    res = nbody.integrate(cell, 45 * memory.PERIOD, n_samples=6000)
    phi = memory.resonant_angle(res)
    years = SYS.years(res["t"])

    t_lib = memory.measured_libration_period(res["t"], phi)
    lib_years = SYS.years(t_lib)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.8), facecolor=GROUND,
                                   gridspec_kw={"width_ratios": [1.55, 1.0]})
    for ax in (axL, axR):
        ax.set_facecolor(GROUND)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.tick_params(colors=DIM)

    # left: phi(t) in years, the librating bit
    axL.axhspan(60 - 12, 60 + 12, color=L4C, alpha=0.06)
    axL.axhline(60, color=DIM, lw=0.8, ls=(0, (4, 3)))
    axL.plot(years, phi, color=L4C, lw=1.6)
    axL.set_xlabel("time  (years)", color=TEXT_W)
    axL.set_ylabel("resonant angle  φ  (deg)", color=TEXT_W)
    axL.set_title("a Jupiter Trojan bit, in real time", color=TEXT_W, fontsize=12.5)
    axL.text(0.03, 0.06,
             f"simulated libration period:  {lib_years:.1f} yr\n"
             f"observed (Jupiter Trojans):  {OBSERVED[0]:.0f}–{OBSERVED[1]:.0f} yr\n"
             f"Jupiter's year:  {SYS.period_years:.2f} yr    ·    "
             f"orbital speed:  {SYS.kmps(1.0):.1f} km/s",
             transform=axL.transAxes, color=GOOD, fontsize=9.5, va="bottom",
             family="monospace")

    # right: the tadpole in the rotating frame, axes in AU
    ref = rotating.to_rotating_frame(res, 2)
    axR.set_aspect("equal")
    axR.plot(SYS.au(ref[:, 0]), SYS.au(ref[:, 1]), color=L4C, lw=1.0, alpha=0.9)
    axR.scatter([0], [0], c=STAR, s=90, zorder=5, edgecolors="none")
    axR.scatter([SYS.au(1 - mu)], [0], c=PLANET, s=45, zorder=5, edgecolors="none")
    l4 = np.array([0.5 - mu, np.sqrt(3) / 2])
    axR.scatter([SYS.au(l4[0])], [SYS.au(l4[1])], marker="+", c=DIM, s=80)
    axR.set_xlabel("x  (AU)", color=TEXT_W)
    axR.set_ylabel("y  (AU)", color=TEXT_W)
    axR.set_title("tadpole around L4", color=TEXT_W, fontsize=12.5)

    fig.suptitle("validated against the sky — Newtonian gravity at the real "
                 "Sun–Jupiter ratio reproduces the observed Trojan libration",
                 color=TEXT_W, fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = DOCS / "validation.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {out}  (libration {lib_years:.1f} yr vs observed "
          f"{OBSERVED[0]:.0f}-{OBSERVED[1]:.0f} yr)")


if __name__ == "__main__":
    main()
