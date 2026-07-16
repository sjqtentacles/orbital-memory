"""THE ANALYTIC BACKBONE, as imagery: theory vs the sim.

Left: the averaged co-orbital potential well V(phi) — a bit librates in it, with
L4 at the floor and L3 the separatrix barrier (3*mu). Right: the closed-form
libration period vs amplitude (quadrature on that well) overlaid on the periods
measured from full n-body integrations. They agree to ~1.5% across the tadpole
range — the sim verifies the theory.

Writes docs/hamiltonian.png.  Usage: python -m demos.hamiltonian
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, DOCS, GOOD, GRID, GROUND, L4C, L5C, TEXT_W
from orbital import hamiltonian as H, memory, nbody

AMPS = [2, 6, 15, 30, 45, 60]


def main():
    DOCS.mkdir(exist_ok=True)
    mu = memory.MU

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.8), facecolor=GROUND)
    for ax in (axL, axR):
        ax.set_facecolor(GROUND)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.tick_params(colors=DIM)

    # left: the potential well V(phi) in units of mu
    phi = np.linspace(3, 357, 1000)
    V = np.array([H.potential(p, mu) for p in phi]) / mu
    axL.plot(phi, V, color=L4C, lw=1.8)
    for x, lab, c in ((60, "L4", L4C), (300, "L5", L5C), (180, "L3\n(separatrix)", DIM)):
        axL.scatter([x], [H.potential(x, mu) / mu], color=c, s=45, zorder=5)
        axL.annotate(lab, (x, H.potential(x, mu) / mu), textcoords="offset points",
                     xytext=(0, 8), color=c, fontsize=9, ha="center")
    axL.axhline(3.0, color=DIM, lw=0.7, ls=(0, (4, 3)))
    axL.text(20, 3.15, "barrier = 3μ", color=DIM, fontsize=9)
    axL.set_ylim(-0.5, 6)
    axL.set_xlabel("resonant angle  φ  (deg)", color=TEXT_W)
    axL.set_ylabel("potential  V(φ) / μ", color=TEXT_W)
    axL.set_title("the averaged co-orbital well", color=TEXT_W, fontsize=12.5)

    # right: period(amp) analytic curve vs measured points
    amp_fine = np.linspace(0.5, 66, 120)
    Ta = [H.libration_period(a, mu) for a in amp_fine]
    axR.plot(amp_fine, Ta, color=GOOD, lw=1.8, label="analytic (Hamiltonian)")
    meas = []
    for a in AMPS:
        res = nbody.integrate(memory.make_cell("L4", mu=mu, libration_deg=a),
                              60 * memory.PERIOD, n_samples=5000)
        meas.append(memory.measured_libration_period(res["t"], memory.resonant_angle(res)))
    axR.scatter(AMPS, meas, color=L4C, s=55, zorder=5, label="measured (n-body)")
    axR.set_xlabel("libration amplitude  (deg)", color=TEXT_W)
    axR.set_ylabel("libration period  (nondim, planet 2π)", color=TEXT_W)
    axR.set_title("period vs amplitude — theory verified by sim", color=TEXT_W, fontsize=12.5)
    leg = axR.legend(facecolor=GROUND, edgecolor=GRID, labelcolor=TEXT_W, fontsize=9)

    fig.suptitle("the analytic backbone — closed forms matched to the simulation",
                 color=TEXT_W, fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = DOCS / "hamiltonian.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    err = max(abs(H.libration_period(a, mu) / m - 1) for a, m in zip(AMPS, meas))
    print(f"wrote {out}  (max period error vs sim: {err*100:.1f}%)")


if __name__ == "__main__":
    main()
