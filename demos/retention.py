"""RETENTION, as imagery: how long a bit lasts vs how deep it is.

Escape time (orbits) as a function of seeded libration amplitude. Below the
chaotic separatrix layer the bit is a censored survivor — Nekhoroshev-stable,
lasting far beyond any feasible run (shown as the shaded 'stable' floor). As the
amplitude climbs toward the L3 separatrix (analytic ~78 deg), the escape time
collapses. Honest: only the escape band is integrated; the deep regime's true
lifetime is theory, and real Trojans are metastable at 10 kyr-100 Myr.

Writes docs/retention.png.  Usage: python -m demos.retention
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, DOCS, ERASE, GOOD, GRID, GROUND, L4C, TEXT_W
from orbital import hamiltonian, retention

HORIZON = 300
AMPS = [50, 60, 66, 70, 74, 78, 82, 86]


def main():
    DOCS.mkdir(exist_ok=True)
    scan = retention.escape_scan(AMPS, max_orbits=HORIZON)
    sep = hamiltonian.separatrix_amplitude()

    fig, ax = plt.subplots(figsize=(8.4, 5.0), facecolor=GROUND)
    ax.set_facecolor(GROUND)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.tick_params(colors=DIM)

    # censored survivors: draw at the horizon with an up-arrow
    esc_x, esc_y, cen_x = [], [], []
    for r in scan:
        if r["censored"]:
            cen_x.append(r["amp"])
        else:
            esc_x.append(r["amp"]); esc_y.append(r["escape_orbits"])

    ax.axhspan(HORIZON, HORIZON * 1.5, color=GOOD, alpha=0.08)
    ax.text(52, HORIZON * 1.22, "censored survivors — Nekhoroshev-stable, "
            "lifetime ≫ integration (theory)", color=GOOD, fontsize=9, va="center")
    for x in cen_x:
        ax.annotate("", xy=(x, HORIZON * 1.4), xytext=(x, HORIZON),
                    arrowprops=dict(arrowstyle="-|>", color=GOOD, lw=1.6))
        ax.scatter([x], [HORIZON], color=GOOD, s=45, zorder=5)
    ax.scatter(esc_x, esc_y, color=ERASE, s=55, zorder=5, label="escaped (integrated)")
    ax.plot(esc_x, esc_y, color=ERASE, lw=1.0, alpha=0.6)

    ax.axvline(sep, color=TEXT_W, lw=1.3, ls=(0, (5, 3)))
    ax.text(sep + 0.6, HORIZON * 0.5, f"L3 separatrix\n(analytic {sep:.0f}°)",
            color=TEXT_W, fontsize=9, va="center")

    ax.set_xlabel("seeded libration amplitude  (deg)", color=TEXT_W)
    ax.set_ylabel("escape time  (orbits)", color=TEXT_W)
    ax.set_ylim(0, HORIZON * 1.5)
    ax.set_title("how long does a bit last? — deep bits survive, "
                 "near-separatrix bits escape fast",
                 color=TEXT_W, fontsize=12.5, pad=12)
    fig.tight_layout()
    out = DOCS / "retention.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    n_esc = len(esc_x)
    print(f"wrote {out}  ({n_esc} escaped, {len(cen_x)} censored survivors)")


if __name__ == "__main__":
    main()
