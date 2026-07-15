"""SATURN, as imagery: the onset of Trojan erosion.

Two lanes, same Jupiter Trojan bit, same clock. Left: real Saturn added as a
fourth body — its periodic tug pumps the Trojan's libration wider and wider.
Right: no Saturn, the bit sits nearly still. Full erosion (shaking Trojans
loose) is a Gyr process; what a short run shows is its first sign, and that is
what the panels claim.

Writes docs/erosion.gif.  Usage: python -m demos.saturn_erosion
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import (DIM, DOCS, ERASE, GOOD, GROUND, L4C, PLANET, STAR,
                         optimize_gif)
from orbital import memory, nbody, rotating, units

FRAMES = 150
PERIODS = 130
SYS = units.SUN_JUPITER


def main():
    DOCS.mkdir(exist_ok=True)
    t_end = PERIODS * memory.PERIOD
    with_sat = nbody.integrate(
        memory.sun_jupiter_saturn_cell("L4", libration_deg=10.0, with_saturn=True),
        t_end, n_samples=FRAMES)
    ctrl = nbody.integrate(
        memory.sun_jupiter_saturn_cell("L4", libration_deg=10.0, with_saturn=False),
        t_end, n_samples=FRAMES)
    moon_s = rotating.to_rotating_frame(with_sat, 2)
    moon_c = rotating.to_rotating_frame(ctrl, 2)
    phi_s = memory.resonant_angle(with_sat)
    phi_c = memory.resonant_angle(ctrl)
    l4 = np.array([0.5 - memory.MU_SUN_JUPITER, np.sqrt(3) / 2])
    years = SYS.years(with_sat["t"])

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6), facecolor=GROUND)
    lanes = []
    for ax, title in zip(axes, ("WITH Saturn", "no Saturn (control)")):
        ax.set_facecolor(GROUND); ax.set_aspect("equal")
        ax.set_xlim(-0.2, 1.15); ax.set_ylim(0.1, 1.25); ax.axis("off")
        ax.set_title(title, color="w", fontsize=12, family="monospace")
        trail = LineCollection([], zorder=4); ax.add_collection(trail)
        head = ax.scatter([], [], s=55, zorder=6, edgecolors="none")
        txt = ax.text(0.04, 0.97, "", transform=ax.transAxes, fontsize=11,
                      va="top", weight="bold", family="monospace")
        lanes.append((ax, trail, head, txt))

    def base(ax):
        th = np.linspace(0.3, 1.9, 120)
        ax.plot(np.cos(th), np.sin(th), color="#1a2340", lw=0.8, zorder=1)
        ax.plot(1 - memory.MU_SUN_JUPITER, 0, marker="o", ms=6, color=PLANET,
                zorder=5)
        ax.plot(*l4, "x", color=L4C, ms=8, mew=1.6, zorder=5)

    def fade(lc, pts, f, color, span=60):
        seg = pts[max(0, f - span):f + 1]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            aa = np.linspace(0.05, 1, len(segs))
            rgb = np.array(matplotlib.colors.to_rgb(color))
            lc.set_segments(segs)
            lc.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), aa]))
            lc.set_linewidth(2.0)

    def amp_to(phi, f):
        s = phi[:f + 1] - 60.0
        return float(np.max(np.abs(s))) if f > 0 else 0.0

    def draw(f):
        for (ax, trail, head, txt), moon, phi, col in (
                (lanes[0], moon_s, phi_s, ERASE),
                (lanes[1], moon_c, phi_c, GOOD)):
            for ln in list(ax.lines):
                ln.remove()
            base(ax)
            fade(trail, moon, f, col)
            head.set_offsets([moon[f]]); head.set_color(col)
            txt.set_text(f"amp {amp_to(phi, f):4.1f}°\n{years[f]:4.0f} yr")
            txt.set_color(col)
        fig.suptitle("Saturn erodes the Trojan swarm — the onset, over "
                     f"{years[-1]:.0f} years (full erosion takes Gyr)",
                     color="w", fontsize=12, y=0.05)
        return []

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=55, blit=False)
    out = DOCS / "erosion.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
