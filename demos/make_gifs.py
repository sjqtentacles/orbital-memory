"""Render the 2D memory story as a GIF: the bit holds, then a kick erases it.

Rotating-frame view. A particle librates in the L4 tadpole (bit = 1) for a few
orbits, then a super-threshold kick pushes it across the separatrix and it
circulates away — the bit is lost. Shows the noise margin as motion.

Writes docs/flipflop_2d.gif.  Usage: python -m demos.make_gifs
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import DIM, DOCS, ERASE, GROUND, L4C, L5C, PLANET, STAR, optimize_gif
from orbital import memory, nbody, rotating

HOLD_ORB, POST_ORB, FPO = 7, 14, 9  # orbits before/after kick, frames per orbit


def main():
    DOCS.mkdir(exist_ok=True)
    # hold, then kick across the separatrix, then watch it circulate
    cell = memory.make_cell("L4", libration_deg=6.0)
    pre = nbody.integrate(cell, HOLD_ORB * memory.PERIOD, n_samples=HOLD_ORB * FPO)
    kicked = memory.kicked_cell(0.05, libration_deg=6.0,
                                settle_periods=HOLD_ORB)
    # t_end is ABSOLUTE time: the post-kick leg spans POST_ORB full orbits
    post = nbody.integrate(kicked, (HOLD_ORB + POST_ORB) * memory.PERIOD,
                           n_samples=POST_ORB * FPO,
                           t0=HOLD_ORB * memory.PERIOD)

    p_pre = rotating.to_rotating_frame(pre, 2)
    p_post = rotating.to_rotating_frame(post, 2)
    path = np.vstack([p_pre, p_post])
    n_hold = len(p_pre)
    total = len(path)

    star, planet = memory.primaries()
    l4, _ = memory.lagrange_state(memory.MU, "L4")
    l5, _ = memory.lagrange_state(memory.MU, "L5")

    fig, ax = plt.subplots(figsize=(6.4, 6.0), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.3, 1.4); ax.axis("off")

    def base():
        ax.plot(0, 0, "+", color="#39445f", ms=9)
        for gx in np.linspace(-1.2, 1.2, 7):
            ax.axvline(gx, color="#121a2c", lw=0.7, zorder=0)
            ax.axhline(gx, color="#121a2c", lw=0.7, zorder=0)
        ax.plot(*star["pos"], "o", color=STAR, ms=16, zorder=5)
        ax.plot(*planet["pos"], "o", color=PLANET, ms=8, zorder=5)
        for pt, lab, col in ((l4, "L4 = 1", L4C), (l5, "L5 = 0", L5C)):
            ax.plot(*pt, "x", color=col, ms=9, mew=2)
            ax.annotate(lab, pt, textcoords="offset points", xytext=(9, 6),
                        color=col, fontsize=11)

    trail = LineCollection([], zorder=3)
    ax.add_collection(trail)
    head = ax.scatter([], [], s=80, zorder=6, edgecolors="none")
    glow = ax.scatter([], [], s=320, alpha=0.25, zorder=5, edgecolors="none")
    state = ax.text(0.03, 0.96, "", transform=ax.transAxes, fontsize=15,
                    va="top", weight="bold", family="monospace")
    ax.text(0.5, 0.02, "rotating frame · a Trojan memory bit and its noise margin",
            transform=ax.transAxes, color=DIM, fontsize=9, ha="center")

    def draw(f):
        # redraw base each frame cheaply by clearing artists we own
        for c in list(ax.collections):
            if c not in (trail, head, glow):
                c.remove()
        for ln in list(ax.lines):
            ln.remove()
        for tx in list(ax.texts):
            if tx not in (state,):
                tx.remove()
        base()
        erased = f >= n_hold
        col = ERASE if erased else L4C
        lo = max(0, f - 90)
        seg = path[lo:f + 1]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            a = np.linspace(0.05, 1, len(segs))
            rgb = np.array(matplotlib.colors.to_rgb(col))
            trail.set_segments(segs)
            trail.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), a]))
            trail.set_linewidth(2.6)
        head.set_offsets([path[f]]); head.set_color(col)
        glow.set_offsets([path[f]]); glow.set_color(col)
        if erased:
            state.set_text("ERASED  (bit lost)"); state.set_color(ERASE)
        else:
            state.set_text("HOLDING  bit = 1"); state.set_color(L4C)
        return []

    anim = FuncAnimation(fig, draw, frames=total, interval=55, blit=False)
    out = DOCS / "flipflop_2d.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
