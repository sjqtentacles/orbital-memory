"""ERASE and the rewrite wall, as imagery.

One continuous run: a bit is written (pinch), stored, then the erase pulse
shrinks the secondary back down and the separatrix releases the moonlet into
a wide blank horseshoe — visibly WIDER than the medium it was written from.
That inflation is the rewrite wall: the released horseshoe (da ~ 0.07) is
past the recapture ceiling (~0.045), so this memory is write-once/erase-once.

Writes docs/erase.gif.  Usage: python -m demos.rewrite_demo
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import (DIM, DOCS, ERASE, GROUND, HORSE, L4C, L5C, PLANET,
                         STAR, optimize_gif)
from orbital import memory, rewrite, theory, write

FRAMES = 170


def main():
    DOCS.mkdir(exist_ok=True)
    run = rewrite.erased_prefix("1", n_samples=FRAMES)
    m = run["markers"]
    da = rewrite.horseshoe_offset(run)
    print(f"released horseshoe offset da = {da:.3f} "
          f"(captured from {write.DA}; ceiling {rewrite.RECAPTURE_CEILING})")

    xy, t, mus = run["xy"], run["t"], run["mu"]
    l4 = theory.lagrange_points(memory.MU)["L4"]
    l5 = theory.lagrange_points(memory.MU)["L5"]

    fig, ax = plt.subplots(figsize=(6.4, 6.2), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.55); ax.axis("off")

    trail = LineCollection([], zorder=4); ax.add_collection(trail)
    head = ax.scatter([], [], s=70, zorder=6, edgecolors="none")
    planet_dot = ax.scatter([], [], s=40, color=PLANET, zorder=5,
                            edgecolors="none")
    state_txt = ax.text(0.03, 0.97, "", transform=ax.transAxes, fontsize=13,
                        va="top", weight="bold", family="monospace")
    mu_txt = ax.text(0.03, 0.905, "", transform=ax.transAxes, fontsize=9.5,
                     va="top", color=DIM, family="monospace")

    def draw(f):
        for ln in list(ax.lines):
            ln.remove()
        for txt in list(ax.texts):
            if txt not in (state_txt, mu_txt):
                txt.remove()
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(np.cos(th), np.sin(th), color="#1a2340", lw=0.9, zorder=1)
        ax.plot(0, 0, marker="o", ms=15, color=STAR, zorder=5)
        for pt, lab, col in ((l4, "L4", L4C), (l5, "L5", L5C)):
            ax.plot(*pt, "x", color=col, ms=8, mew=1.8, zorder=5)
            ax.annotate(lab, pt, xytext=(8, 6), textcoords="offset points",
                        color=col, fontsize=10)
        now = t[f]
        if now < m["write1"] + write.T_RAMP:
            phase, col = "WRITE · pinch in progress", "#ffe08a"
        elif now < m["erase"]:
            phase, col = "STORED · bit = 1", L4C
        elif now < m["erase"] + rewrite.T_ERASE:
            phase, col = "ERASE · secondary shrinking", ERASE
        else:
            phase, col = "BLANK · released horseshoe (inflated)", HORSE
        state_txt.set_text(phase); state_txt.set_color(col)
        mu_txt.set_text(f"mass ratio μ = {mus[f]:.5f}")
        lo = max(0, f - 55)
        seg = xy[lo:f + 1]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            aa = np.linspace(0.05, 1, len(segs))
            rgb = np.array(matplotlib.colors.to_rgb(col))
            trail.set_segments(segs)
            trail.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), aa]))
            trail.set_linewidth(2.4)
        head.set_offsets([xy[f]]); head.set_color(col)
        planet_dot.set_offsets([[1 - mus[f], 0]])
        planet_dot.set_sizes([40 + 60 * (mus[f] / memory.MU)])
        ax.text(0.5, 0.015, "erasing frees the bit — but inflates the medium "
                "past recapture: write-once memory",
                transform=ax.transAxes, color=DIM, fontsize=8.5, ha="center")
        return []

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=60, blit=False)
    out = DOCS / "erase.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
