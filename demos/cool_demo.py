"""COOL, as imagery: a wide written bit tightened by station-keeping burns.

Rotating-frame animation: the freshly written ±66° tadpole sweeps its wide
arc; each tangential burn (flash) shrinks the radial offset; the trail
tightens into a deep ~25° bit parked by L4. An inset counts the burns and
tracks the libration amplitude.

Writes docs/cool.gif.  Usage: python -m demos.cool_demo
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import (DIM, DOCS, GOOD, GROUND, HORSE, L4C, PLANET, STAR,
                         optimize_gif)
from orbital import cool, memory, theory, write

FRAMES = 170


def main():
    DOCS.mkdir(exist_ok=True)
    wide = write.write_bit("1")
    cooled = cool.cool(wide["yf"], t0=wide["t"][-1])
    label, center, amp = memory.classify(cooled["phi"][-800:])
    assert label == "L4", f"cool demo produced {label}"
    print(f"cooled: {cooled['n_kicks']} burns -> amp {amp:.1f} deg")

    # animate the cooling run (subsample to FRAMES)
    idx = np.linspace(0, len(cooled["t"]) - 1, FRAMES).astype(int)
    t_f = cooled["t"][idx]
    xy_f = cooled["xy"][idx]
    kick_ts = cooled["kick_times"]

    l4 = theory.lagrange_points(memory.MU)["L4"]

    fig, ax = plt.subplots(figsize=(6.4, 6.2), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.45); ax.axis("off")

    trail = LineCollection([], zorder=4); ax.add_collection(trail)
    head = ax.scatter([], [], s=70, zorder=6, edgecolors="none")
    burn = ax.scatter([], [], s=500, zorder=5, edgecolors="none",
                      color=HORSE, alpha=0.0)
    state_txt = ax.text(0.03, 0.97, "", transform=ax.transAxes, fontsize=13.5,
                        va="top", weight="bold", family="monospace")
    sub_txt = ax.text(0.03, 0.90, "", transform=ax.transAxes, fontsize=9.5,
                      va="top", color=DIM, family="monospace")

    def draw(f):
        for ln in list(ax.lines):
            ln.remove()
        for txt in list(ax.texts):
            if txt not in (state_txt, sub_txt):
                txt.remove()
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(np.cos(th), np.sin(th), color="#1a2340", lw=0.9, zorder=1)
        ax.plot(0, 0, marker="o", ms=15, color=STAR, zorder=5)
        ax.plot(1 - memory.MU, 0, marker="o", ms=7, color=PLANET, zorder=5)
        ax.plot(*l4, "x", color=L4C, ms=8, mew=1.8, zorder=5)
        ax.annotate("L4", l4, xytext=(8, 6), textcoords="offset points",
                    color=L4C, fontsize=10)
        now = t_f[f]
        burns_done = int(np.sum(kick_ts <= now))
        cooling_on = len(kick_ts) and now <= kick_ts[-1] + memory.PERIOD
        recent_burn = len(kick_ts) and np.any(
            np.abs(kick_ts - now) < 1.5 * memory.PERIOD)
        col = HORSE if cooling_on else GOOD
        if burns_done == 0:
            state_txt.set_text("WIDE BIT · fresh from the write pinch")
            state_txt.set_color(L4C)
        elif cooling_on:
            state_txt.set_text(f"COOLING · burn {burns_done}/{cooled['n_kicks']}")
            state_txt.set_color(HORSE)
        else:
            state_txt.set_text("DEEP BIT · cooled, still reads 1")
            state_txt.set_color(GOOD)
        sub_txt.set_text(f"station-keeping burns, |dv| ≤ {cool.CAP} "
                         f"(erase kick: {memory.ERASE_KICK})")
        lo = max(0, f - 60)
        seg = xy_f[lo:f + 1]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            aa = np.linspace(0.05, 1, len(segs))
            rgb = np.array(matplotlib.colors.to_rgb(
                L4C if not cooling_on else HORSE))
            trail.set_segments(segs)
            trail.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), aa]))
            trail.set_linewidth(2.4)
        head.set_offsets([xy_f[f]])
        head.set_color(L4C if not cooling_on else HORSE)
        burn.set_offsets([xy_f[f]])
        burn.set_alpha(0.35 if recent_burn and cooling_on else 0.0)
        ax.text(0.5, 0.015, "cooling a wide bit — station-keeping burns "
                "damp the libration (rotating frame)",
                transform=ax.transAxes, color=DIM, fontsize=8.5, ha="center")
        return []

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=55, blit=False)
    out = DOCS / "cool.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
