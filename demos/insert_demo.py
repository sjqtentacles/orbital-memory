"""WRITE, as imagery: a bit placed by orbit insertion.

One rotating-frame panel: the body coasts in on a co-orbital transfer, an
insertion burn (flash) circularizes it onto the L4 tadpole, and it settles
into a librating bit. The burn is a real delta-v, reported in m/s at the real
Sun-Jupiter scale — the honest cost of placing a bit, no grown planet.

Writes docs/insert.gif.  Usage: python -m demos.insert_demo
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import (DIM, DOCS, GOOD, GROUND, HORSE, L4C, PLANET, STAR,
                         optimize_gif)
from orbital import memory, nbody, rotating, units, write

APPROACH = 40                 # approach frames
LIBRATE = 130                 # libration frames
T_APPROACH = 6 * memory.PERIOD
T_LIBRATE = 30 * memory.PERIOD


def main():
    DOCS.mkdir(exist_ok=True)
    mu = memory.MU
    cell, dv = write.insert("1", amplitude_deg=8.0, mu=mu)
    dv_mps = units.SUN_JUPITER.mps(dv)

    # post-burn: the tadpole (forward run)
    post = nbody.integrate(cell, T_LIBRATE, n_samples=LIBRATE)
    moon_post = rotating.to_rotating_frame(post, 2)

    # pre-burn: same position, slower transfer-apoapsis velocity, run BACKWARD
    P = np.array(cell[2]["pos"])
    that = np.array([-P[1], P[0]]) / np.linalg.norm(P)
    v_apo = np.linalg.norm(cell[2]["vel"]) - dv
    pre_cell = [dict(cell[0]), dict(cell[1]),
                {"m": 0.0, "pos": P.tolist(), "vel": (v_apo * that).tolist()}]
    back = nbody.integrate(pre_cell, -T_APPROACH, n_samples=APPROACH, t0=0.0)
    moon_pre = rotating.to_rotating_frame(back, 2)[::-1]   # forward playback

    l4 = np.array([0.5 - mu, np.sqrt(3) / 2])
    fig, ax = plt.subplots(figsize=(6.6, 6.4), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.4, 1.6); ax.axis("off")
    trail = LineCollection([], zorder=4); ax.add_collection(trail)
    head = ax.scatter([], [], s=70, zorder=6, edgecolors="none")
    burn = ax.scatter([], [], s=0, zorder=7, color=GOOD, edgecolors="none")
    txt = ax.text(0.04, 0.96, "", transform=ax.transAxes, fontsize=12,
                  va="top", weight="bold", family="monospace")
    sub = ax.text(0.04, 0.05, "", transform=ax.transAxes, fontsize=10,
                  va="bottom", color=DIM, family="monospace")

    def base():
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(th), np.sin(th), color="#1a2340", lw=0.8, zorder=1)
        ax.plot(0, 0, marker="o", ms=15, color=STAR, zorder=5)
        ax.plot(1 - mu, 0, marker="o", ms=7, color=PLANET, zorder=5)
        ax.plot(*l4, "x", color=L4C, ms=8, mew=1.6, zorder=5)

    TOTAL = APPROACH + LIBRATE

    def frame_xy(f):
        return moon_pre[f] if f < APPROACH else moon_post[f - APPROACH]

    def draw(f):
        for ln in list(ax.lines):
            ln.remove()
        base()
        pts = np.vstack([moon_pre[:min(f, APPROACH)],
                         moon_post[:max(0, f - APPROACH)]]) if f else moon_pre[:1]
        lo = max(0, len(pts) - 55)
        seg = pts[lo:]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            aa = np.linspace(0.05, 1, len(segs))
            col = HORSE if f < APPROACH else L4C
            rgb = np.array(matplotlib.colors.to_rgb(col))
            trail.set_segments(segs)
            trail.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), aa]))
            trail.set_linewidth(2.2)
        p = frame_xy(f)
        head.set_offsets([p]); head.set_color(HORSE if f < APPROACH else L4C)
        # burn flash near the insertion frame
        d = f - APPROACH
        burn.set_sizes([max(0, 900 - 240 * abs(d))] if abs(d) < 4 else [0])
        burn.set_offsets([moon_post[0]])
        if f < APPROACH:
            txt.set_text("APPROACH"); txt.set_color(HORSE)
            sub.set_text("co-orbital transfer toward L4")
        elif abs(d) < 4:
            txt.set_text("INSERTION BURN"); txt.set_color(GOOD)
            sub.set_text(f"Δv = {dv_mps:.0f} m/s")
        else:
            txt.set_text("bit = 1  (L4)"); txt.set_color(L4C)
            sub.set_text(f"written · insertion cost {dv_mps:.0f} m/s")
        fig.suptitle("WRITE by orbit insertion — placing a bit, no grown planet",
                     color="w", fontsize=12.5, y=0.04)
        return []

    anim = FuncAnimation(fig, draw, frames=TOTAL, interval=55, blit=False)
    out = DOCS / "insert.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB, Δv={dv_mps:.0f} m/s)")


if __name__ == "__main__":
    main()
