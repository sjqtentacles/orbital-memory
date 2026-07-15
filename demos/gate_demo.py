"""The capstone, as imagery: a flyby (logic) acts on a stored bit (memory).

Two synchronized lanes, same cell, same clock — the only difference is the
bullet. Left: the aimed bullet streaks through and the moonlet is flung off
its island (bit erased). Right: no bullet, the bit keeps librating.

Writes docs/gate.gif.  Usage: python -m demos.gate_demo
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from demos.style import (DIM, DOCS, ERASE, GOOD, GROUND, HORSE, L4C, PLANET,
                         STAR, optimize_gif)
from orbital import gate, memory, nbody, rotating, theory

FRAMES = 160
T_SHOW = 22 * memory.PERIOD   # animation horizon (erase verdict is visible early)


def main():
    DOCS.mkdir(exist_ok=True)
    cell = memory.make_cell("L4", libration_deg=6.0)
    psi, achieved = gate.aim(cell)
    print(f"aimed: miss {achieved:.4f}")
    shot = gate.fire(cell, psi, t_end=T_SHOW, n_samples=FRAMES)
    ctrl = nbody.integrate([dict(b) for b in cell], T_SHOW, n_samples=FRAMES)
    verdict = gate.fire(cell, psi)  # full-length verdict for the caption
    lab, _, _ = memory.classify(gate.resonant_angle_com(verdict))
    assert lab == "erased"

    # rotating-frame views
    moon_shot = rotating.to_rotating_frame(shot, 2)
    bullet_shot = rotating.to_rotating_frame(shot, 3)
    moon_ctrl = rotating.to_rotating_frame(ctrl, 2)
    l4 = theory.lagrange_points(memory.MU)["L4"]

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6), facecolor=GROUND)
    artists = []
    for ax, title in zip(axes, ("bullet PRESENT", "bullet ABSENT")):
        ax.set_facecolor(GROUND); ax.set_aspect("equal")
        ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.5, 1.7); ax.axis("off")
        ax.set_title(title, color="w", fontsize=12, family="monospace")
        trail = LineCollection([], zorder=4); ax.add_collection(trail)
        btrail = LineCollection([], zorder=3); ax.add_collection(btrail)
        head = ax.scatter([], [], s=60, zorder=6, edgecolors="none")
        bhead = ax.scatter([], [], s=45, zorder=5, edgecolors="none",
                           color=HORSE)
        txt = ax.text(0.04, 0.96, "", transform=ax.transAxes, fontsize=11.5,
                      va="top", weight="bold", family="monospace")
        artists.append((ax, trail, btrail, head, bhead, txt))

    def lane_base(ax):
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(th), np.sin(th), color="#1a2340", lw=0.8, zorder=1)
        ax.plot(0, 0, marker="o", ms=13, color=STAR, zorder=5)
        ax.plot(1 - memory.MU, 0, marker="o", ms=6, color=PLANET, zorder=5)
        ax.plot(*l4, "x", color=L4C, ms=7, mew=1.6, zorder=5)

    def fading(lc, pts, f, color, span=45):
        lo = max(0, f - span)
        seg = pts[lo:f + 1]
        if len(seg) >= 2:
            segs = np.stack([seg[:-1], seg[1:]], axis=1)
            aa = np.linspace(0.05, 1, len(segs))
            rgb = np.array(matplotlib.colors.to_rgb(color))
            lc.set_segments(segs)
            lc.set_color(np.column_stack([np.tile(rgb, (len(segs), 1)), aa]))
            lc.set_linewidth(2.2)

    enc_i = int(np.argmin(np.linalg.norm(shot["traj"][3] - shot["traj"][2],
                                         axis=1)))

    def draw(f):
        for ax, trail, btrail, head, bhead, txt in artists:
            for ln in list(ax.lines):
                ln.remove()
            lane_base(ax)
        # lane 1: shot
        ax, trail, btrail, head, bhead, txt = artists[0]
        col = ERASE if f > enc_i else L4C
        fading(trail, moon_shot, f, col)
        fading(btrail, bullet_shot, f, HORSE, span=18)
        head.set_offsets([moon_shot[f]]); head.set_color(col)
        r_b = np.linalg.norm(shot["traj"][3][f])
        bhead.set_offsets([bullet_shot[f]])
        bhead.set_alpha(1.0 if r_b < 2.2 else 0.0)
        txt.set_text("bit = 1" if f <= enc_i else "ERASED")
        txt.set_color(L4C if f <= enc_i else ERASE)
        # lane 2: control
        ax, trail, btrail, head, bhead, txt = artists[1]
        fading(trail, moon_ctrl, f, L4C)
        head.set_offsets([moon_ctrl[f]]); head.set_color(L4C)
        txt.set_text("bit = 1"); txt.set_color(GOOD if f > enc_i else L4C)
        fig.suptitle("conditional erase — a flyby (logic) acting on a stored "
                     "bit (memory)", color="w", fontsize=12.5, y=0.045)
        return []

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=55, blit=False)
    out = DOCS / "gate.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
