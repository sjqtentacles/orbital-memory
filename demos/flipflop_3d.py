"""3D view of the gravitational memory: an inclined Trojan bit, rotating camera.

Same bit as the flat demo (L4 = 1, L5 = 0), but the test particle is inclined,
so in the rotating frame it traces a tadpole around its Lagrange point AND bobs
up and down through the orbital plane once per orbit — a genuinely 3D "wavy
banana." Renders docs/orbital_3d.gif with a slowly orbiting camera.

Usage: python -m demos.flipflop_3d
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from demos.style import DIM, DOCS, GROUND, L4C, L5C, PLANET, STAR, optimize_gif
from orbital import memory, nbody

ORBITS = 15
FRAMES = 130


def rot_frame_3d(res, body):
    t = res["t"]
    xy = res["traj"][body]
    c, s = np.cos(-t), np.sin(-t)
    return np.column_stack([c * xy[:, 0] - s * xy[:, 1],
                            s * xy[:, 0] + c * xy[:, 1], xy[:, 2]])


def main():
    DOCS.mkdir(exist_ok=True)
    runs = {}
    for st in ("L4", "L5"):
        b = memory.make_cell_3d(st, libration_deg=10.0, inclination_z=0.16)
        res = nbody.integrate(b, ORBITS * memory.PERIOD,
                              n_samples=FRAMES, rtol=1e-10, atol=1e-11)
        runs[st] = rot_frame_3d(res, 2)
    star, planet = memory.primaries()
    l4, _ = memory.lagrange_state(memory.MU, "L4")
    l5, _ = memory.lagrange_state(memory.MU, "L5")
    print(f"integrated {ORBITS} orbits for L4 and L5 (3D)")

    fig = plt.figure(figsize=(7.2, 6.4), facecolor=GROUND)
    ax = fig.add_subplot(projection="3d")
    ax.set_facecolor(GROUND)

    gx, gy = np.meshgrid(np.linspace(-1.3, 1.3, 2), np.linspace(-1.3, 1.3, 2))

    def draw(frame):
        ax.clear()
        ax.set_facecolor(GROUND)
        ax.set_axis_off()
        ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_zlim(-0.6, 0.6)
        ax.set_box_aspect((1, 1, 0.55))
        # orbital plane
        ax.plot_surface(gx, gy, np.zeros_like(gx), color="#152036",
                        alpha=0.30, zorder=0, shade=False)
        # equilateral-triangle guide lines (star-planet-L point)
        for L, col in ((l4, L4C), (l5, L5C)):
            for A in (star["pos"], planet["pos"]):
                ax.plot([A[0], L[0]], [A[1], L[1]], [0, 0], color=col,
                        lw=0.6, alpha=0.22)
        ax.plot([star["pos"][0], planet["pos"][0]],
                [star["pos"][1], planet["pos"][1]], [0, 0], color=DIM, lw=0.6, alpha=0.3)
        # bodies
        ax.scatter(*star["pos"], 0, color=STAR, s=260, edgecolors="none", zorder=5)
        ax.scatter(*planet["pos"], 0, color=PLANET, s=70, edgecolors="none", zorder=5)
        # tadpole trails + live particle, bobbing through the plane
        for st, col in (("L4", L4C), ("L5", L5C)):
            p = runs[st]
            lo = max(0, frame - 70)
            ax.plot(p[lo:frame + 1, 0], p[lo:frame + 1, 1], p[lo:frame + 1, 2],
                    color=col, lw=1.8, alpha=0.9)
            ax.scatter(p[frame, 0], p[frame, 1], p[frame, 2],
                       color=col, s=60, edgecolors="none", zorder=6)
            Lp = l4 if st == "L4" else l5
            ax.text(Lp[0], Lp[1], 0.30, f"{st} = {'1' if st == 'L4' else '0'}",
                    color=col, fontsize=11, ha="center")
        ax.text2D(0.5, 0.95, "Orbital Memory · an inclined Trojan bit in 3D",
                  transform=ax.transAxes, color="w", fontsize=12, ha="center")
        ax.text2D(0.5, 0.05,
                  "rotating frame — the bit is which corner it librates in; "
                  "it bobs through the plane each orbit",
                  transform=ax.transAxes, color=DIM, fontsize=8.5, ha="center")
        ax.view_init(elev=20 + 6 * np.sin(frame / FRAMES * 2 * np.pi),
                     azim=-70 + 200 * frame / FRAMES)
        return []

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=60, blit=False)
    out = DOCS / "orbital_3d.gif"
    anim.save(out, writer=PillowWriter(fps=20), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
