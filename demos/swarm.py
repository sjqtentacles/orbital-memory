"""SWARM, as imagery: the real Greek and Trojan clouds.

Real Jupiter carries ~13,000 known Trojans in two clouds — the Greeks leading
at L4, the Trojans trailing at L5. Here a cloud of massless particles is seeded
around both islands at the real Sun-Jupiter mass ratio and advanced together in
the co-rotating field (a vectorized test-particle integrator — the full O(N^2)
engine can't carry a cloud). They librate as two lobes, exactly the swarm shape
the sky shows. Axes in AU.

Writes docs/swarm.gif.  Usage: python -m demos.swarm
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from demos.style import DIM, DOCS, GROUND, L4C, L5C, PLANET, STAR, optimize_gif
from orbital import memory, units

MU = memory.MU_SUN_JUPITER
SYS = units.SUN_JUPITER
N_PER = 500                    # particles per cloud
FRAMES = 150
T_END = 120.0                  # ~1.5 libration periods (nondim)
STEPS = 3000


def _accel(X, Y, VX, VY):
    r1 = np.hypot(X + MU, Y); r2 = np.hypot(X - (1 - MU), Y)
    ax = X + 2 * VY - (1 - MU) * (X + MU) / r1**3 - MU * (X - (1 - MU)) / r2**3
    ay = Y - 2 * VX - (1 - MU) * Y / r1**3 - MU * Y / r2**3
    return ax, ay


def _seed(point):
    """A cloud around L4/L5: spread of libration amplitude + small radial jitter."""
    rng = np.random.RandomState(0 if point == "L4" else 1)
    sy = 1.0 if point == "L4" else -1.0
    amp = rng.uniform(2, 34, N_PER) * sy            # deg along the ring
    da = rng.normal(0, 0.006, N_PER)                # radial jitter
    ang = np.radians(60.0 * sy + amp)
    r = (1.0 + da)
    X = r * np.cos(ang); Y = r * np.sin(ang)
    # corotation (rest in the rotating frame) is the tadpole turning point
    return X, Y, np.zeros(N_PER), np.zeros(N_PER)


def main():
    DOCS.mkdir(exist_ok=True)
    seeds = [_seed("L4"), _seed("L5")]
    X = np.concatenate([s[0] for s in seeds])
    Y = np.concatenate([s[1] for s in seeds])
    VX = np.concatenate([s[2] for s in seeds])
    VY = np.concatenate([s[3] for s in seeds])
    colors = [L4C] * N_PER + [L5C] * N_PER

    dt = T_END / STEPS
    keep = np.linspace(0, STEPS, FRAMES).astype(int)
    snaps = []
    for step in range(STEPS + 1):
        if step in keep:
            snaps.append((X.copy(), Y.copy()))
        ax1, ay1 = _accel(X, Y, VX, VY)
        Xh, Yh = X + 0.5 * dt * VX, Y + 0.5 * dt * VY
        VXh, VYh = VX + 0.5 * dt * ax1, VY + 0.5 * dt * ay1
        ax2, ay2 = _accel(Xh, Yh, VXh, VYh)
        X = X + dt * VXh; Y = Y + dt * VYh
        VX = VX + dt * ax2; VY = VY + dt * ay2

    fig, ax = plt.subplots(figsize=(6.6, 6.4), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    lim = SYS.au(1.35)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim * 0.92, lim * 1.05); ax.axis("off")
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(SYS.au(np.cos(th)), SYS.au(np.sin(th)), color="#1a2340", lw=0.8)
    ax.scatter([0], [0], c=STAR, s=130, zorder=5, edgecolors="none")
    ax.scatter([SYS.au(1 - MU)], [0], c=PLANET, s=40, zorder=5, edgecolors="none")
    x0, y0 = snaps[0]
    scat = ax.scatter(SYS.au(x0), SYS.au(y0), s=5, c=colors, alpha=0.7,
                      edgecolors="none")
    ax.text(SYS.au(0.5), SYS.au(0.92), "Greeks · L4", color=L4C, fontsize=10,
            ha="center")
    ax.text(SYS.au(0.5), SYS.au(-0.98), "Trojans · L5", color=L5C, fontsize=10,
            ha="center")
    ax.text(0.5, 1.0, f"~13,000 real Jupiter Trojans, {N_PER*2} simulated  "
            f"·  a = {SYS.au(1.0):.2f} AU", transform=ax.transAxes,
            color=DIM, fontsize=9.5, ha="center", va="top")

    def draw(f):
        xs, ys = snaps[f]
        scat.set_offsets(np.column_stack([SYS.au(xs), SYS.au(ys)]))
        fig.suptitle("Jupiter's Trojan swarm — two clouds of bits, "
                     "60° ahead of and behind the planet",
                     color="w", fontsize=11.5, y=0.05)
        return [scat]

    anim = FuncAnimation(fig, draw, frames=FRAMES, interval=55, blit=False)
    out = DOCS / "swarm.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    optimize_gif(out)
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
