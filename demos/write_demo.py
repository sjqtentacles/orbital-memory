"""The WRITE operation, as imagery.

1. docs/timing.png — the write timing diagram (the cell's datasheet): pulse
   delay vs written value. Same blank, same pulse; timing alone selects the
   bit, with erased guard bands between the write windows.
2. docs/write.gif — one full write of bit 1, animated in the rotating frame:
   blank horseshoe -> growth pulse (the secondary visibly grows) -> pinch ->
   captured tadpole at L4.

Usage:
    python -m demos.write_demo timing   # ~2 min, 14 short runs
    python -m demos.write_demo gif      # ~1 min
"""

import pathlib
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

from orbital import memory, rotating, theory, write

DOCS = pathlib.Path(__file__).resolve().parent.parent / "docs"
GROUND = "#0a0e17"
STAR = "#ffd166"
PLANET = "#9fb0dd"
L4C = "#54d1ff"
L5C = "#ff8fa3"
HORSE = "#ffb454"
DIM = "#7c88a8"
T = memory.PERIOD


def timing():
    delays = np.arange(0, 40, 3)
    outcomes = []
    for d in delays:
        ramp = rotating.smooth_ramp(write.MU0, memory.MU,
                                    t_ramp=write.T_RAMP, t0=d * T)
        run = rotating.integrate(write.blank(), d * T + write.T_RAMP + 45 * T,
                                 n_samples=1300, mu=ramp, rtol=5e-8, atol=5e-9)
        lab, c, a = memory.classify(run["phi"][-450:])
        outcomes.append(lab)
        print(f"  delay {d:2d} orbits -> {lab}", flush=True)

    fig, ax = plt.subplots(figsize=(9.2, 3.4), facecolor=GROUND)
    ax.set_facecolor(GROUND)
    colors = {"L4": L4C, "L5": L5C, "erased": "#39445f"}
    labels = {"L4": "1", "L5": "0", "erased": "–"}
    for d, lab in zip(delays, outcomes):
        ax.bar(d, 1.0, width=2.8, color=colors[lab], edgecolor=GROUND)
        ax.text(d, 0.5, labels[lab], ha="center", va="center",
                color="#0a0e17" if lab != "erased" else DIM,
                fontsize=13, weight="bold", family="monospace")
    ax.set_xlabel("write-pulse delay  (orbits after blank preparation)",
                  color="w")
    ax.set_yticks([])
    ax.set_xlim(-2.5, 40)
    ax.tick_params(colors="#889")
    for s in ax.spines.values():
        s.set_color("#232c47")
    ax.set_title("write timing diagram — same pulse, same blank; "
                 "timing selects the bit", color="w", fontsize=12)
    ax.text(0.995, -0.34, "cyan = writes 1 (L4)   pink = writes 0 (L5)   "
            "grey = guard band (pinch mid-transit, stays erased)",
            transform=ax.transAxes, color=DIM, fontsize=8.5, ha="right")
    fig.tight_layout()
    fig.savefig(DOCS / "timing.png", dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {DOCS / 'timing.png'}")


def gif(frames=150):
    run = write.write_bit("1", n_samples=frames, rtol=1e-8, atol=1e-9)
    # frames are coarse (~1/orbit): read the post-pulse tail only
    bit, center, amp = write.read(run, window=int(frames * 0.25))
    assert bit == "1", f"write demo produced {bit}"
    xy, mus, t = run["xy"], run["mu"], run["t"]
    t_pulse0 = write.WRITE_DELAY["1"]
    t_pulse1 = t_pulse0 + write.T_RAMP

    fig, ax = plt.subplots(figsize=(6.4, 6.2), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.45); ax.axis("off")

    trail = LineCollection([], zorder=4); ax.add_collection(trail)
    head = ax.scatter([], [], s=70, zorder=6, edgecolors="none")
    planet_dot = ax.scatter([], [], s=40, color=PLANET, zorder=5,
                            edgecolors="none")
    state_txt = ax.text(0.03, 0.97, "", transform=ax.transAxes, fontsize=13.5,
                        va="top", weight="bold", family="monospace")
    mu_txt = ax.text(0.03, 0.90, "", transform=ax.transAxes, fontsize=9.5,
                     va="top", color=DIM, family="monospace")

    l4 = theory.lagrange_points(memory.MU)["L4"]
    l5 = theory.lagrange_points(memory.MU)["L5"]

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
        if now < t_pulse0:
            phase, col = "BLANK · horseshoe moonlet", HORSE
        elif now < t_pulse1:
            phase, col = "WRITE PULSE · secondary growing", "#ffe08a"
        else:
            phase, col = "WRITTEN · bit = 1  (tadpole @ L4)", L4C
        state_txt.set_text(phase); state_txt.set_color(col)
        mu_txt.set_text(f"mass ratio μ = {mus[f]:.5f}")
        # trail colored by phase
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
        ax.text(0.5, 0.015, "writing a bit by growing the planet — "
                "the horseshoe pinch (rotating frame)",
                transform=ax.transAxes, color=DIM, fontsize=8.5, ha="center")
        return []

    anim = FuncAnimation(fig, draw, frames=frames, interval=60, blit=False)
    out = DOCS / "write.gif"
    anim.save(out, writer=PillowWriter(fps=18), dpi=100,
              savefig_kwargs={"facecolor": GROUND})
    plt.close(fig)
    tool = "magick" if shutil.which("magick") else None
    if tool:
        tmp = out.with_suffix(".opt.gif")
        try:
            subprocess.run([tool, str(out), "-layers", "optimize", "-fuzz",
                            "3%", str(tmp)], check=True, capture_output=True)
            if tmp.stat().st_size < out.stat().st_size:
                tmp.replace(out)
            elif tmp.exists():
                tmp.unlink()
        except subprocess.CalledProcessError:
            tmp.exists() and tmp.unlink()
    print(f"wrote {out}  ({out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    DOCS.mkdir(exist_ok=True)
    if which in ("timing", "both"):
        timing()
    if which in ("gif", "both"):
        gif()
