"""REGISTER, as imagery: more than one bit.

A nibble written into four co-orbital cells at spaced radii — each host a light
secondary, each carrying a Trojan at its L4 (cyan, bit 1) or L5 (pink, bit 0).
The inertial snapshot shows the concentric host rings with their Trojan clumps
leading or trailing; the register reads back the word it was written. Capacity
is bounded by mutual stability (hosts closer than ~1.7x radius ratio go chaotic).

Writes docs/register.png.  Usage: python -m demos.register
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import DIM, DOCS, GOOD, GRID, GROUND, L4C, L5C, PLANET, STAR, TEXT_W
from orbital import memory, register

BITS = "1101"


def main():
    DOCS.mkdir(exist_ok=True)
    run = register.hold_register(BITS)
    got = register.read_register(run)
    rs = register.radii(len(BITS))

    fig, ax = plt.subplots(figsize=(7.2, 7.2), facecolor=GROUND)
    ax.set_facecolor(GROUND); ax.set_aspect("equal")
    lim = rs[-1] * 1.15
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.axis("off")

    th = np.linspace(0, 2 * np.pi, 300)
    ax.scatter([0], [0], c=STAR, s=180, zorder=6, edgecolors="none")

    # show the last ~1 orbit of each host so the Trojan clump is a compact arc
    t = run["t"]
    for i, (r, b) in enumerate(zip(rs, BITS)):
        col = L4C if b == "1" else L5C
        ax.plot(r * np.cos(th), r * np.sin(th), color=GRID, lw=0.8, zorder=1)
        orbit = 2 * np.pi * r ** 1.5
        idx = int(np.searchsorted(t, max(t[-1] - orbit, 0.0)))
        planet = run["traj"][1 + 2 * i][idx:]
        troj = run["traj"][2 + 2 * i][idx:]
        ax.plot(planet[:, 0], planet[:, 1], color=PLANET, lw=1.0, alpha=0.5, zorder=2)
        ax.plot(troj[:, 0], troj[:, 1], color=col, lw=2.4, alpha=0.9, zorder=4)
        ax.scatter([planet[-1, 0]], [planet[-1, 1]], c=PLANET, s=40, zorder=5,
                   edgecolors="none")
        ax.scatter([troj[-1, 0]], [troj[-1, 1]], c=col, s=55, zorder=5,
                   edgecolors="none")
        # bit label just outside the ring
        ax.text(0, r + lim * 0.03, b, color=col, fontsize=13, ha="center",
                va="bottom", weight="bold", family="monospace")

    ax.text(0.5, 1.02, f"wrote {BITS}   →   read {got}",
            transform=ax.transAxes, color=(GOOD if got == BITS else "#ff6b5b"),
            fontsize=14, ha="center", va="bottom", weight="bold",
            family="monospace")
    ax.text(0.5, -0.02, "cyan = L4 (bit 1, leading)   ·   pink = L5 (bit 0, trailing)"
            "   ·   capacity bounded by host stability (spacing ≳ 1.7)",
            transform=ax.transAxes, color=DIM, fontsize=9.5, ha="center", va="top")
    fig.suptitle("a register — four co-orbital bits at spaced radii",
                 color=TEXT_W, fontsize=13, y=0.97)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    out = DOCS / "register.png"
    fig.savefig(out, dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {out}  ({BITS} -> {got})")


if __name__ == "__main__":
    main()
