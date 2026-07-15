"""The memory's phase-space anatomy, as pictures.

1. docs/landscape.png — the rotating-frame energy landscape: the effective
   potential, the zero-velocity curve at the held bit's Jacobi constant, the
   five Lagrange points, and the two tadpole orbits parked in the L4/L5 bowls.
2. docs/anatomy.png — the co-orbital orbit families that ARE the memory's
   states: tadpole at L4 (bit 1), tadpole at L5 (bit 0), a horseshoe (the
   erased band between them), and free circulation (blank medium).

Usage: python -m demos.landscape
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demos.style import (CIRC, DIM, DOCS, GROUND, HORSE, L4C, L5C, PLANET,
                         STAR)
from orbital import memory, rotating, theory, write

MU = memory.MU


def _fig_ax(title):
    fig, ax = plt.subplots(figsize=(7.6, 7.0), facecolor=GROUND)
    ax.set_facecolor(GROUND)
    ax.set_aspect("equal")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#232c47")
    ax.set_title(title, color="w", fontsize=12.5, pad=12)
    return fig, ax


def _primaries(ax):
    star, planet = memory.primaries()
    ax.plot(*star["pos"], "o", color=STAR, ms=17, zorder=6)
    ax.plot(*planet["pos"], "o", color=PLANET, ms=8, zorder=6)
    ax.annotate("primary", star["pos"], xytext=(0, -16), ha="center",
                textcoords="offset points", color=DIM, fontsize=8.5)
    ax.annotate("secondary", planet["pos"], xytext=(0, -14), ha="center",
                textcoords="offset points", color=DIM, fontsize=8.5)


def landscape():
    fig, ax = _fig_ax("the energy landscape of the bit — effective potential,"
                      "\nzero-velocity curve at the held bit's Jacobi constant")
    g = np.linspace(-1.45, 1.45, 900)
    X, Y = np.meshgrid(g, g)
    C_grid = 2 * theory.effective_potential(X, Y, mu=MU)  # C at zero velocity

    # background: shaded potential (arcsinh stretch, clipped near primaries)
    Z = np.clip(C_grid, 2.9, 3.35)
    ax.contourf(X, Y, -np.arcsinh((Z - 3.0) * 40), levels=60,
                cmap="twilight_shifted", alpha=0.75, zorder=0)
    # a few faint C contours for texture
    ax.contour(X, Y, C_grid, levels=np.linspace(2.999, 3.20, 9),
               colors="#2a3557", linewidths=0.5, zorder=1)

    # the forbidden region for the HELD bit: C_J of a real held tadpole
    run = rotating.integrate(rotating.lagrange_tadpole("L4", 6.0),
                             20 * memory.PERIOD, n_samples=400)
    C_held = float(np.mean(run["jacobi"]))
    ax.contourf(X, Y, C_grid, levels=[C_held, 1e9], colors=["#000000"],
                alpha=0.42, zorder=2)
    ax.contour(X, Y, C_grid, levels=[C_held], colors=["w"],
               linewidths=1.2, zorder=3)

    # Lagrange points + the two stored states
    pts = theory.lagrange_points(MU)
    for name, (px, py) in pts.items():
        col = {"L4": L4C, "L5": L5C}.get(name, DIM)
        ax.plot(px, py, "x", color=col, ms=8, mew=1.8, zorder=6)
        dy = 0.09 if name != "L3" else 0.11
        ax.annotate(name, (px, py + dy), color=col, ha="center",
                    fontsize=10, zorder=6)
    for st, col in (("L4", L4C), ("L5", L5C)):
        r = rotating.integrate(rotating.lagrange_tadpole(st, 6.0),
                               25 * memory.PERIOD, n_samples=1200)
        ax.plot(r["xy"][:, 0], r["xy"][:, 1], color=col, lw=1.6,
                alpha=0.95, zorder=5)
    _primaries(ax)

    ax.text(0.02, 0.03,
            "white curve: zero-velocity boundary at the stored bit's $C_J$\n"
            "shaded: classically forbidden — the moat around the memory",
            transform=ax.transAxes, color=DIM, fontsize=8.5, va="bottom")
    fig.tight_layout()
    fig.savefig(DOCS / "landscape.png", dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {DOCS / 'landscape.png'}  (C_held = {C_held:.6f}, "
          f"C_L4 = {theory.C_L4(MU):.6f})")


def anatomy():
    """Drawn at the write-time mass ratio write.MU0, where all four states
    coexist cleanly (at the full cell mu the horseshoe band sits in the
    chaotic layer — which is exactly why the write STARTS at mu0)."""
    mu0 = write.MU0
    fig, ax = _fig_ax("the four states of the medium — "
                      "tadpoles (bits), horseshoe, circulation\n"
                      r"(drawn at the write-time mass ratio $\mu_0 = 3\times10^{-4}$)")
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="#1c2440", lw=1.0, zorder=1)

    T = memory.PERIOD
    for st, col, lab in (("L4", L4C, "tadpole @ L4  ·  bit 1"),
                         ("L5", L5C, "tadpole @ L5  ·  bit 0")):
        r = rotating.integrate(rotating.lagrange_tadpole(st, 10.0, mu=mu0),
                               60 * T, n_samples=2400, mu=mu0)
        ax.plot(r["xy"][:, 0], r["xy"][:, 1], color=col, lw=1.6, label=lab,
                zorder=5)
    # erased: a horseshoe — sweeps through L3, turns around near the secondary
    # (drawn at the blank medium's exact offset write.DA)
    hs = rotating.integrate(rotating.circular_coorbital(180.0, da=write.DA),
                            150 * T, n_samples=5000, mu=mu0)
    ax.plot(hs["xy"][:, 0], hs["xy"][:, 1], color=HORSE, lw=0.9, alpha=0.9,
            label="horseshoe  ·  erased", zorder=4)
    # blank: free circulation, outside the co-orbital zone
    circ = rotating.integrate(rotating.circular_coorbital(90.0, da=0.06),
                              45 * T, n_samples=2200, mu=mu0)
    ax.plot(circ["xy"][:, 0], circ["xy"][:, 1], color=CIRC, lw=0.8,
            alpha=0.85, ls=(0, (4, 3)), label="circulation  ·  blank",
            zorder=3)

    pts = theory.lagrange_points(mu0)
    for name in ("L4", "L5"):
        px, py = pts[name]
        col = L4C if name == "L4" else L5C
        ax.plot(px, py, "x", color=col, ms=8, mew=1.8, zorder=6)
    _primaries(ax)
    ax.legend(loc="upper left", facecolor="#121a2e", edgecolor="#26304f",
              labelcolor="w", fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(DOCS / "anatomy.png", dpi=150, facecolor=GROUND)
    plt.close(fig)
    print(f"wrote {DOCS / 'anatomy.png'}")


def main():
    DOCS.mkdir(exist_ok=True)
    landscape()
    anatomy()


if __name__ == "__main__":
    main()
