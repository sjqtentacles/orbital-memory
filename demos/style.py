"""Shared visual system for every demo: one palette, one dark-figure factory,
one GIF-optimize pass. Keeps the imagery consistent and the code DRY."""

import pathlib
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DOCS = pathlib.Path(__file__).resolve().parent.parent / "docs"

GROUND = "#0a0e17"
PANEL = "#121a2e"
GRID = "#1c2440"
STAR = "#ffd166"
PLANET = "#9fb0dd"
L4C = "#54d1ff"
L5C = "#ff8fa3"
HORSE = "#ffb454"
CIRC = "#8b93b8"
ERASE = "#ff6b5b"
GOOD = "#58e0a8"
DIM = "#7c88a8"
TEXT_W = "w"


def dark_fig(figsize, title=None, aspect_equal=True):
    """A GROUND-styled figure/axes pair with the house look."""
    fig, ax = plt.subplots(figsize=figsize, facecolor=GROUND)
    ax.set_facecolor(GROUND)
    if aspect_equal:
        ax.set_aspect("equal")
    if title:
        ax.set_title(title, color=TEXT_W, fontsize=12.5, pad=12)
    return fig, ax


def optimize_gif(path):
    """Shrink a GIF in place with ImageMagick if available (v7 `magick` or
    v6 `convert`); silently no-op otherwise."""
    tool = shutil.which("magick") or shutil.which("convert")
    if not tool:
        return
    tmp = path.with_suffix(".opt.gif")
    try:
        subprocess.run([tool, str(path), "-layers", "optimize", "-fuzz", "3%",
                        str(tmp)], check=True, capture_output=True)
        if tmp.exists() and tmp.stat().st_size < path.stat().st_size:
            tmp.replace(path)
        elif tmp.exists():
            tmp.unlink()
    except subprocess.CalledProcessError:
        if tmp.exists():
            tmp.unlink()
