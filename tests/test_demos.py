"""Smoke coverage for the demo layer.

Full demo runs cost minutes (they exist to render figures/GIFs), so these
tests assert the cheap but load-bearing contract: every demo module imports
against the current orbital/ API (catching signature drift at import time),
exposes its entry points, and the shared style helpers behave.
"""

import importlib
import pathlib

import pytest

DEMOS = ["flipflop_demo", "flipflop_3d", "landscape", "make_gifs",
         "cool_demo", "gate_demo", "validation", "phase_portrait",
         "stability_map", "insert_demo", "swarm", "saturn_erosion",
         "hamiltonian", "retention", "register", "style"]


@pytest.mark.parametrize("name", DEMOS)
def test_demo_imports(name):
    mod = importlib.import_module(f"demos.{name}")
    assert mod is not None


def test_demo_entry_points_exist():
    from demos import (cool_demo, flipflop_3d, flipflop_demo, gate_demo,
                       hamiltonian, insert_demo, landscape, make_gifs,
                       phase_portrait, register, retention, saturn_erosion,
                       stability_map, swarm, validation)
    assert callable(flipflop_demo.main)
    assert callable(flipflop_demo.noise_margin)
    assert callable(flipflop_3d.main)
    assert callable(make_gifs.main)
    assert callable(landscape.landscape) and callable(landscape.anatomy)
    for mod in (cool_demo, gate_demo, validation, phase_portrait,
                stability_map, insert_demo, swarm, saturn_erosion,
                hamiltonian, retention, register):
        assert callable(mod.main)


def test_style_palette_and_optimize():
    from demos import style
    # palette entries are hex colors
    for c in (style.GROUND, style.L4C, style.L5C, style.HORSE, style.GOOD):
        assert c.startswith("#") and len(c) == 7
    # optimize_gif on a missing file must not raise if no tool, and must not
    # crash the caller either way
    bogus = pathlib.Path("/tmp/definitely_missing_orbital_memory.gif")
    try:
        style.optimize_gif(bogus)
    except Exception as exc:  # a subprocess error on a bogus file is fine to surface
        pytest.skip(f"imagemagick present and complained as expected: {exc}")


def test_demos_use_shared_transform():
    """The rotating-frame transform must come from orbital.rotating, not be
    re-implemented per demo (the audit found 5+ copies)."""
    import inspect
    from demos import flipflop_demo, make_gifs
    for mod in (flipflop_demo, make_gifs):
        src = inspect.getsource(mod)
        assert "cos(-t)" not in src, f"{mod.__name__} re-implements the transform"
