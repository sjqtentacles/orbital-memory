"""Integrator invariants — the numerical bedrock every memory claim rests on.

Energy/momentum conservation, a closed-form Kepler check (a moon on a circular
orbit), determinism, and 2D/3D consistency.
"""

import numpy as np
import pytest

from orbital import memory, nbody


def _momentum(res):
    m = res["masses"][:, None]
    return (m * res["vel"][:, -1, :]).sum(axis=0)


def _com(res, frame):
    m = res["masses"]
    pos = res["traj"][:, frame, :]
    return (m[:, None] * pos).sum(axis=0) / m.sum()


class TestKeplerMoon:
    """A light moon on a circular orbit around a heavy body: the mental model
    the project is named for. Tests the integrator against Kepler's laws."""

    def _moon(self, a=1.0, M=1.0):
        v = np.sqrt(nbody.G * M / a)  # circular speed
        return [{"m": M, "pos": [0.0, 0.0], "vel": [0.0, 0.0]},
                {"m": 1e-9, "pos": [a, 0.0], "vel": [0.0, v]}]

    def test_circular_orbit_stays_circular(self):
        res = nbody.integrate(self._moon(a=1.0), 6 * np.pi, n_samples=400)
        r = np.linalg.norm(res["traj"][1] - res["traj"][0], axis=1)
        assert np.std(r) < 1e-4  # radius constant

    def test_period_matches_keplers_third_law(self):
        a, M = 1.3, 1.0
        T = 2 * np.pi * np.sqrt(a ** 3 / (nbody.G * M))
        res = nbody.integrate(self._moon(a, M), T, n_samples=2000)
        start = np.array([a, 0.0])
        end = res["traj"][1, -1] - res["traj"][0, -1]
        assert np.linalg.norm(end - start) < 5e-3  # back to start after one period

    def test_energy_conserved(self):
        res = nbody.integrate(self._moon(), 20 * np.pi, n_samples=1000)
        assert res["energy_drift"] < 1e-10


class TestConservation:
    def test_energy_2d_cell(self):
        res = nbody.integrate(memory.make_cell("L4"), 40 * memory.PERIOD,
                              n_samples=2000)
        assert res["energy_drift"] < 1e-9

    def test_energy_3d_cell(self):
        res = nbody.integrate(memory.make_cell_3d("L4"), 40 * memory.PERIOD,
                              n_samples=2000)
        assert res["energy_drift"] < 1e-9

    def test_momentum_conserved_and_zero(self):
        res = nbody.integrate(memory.make_cell("L4"), 30 * memory.PERIOD,
                              n_samples=1500)
        assert np.allclose(_momentum(res), [0.0, 0.0], atol=1e-10)

    def test_barycenter_pinned_at_origin(self):
        """resonant_angle uses the origin as the barycenter — it must stay there."""
        res = nbody.integrate(memory.make_cell("L4"), 30 * memory.PERIOD,
                              n_samples=1500)
        for f in (0, 750, -1):
            assert np.linalg.norm(_com(res, f)) < 1e-10


class TestMasslessParticle:
    def test_no_back_reaction_on_primaries(self):
        """A m=0 test particle must not perturb the star+planet two-body orbit."""
        with_p = nbody.integrate(memory.make_cell("L4"), 20 * memory.PERIOD,
                                 n_samples=800)
        star, planet = memory.primaries()
        without = nbody.integrate([star, planet], 20 * memory.PERIOD, n_samples=800)
        assert np.allclose(with_p["traj"][:2], without["traj"], atol=1e-9)


class TestDimensionAgnostic:
    def test_2d_equals_3d_in_plane(self):
        """A 3D run with z=0, vz=0 must reproduce the 2D run exactly (xy)."""
        b2 = memory.make_cell("L4", libration_deg=6.0)
        b3 = [dict(x, pos=x["pos"] + [0.0], vel=x["vel"] + [0.0]) for x in b2]
        r2 = nbody.integrate(b2, 15 * memory.PERIOD, n_samples=600)
        r3 = nbody.integrate(b3, 15 * memory.PERIOD, n_samples=600)
        assert r2["dim"] == 2 and r3["dim"] == 3
        assert np.allclose(r2["traj"], r3["traj"][:, :, :2], atol=1e-8)
        assert np.allclose(r3["traj"][:, :, 2], 0.0, atol=1e-12)


class TestDeterminism:
    def test_repeatable(self):
        a = nbody.integrate(memory.make_cell("L4"), 12 * memory.PERIOD, n_samples=500)
        b = nbody.integrate(memory.make_cell("L4"), 12 * memory.PERIOD, n_samples=500)
        assert np.array_equal(a["traj"], b["traj"])


class TestTimeReversal:
    def test_conservative_run_retraces(self):
        """No drag: flip velocities at the end, integrate back, land on start."""
        cell = memory.make_cell("L4", libration_deg=6.0)
        fwd = nbody.integrate(cell, 8 * memory.PERIOD, n_samples=400)
        back = memory.state_to_bodies(fwd)
        for b in back:
            b["vel"] = [-c for c in b["vel"]]
        rev = nbody.integrate(back, 8 * memory.PERIOD, n_samples=400)
        for i, b in enumerate(cell):
            assert np.allclose(rev["traj"][i, -1], b["pos"], atol=1e-5)
