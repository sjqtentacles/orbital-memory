"""The Jacobi constant and Lagrange points — closed-form theory vs simulation.

These are the tests that turn the memory's claims (it holds; a kick past a
margin erases it) into statements checked against the conserved quantity of
the problem, not just against the simulation itself.
"""

import numpy as np
import pytest

from orbital import memory, nbody, theory


class TestLagrangePoints:
    def test_triangular_points_are_equilateral(self):
        pts = theory.lagrange_points(memory.MU)
        star, planet = memory.primaries()
        for name in ("L4", "L5"):
            L = np.array(pts[name])
            assert np.linalg.norm(L - star["pos"]) == pytest.approx(1.0, abs=1e-9)
            assert np.linalg.norm(L - planet["pos"]) == pytest.approx(1.0, abs=1e-9)

    def test_L3_is_opposite_the_secondary(self):
        pts = theory.lagrange_points(memory.MU)
        assert pts["L3"][0] < -1.0 and pts["L3"][1] == 0.0

    def test_L3_sits_at_a_potential_extremum(self):
        """dΩ/dx = 0 at L3 (it's a genuine equilibrium of the rotating frame)."""
        mu = memory.MU
        xL3 = theory.lagrange_points(mu)["L3"][0]
        h = 1e-6
        d = (theory.effective_potential(xL3 + h, 0.0)
             - theory.effective_potential(xL3 - h, 0.0)) / (2 * h)
        assert abs(d) < 1e-4


class TestJacobiConstant:
    def test_conserved_along_held_bit(self):
        res = memory.hold("L4", periods=40, libration_deg=6.0)
        C = theory.jacobi_of(res)
        assert (C.max() - C.min()) < 1e-9   # a strong integrator invariant

    def test_conserved_along_erased_orbit(self):
        cell = memory.make_cell("L4", libration_deg=2.0)
        pre = nbody.integrate(cell, 5 * memory.PERIOD, n_samples=200)
        b = memory.state_to_bodies(pre)
        v = np.array(b[2]["vel"]); u = v / np.linalg.norm(v)
        er = nbody.integrate(memory.kick(b, (u * np.linalg.norm(v) * 0.05).tolist()),
                             40 * memory.PERIOD, n_samples=2000)
        C = theory.jacobi_of(er)
        assert (C.max() - C.min()) < 1e-7

    def test_held_bit_matches_analytic_C_L4(self):
        res = memory.hold("L4", periods=40, libration_deg=6.0)
        assert theory.jacobi_of(res).mean() == pytest.approx(theory.C_L4(memory.MU),
                                                             abs=2e-4)

    def test_triangular_points_share_a_jacobi_constant(self):
        l4 = theory.jacobi_of(memory.hold("L4", periods=40)).mean()
        l5 = theory.jacobi_of(memory.hold("L5", periods=40)).mean()
        assert l4 == pytest.approx(l5, abs=1e-3)
        assert l4 == pytest.approx(theory.C_L4(memory.MU), abs=2e-4)

    def test_kick_lowers_jacobi_below_the_held_value(self):
        """The noise margin, in the conserved quantity: erasing the bit means
        adding energy, which lowers C_J below the held (≈C_L4) value."""
        held = theory.jacobi_of(memory.hold("L4", periods=40)).mean()
        cell = memory.make_cell("L4", libration_deg=2.0)
        pre = nbody.integrate(cell, 5 * memory.PERIOD, n_samples=200)
        b = memory.state_to_bodies(pre)
        v = np.array(b[2]["vel"]); u = v / np.linalg.norm(v)
        er = nbody.integrate(memory.kick(b, (u * np.linalg.norm(v) * 0.05).tolist()),
                             40 * memory.PERIOD, n_samples=2000)
        erased = theory.jacobi_of(er).mean()
        assert erased < held - 1e-3
