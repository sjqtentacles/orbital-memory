"""The memory bit: geometry, equilibrium, retention, readout, and noise margin.

Where possible the simulation is checked against closed-form Lagrange theory
(equilateral geometry, tadpole libration period), not just against itself.
"""

import numpy as np
import pytest

from orbital import memory, nbody, rotating


class TestPrimaries:
    def test_zero_net_momentum(self):
        star, planet = memory.primaries()
        p = star["m"] * np.array(star["vel"]) + planet["m"] * np.array(planet["vel"])
        assert np.allclose(p, [0.0, 0.0])

    def test_barycenter_at_origin(self):
        star, planet = memory.primaries()
        com = star["m"] * np.array(star["pos"]) + planet["m"] * np.array(planet["pos"])
        assert np.allclose(com, [0.0, 0.0])


class TestLagrangeGeometry:
    def test_equilateral_triangle(self):
        star, planet = memory.primaries()
        for pt in ("L4", "L5"):
            L, _ = memory.lagrange_state(memory.MU, pt)
            d_star = np.linalg.norm(L - np.array(star["pos"]))
            d_planet = np.linalg.norm(L - np.array(planet["pos"]))
            d_sp = np.linalg.norm(np.array(star["pos"]) - np.array(planet["pos"]))
            assert d_star == pytest.approx(1.0, abs=1e-12)
            assert d_planet == pytest.approx(1.0, abs=1e-12)
            assert d_sp == pytest.approx(1.0, abs=1e-12)

    def test_L4_leads_L5_trails(self):
        l4, _ = memory.lagrange_state(memory.MU, "L4")
        l5, _ = memory.lagrange_state(memory.MU, "L5")
        assert l4[1] > 0 and l5[1] < 0  # leading (+y) vs trailing (-y)

    def test_exact_point_is_an_equilibrium(self):
        """Seeded with libration_deg=0, the particle should sit still at L4 in
        the rotating frame (a true equilibrium of the restricted problem)."""
        cell = memory.make_cell("L4", libration_deg=0.0)
        res = nbody.integrate(cell, 20 * memory.PERIOD, n_samples=800)
        rp = rotating.to_rotating_frame(res, 2)
        drift = np.linalg.norm(rp - rp[0], axis=1).max()
        assert drift < 1e-3


class TestReadAndHold:
    @pytest.mark.parametrize("state,sign,bit", [("L4", +1, "L4"), ("L5", -1, "L5")])
    def test_holds_and_reads(self, state, sign, bit):
        res = memory.hold(state, periods=80, libration_deg=6.0)
        label, center, amp = memory.classify(res["phi"])
        assert label == bit
        assert np.sign(center) == sign
        assert abs(abs(center) - 60) < 6  # librates about +/-60 deg
        assert res["energy_drift"] < 1e-9

    def test_states_are_mirror_symmetric(self):
        l4 = memory.classify(memory.hold("L4", periods=40)["phi"])
        l5 = memory.classify(memory.hold("L5", periods=40)["phi"])
        assert l4[1] == pytest.approx(-l5[1], abs=2.0)   # opposite centers
        assert l4[2] == pytest.approx(l5[2], abs=2.0)    # similar amplitude

    @pytest.mark.slow
    def test_retention_amplitude_bounded(self):
        """Nonvolatility: over 300 orbits the libration must not grow (no secular
        drift toward the separatrix)."""
        res = memory.hold("L4", periods=300, libration_deg=6.0, n_per_period=30)
        phi = res["phi"]
        early = phi[:len(phi) // 6]
        late = phi[-len(phi) // 6:]
        amp_e = (early.max() - early.min()) / 2
        amp_l = (late.max() - late.min()) / 2
        assert amp_l < amp_e + 3.0
        assert memory.classify(phi)[0] == "L4"


class TestAgainstTheory:
    def test_libration_period_matches_linear_theory(self):
        """Measured tadpole period should match 2*pi/sqrt(27/4 mu) to a few %."""
        res = memory.hold("L4", periods=60, libration_deg=5.0, n_per_period=80)
        phi = res["phi"] - 60.0                       # center on the L4 point
        t = res["t"]
        crossings = t[np.where((phi[:-1] < 0) & (phi[1:] >= 0))[0]]
        assert len(crossings) >= 3
        measured = np.mean(np.diff(crossings))        # one libration period
        theory = memory.libration_period(memory.MU)
        assert measured == pytest.approx(theory, rel=0.08)

    @pytest.mark.parametrize("mu", [0.001, 0.01])
    def test_stable_across_mass_ratios(self, mu):
        """L4 holds for any mu below the Gascheau/Routh limit (~0.0385)."""
        res = memory.hold("L4", periods=40, mu=mu, libration_deg=5.0)
        assert memory.classify(res["phi"])[0] == "L4"


class TestNoiseMargin:
    def _kick_outcome(self, frac):
        res = nbody.integrate(memory.kicked_cell(frac), 60 * memory.PERIOD,
                              n_samples=3000)
        return memory.classify(memory.resonant_angle(res))[0]

    def test_small_kick_preserves_bit(self):
        assert self._kick_outcome(0.02) == "L4"

    def test_kick_at_named_threshold_erases_bit(self):
        """ERASE_KICK is the promoted noise-margin constant: a kick at the
        threshold must actually erase (it is defined as the measured minimum
        erasing kick, in fractions of orbital speed)."""
        assert memory.ERASE_KICK == pytest.approx(0.035)
        assert self._kick_outcome(memory.ERASE_KICK) == "erased"

    def test_large_kick_erases_bit(self):
        assert self._kick_outcome(0.05) == "erased"


class TestDissipation:
    def test_drag_destabilizes_L4(self):
        """The claim the README rests on, committed as a test: L4/L5 are
        Coriolis-stabilized, so velocity drag REMOVES the stabilization —
        a tadpole under weak corotation drag is expelled, while the drag-free
        control holds. Memory here is topological, not dissipative.
        (Short horizon + loose tolerances on purpose: within ~10 orbits the
        dragged moonlet has already spiraled from r=1 to r~0.03 — going
        further just grinds the integrator inside the star's well, and
        precision is irrelevant to the verdict.)"""
        cell = memory.make_cell("L4", libration_deg=15.0)
        held = nbody.integrate(cell, 10 * memory.PERIOD, n_samples=500,
                               rtol=1e-8, atol=1e-9)
        assert memory.classify(memory.resonant_angle(held))[0] == "L4"
        dragged = nbody.integrate(cell, 10 * memory.PERIOD, n_samples=500,
                                  rtol=1e-8, atol=1e-9, drag=(0.05, [2]))
        res = memory.classify(memory.resonant_angle(dragged))
        assert res[0] == "erased"
        # and it truly left the co-orbital ring (spiraled inward)
        r = np.linalg.norm(dragged["traj"][2], axis=1)
        assert r.min() < 0.5


class TestThreeD:
    def test_inclined_cell_holds_bit_and_bobs(self):
        cell = memory.make_cell_3d("L4", libration_deg=6.0, inclination_z=0.16)
        res = nbody.integrate(cell, 40 * memory.PERIOD, n_samples=2400)
        assert memory.classify(memory.resonant_angle(res))[0] == "L4"
        z = res["traj"][2, :, 2]
        assert z.max() > 0.1 and z.min() < -0.1          # genuine out-of-plane bob
        crossings = np.sum(np.diff(np.sign(z)) != 0)     # ~2 per orbit
        assert crossings >= 40

    def test_kick_is_dimension_safe(self):
        """A 2-vector kick on a 3D cell must leave vz intact (regression)."""
        cell = memory.make_cell_3d("L4")
        vz0 = cell[2]["vel"][2]
        kicked = memory.kick(cell, [0.01, 0.0])
        assert len(kicked[2]["vel"]) == 3
        assert kicked[2]["vel"][2] == vz0


class TestClassify:
    def test_circulating_orbit_reads_erased(self):
        # a resonant angle that winds all the way around = not trapped
        phi = np.linspace(0, 10 * 360.0, 500)
        wrapped = (phi + 180) % 360 - 180
        assert memory.classify(wrapped)[0] == "erased"

    def test_tadpole_series_reads_its_island(self):
        t = np.linspace(0, 4 * np.pi, 400)
        l4 = 60 + 20 * np.sin(t)
        l5 = -60 + 20 * np.sin(t)
        assert memory.classify(l4)[0] == "L4"
        assert memory.classify(l5)[0] == "L5"

    def test_wide_tadpole_still_reads_its_island(self):
        """Regression for the audit's worry: a pinch-written tadpole librates
        ±70° around a center near 95°, bottoming ~25° — it never crosses
        conjunction (0°) or L3 (±180°), so it MUST read as a bit, not erased.
        (Real-trajectory coverage of the same fact lives in test_write.)"""
        t = np.linspace(0, 6 * np.pi, 600)
        wide_l4 = 95 + 70 * np.sin(t)          # range [25, 165]
        wide_l5 = -(95 + 70 * np.sin(t))
        assert memory.classify(wide_l4)[0] == "L4"
        assert memory.classify(wide_l5)[0] == "L5"

    def test_horseshoe_series_reads_erased(self):
        """Sweeps through L3 (±180) but never conjunction — still no bit."""
        t = np.linspace(0, 4 * np.pi, 800)
        horseshoe = 180.0 * np.sin(t / 2 + 0.3)
        horseshoe = (horseshoe + 180) % 360 - 180
        assert memory.classify(horseshoe)[0] == "erased"

    def test_parked_bit_is_still_a_bit(self):
        """A constant series is a moonlet sitting AT the point (a valid, if
        idealized, stored bit) — classified by side, not rejected."""
        assert memory.classify(np.full(50, 60.0))[0] == "L4"
        assert memory.classify(np.full(50, -60.0))[0] == "L5"

    def test_too_short_series_is_not_evidence(self):
        # a handful of samples can't distinguish libration from anything
        assert memory.classify(np.array([60.0]))[0] == "erased"
        assert memory.classify(np.full(5, 60.0))[0] == "erased"

    def test_empty_series_raises(self):
        with pytest.raises(ValueError):
            memory.classify(np.array([]))
