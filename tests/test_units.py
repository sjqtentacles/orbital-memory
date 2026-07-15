"""Real physical units — the layer that turns the canonical (G=1, a=1, n=1)
simulation into metres, kilograms, seconds, and validates it against the real
solar system.

The engine stays nondimensional; a ``System`` only carries the three scale
factors (length = a, mass = M_total, time = 1/n) and converts results out.
The point of these tests is that those scales, fed the *real* Sun-Jupiter
numbers, reproduce Jupiter's real orbital period and speed.
"""

import numpy as np
import pytest

from orbital import memory, units


class TestScales:
    def test_sun_jupiter_mu_matches_the_preset(self):
        assert units.SUN_JUPITER.mu == pytest.approx(memory.MU_SUN_JUPITER, rel=2e-3)

    def test_jupiter_orbital_period_is_11_86_years(self):
        assert units.SUN_JUPITER.period_years == pytest.approx(11.86, rel=1e-2)

    def test_jupiter_orbital_speed_is_13_km_per_s(self):
        # the nondimensional circular speed of the primaries is 1 -> vel_scale
        assert units.SUN_JUPITER.kmps(1.0) == pytest.approx(13.06, rel=1e-2)

    def test_length_scale_is_the_semimajor_axis_in_au(self):
        assert units.SUN_JUPITER.au(1.0) == pytest.approx(5.2044, rel=1e-3)

    def test_saturn_moon_ratios_are_order_of_magnitude_real(self):
        assert units.SATURN_TETHYS.mu == pytest.approx(memory.MU_SATURN_TETHYS, rel=0.1)
        assert units.SATURN_DIONE.mu == pytest.approx(memory.MU_SATURN_DIONE, rel=0.1)


class TestConversionRoundTrips:
    def test_time_round_trips(self):
        s = units.SUN_JUPITER
        for t_nd in (1.0, memory.PERIOD, 137.0):
            assert s.nd_time(s.years(t_nd)) == pytest.approx(t_nd, rel=1e-12)

    def test_length_round_trips(self):
        s = units.SUN_JUPITER
        assert s.nd_length(s.au(0.37)) == pytest.approx(0.37, rel=1e-12)

    def test_velocity_round_trips(self):
        s = units.SUN_JUPITER
        assert s.nd_vel(s.mps(0.021)) == pytest.approx(0.021, rel=1e-12)

    def test_one_nondim_period_is_one_jupiter_year(self):
        s = units.SUN_JUPITER
        assert s.years(memory.PERIOD) == pytest.approx(s.period_years, rel=1e-12)


class TestLibrationInYears:
    def test_trojan_libration_period_lands_near_148_years(self):
        """Linear theory: T_lib = P / sqrt(27/4 mu). At the real Sun-Jupiter
        ratio that is ~148 yr — the observed Jupiter-Trojan libration period.
        This is the closed-form check; test_validation.py runs the sim."""
        mu = memory.MU_SUN_JUPITER
        t_lib_nd = memory.libration_period(mu)
        assert units.SUN_JUPITER.years(t_lib_nd) == pytest.approx(148.0, rel=0.05)
