"""Does this actually work? — the simulation against the real sky.

The cell is run at the REAL Sun-Jupiter mass ratio and the result converted to
physical years with orbital.units. The headline assertion: the simulated bit
librates with a ~148-year period, which is the period actually observed for
Jupiter's Trojan asteroids. Nothing here is fit to that number — it falls out
of Newtonian gravity at the real mass ratio.
"""

import numpy as np
import pytest

from orbital import memory, nbody, theory, units


def _run_real_trojan(libration_deg=10.0, periods=45, n_samples=6000, mu=None):
    mu = memory.MU_SUN_JUPITER if mu is None else mu
    res = nbody.integrate(memory.make_cell("L4", mu=mu, libration_deg=libration_deg),
                          periods * memory.PERIOD, n_samples=n_samples)
    res["phi"] = memory.resonant_angle(res)
    return res


class TestJupiterTrojanPeriod:
    def test_simulated_libration_period_matches_observed_148_years(self):
        """The money test: sim libration period at the real ratio, in years."""
        res = _run_real_trojan()
        t_lib_nd = memory.measured_libration_period(res["t"], res["phi"])
        years = units.SUN_JUPITER.years(t_lib_nd)
        assert years == pytest.approx(148.0, rel=0.05)

    def test_simulation_agrees_with_linear_theory(self):
        res = _run_real_trojan()
        measured = memory.measured_libration_period(res["t"], res["phi"])
        assert measured == pytest.approx(memory.libration_period(memory.MU_SUN_JUPITER),
                                         rel=0.02)

    def test_the_bit_reads_and_is_centered_on_L4(self):
        res = _run_real_trojan()
        label, center, amp = memory.classify(res["phi"])
        assert label == "L4"
        assert center == pytest.approx(60.0, abs=3.0)


class TestStabilityAtTheRealRatio:
    def test_bit_does_not_degrade(self):
        """A memory must not drift. The Jacobi constant of the test particle has
        no SECULAR trend (its small ripple just breathes with the libration):
        the first- and last-third means agree to well below the ripple."""
        res = _run_real_trojan(periods=45)
        C = theory.jacobi_of(res)
        n = len(C) // 3
        assert abs(C[:n].mean() - C[-n:].mean()) < 1e-6

    def test_amplitude_stays_bounded(self):
        res = _run_real_trojan(periods=45)
        s = res["phi"] - memory.classify(res["phi"])[1]
        assert np.max(np.abs(s)) < 15.0        # seeded at 10 deg, no runaway


class TestObservedAmplitudeSpread:
    @pytest.mark.parametrize("amp_deg", [5.0, 15.0, 25.0])
    def test_real_trojan_amplitudes_still_store_the_bit(self, amp_deg):
        """Real Trojans librate with a spread of amplitudes (tens of degrees).
        Across that observed range the cell still classifies as a stored bit,
        not a horseshoe."""
        res = _run_real_trojan(libration_deg=amp_deg, periods=40)
        assert memory.classify(res["phi"])[0] == "L4"


class TestEarthTrojanHarderCase:
    def test_earth_trojan_period_is_longer_and_honestly_off_linear_theory(self):
        """2010 TK7 is a large-amplitude, partly-horseshoe Earth Trojan; linear
        theory (~222 yr) underestimates its real ~395-yr libration. We assert
        only the honest facts: the linear period at Earth's ratio is a few
        hundred years and far longer than Jupiter's, and a small-amplitude seed
        does store the bit. The real-object discrepancy is documented, not faked."""
        t_lib_earth = units.SUN_EARTH.years(memory.libration_period(memory.MU_SUN_EARTH))
        assert 150.0 < t_lib_earth < 300.0
        assert t_lib_earth > units.SUN_JUPITER.years(
            memory.libration_period(memory.MU_SUN_JUPITER))
        res = _run_real_trojan(libration_deg=8.0, periods=30, mu=memory.MU_SUN_EARTH)
        assert memory.classify(res["phi"])[0] == "L4"
