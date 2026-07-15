"""Contract for the eccentric (elliptic restricted) cell.

Real planets are on eccentric orbits — Jupiter's is e ~ 0.0489. The bit must
survive when the planet breathes in and out once per year. There is no exact
Jacobi integral in the elliptic problem (the frame pulsates), so the honest
invariants are: the bit is retained, its libration stays bounded, and the
motion breathes at the planet's orbital frequency.
"""

import numpy as np
import pytest

from orbital import memory, nbody


def _run(e, state="L4", libration_deg=6.0, periods=60, n_samples=4000):
    cell = memory.make_cell_eccentric(state, e=e, libration_deg=libration_deg)
    res = nbody.integrate(cell, periods * memory.PERIOD, n_samples=n_samples)
    res["phi"] = memory.resonant_angle(res)
    res["sep"] = np.linalg.norm(res["traj"][0] - res["traj"][1], axis=1)
    return res


class TestReducesToCircular:
    def test_zero_eccentricity_has_constant_separation(self):
        res = _run(0.0, periods=20)
        assert np.ptp(res["sep"]) < 1e-6

    def test_zero_eccentricity_matches_circular_cell(self):
        ecc = memory.classify(_run(0.0, periods=40)["phi"])
        circ = memory.classify(memory.hold("L4", periods=40, libration_deg=6.0)["phi"])
        assert ecc[0] == circ[0] == "L4"
        assert ecc[2] == pytest.approx(circ[2], abs=1.0)


class TestBitSurvivesEccentricity:
    @pytest.mark.parametrize("state", ["L4", "L5"])
    def test_real_jupiter_eccentricity_holds_the_bit(self, state):
        assert memory.classify(_run(0.0489, state=state)["phi"])[0] == state

    def test_separation_breathes_at_one_plus_minus_e(self):
        res = _run(0.0489, periods=20)
        assert res["sep"].min() == pytest.approx(1 - 0.0489, abs=2e-3)
        assert res["sep"].max() == pytest.approx(1 + 0.0489, abs=2e-3)

    def test_libration_stays_bounded(self):
        _, _, amp = memory.classify(_run(0.0489, periods=80)["phi"])
        assert amp < 25.0                         # seeded 6 deg, no runaway


class TestBreathingIsAtOrbitalFrequency:
    def test_separation_period_is_one_orbit(self):
        """The planet returns to periapsis once per orbital period (2*pi)."""
        res = _run(0.0489, periods=20, n_samples=6000)
        sep = res["sep"] - res["sep"].mean()
        t = res["t"]
        crossings = t[np.where(np.diff(np.sign(sep)) > 0)[0]]  # rising, once/orbit
        period = np.mean(np.diff(crossings))
        assert period == pytest.approx(memory.PERIOD, rel=0.03)
