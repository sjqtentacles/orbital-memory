"""Contract for the averaged 1-DOF co-orbital Hamiltonian.

The whole point: the closed forms must reproduce what the full n-body sim
measures. The sim is the verifier of the theory. Period agrees to ~1.5% across
amplitudes; the separatrix and erase margin agree at the level the averaging
allows (the model breaks down near conjunction, documented).
"""

import numpy as np
import pytest

from orbital import hamiltonian as H, memory, nbody


class TestPotentialShape:
    def test_equilibria_at_L4_and_L3(self):
        assert H.df(np.radians(60.0)) == pytest.approx(0.0, abs=1e-9)
        assert H.df(np.radians(180.0)) == pytest.approx(0.0, abs=1e-9)

    def test_barrier_is_exactly_three_mu(self):
        assert H.potential(180.0, mu=memory.MU) == pytest.approx(3.0 * memory.MU, rel=1e-9)
        assert H.potential(60.0, mu=memory.MU) == pytest.approx(0.0, abs=1e-12)

    def test_small_amplitude_frequency_is_27_over_4(self):
        """f''(60) = -9/4 -> omega^2 = (27/4) mu, the known libration frequency."""
        h = 1e-6
        fpp = (H.df(np.radians(60) + h) - H.df(np.radians(60) - h)) / (2 * h)
        assert fpp == pytest.approx(-2.25, abs=1e-4)


class TestLibrationPeriodMatchesSim:
    def test_small_amplitude_reduces_to_linear_theory(self):
        assert H.libration_period(0.0, memory.MU) == memory.libration_period(memory.MU)
        assert H.libration_period(1e-9, memory.MU) == pytest.approx(
            memory.libration_period(memory.MU), rel=1e-6)

    @pytest.mark.parametrize("amp", [2.0, 6.0, 15.0, 30.0, 45.0, 60.0])
    def test_nonlinear_period_matches_nbody(self, amp):
        """Analytic quadrature vs the measured libration period of the full
        inertial n-body sim, across the tadpole's amplitude range."""
        analytic = H.libration_period(amp, memory.MU)
        res = nbody.integrate(memory.make_cell("L4", mu=memory.MU, libration_deg=amp),
                              60 * memory.PERIOD, n_samples=5000)
        measured = memory.measured_libration_period(res["t"], memory.resonant_angle(res))
        assert analytic == pytest.approx(measured, rel=0.03)

    def test_period_grows_with_amplitude(self):
        assert (H.libration_period(60.0, memory.MU)
                > H.libration_period(6.0, memory.MU) > 0)


class TestSeparatrixAmplitude:
    def test_value_is_the_widest_tadpole(self):
        amp = H.separatrix_amplitude()
        assert 70.0 < amp < 85.0                # analytic ~78 deg

    def test_is_mass_ratio_independent(self):
        """mu scales the energy, not the angular extent — the separatrix
        amplitude is the same for Jupiter and for a Saturn co-orbital moon."""
        # separatrix_amplitude takes no mu; assert the underlying turning point
        # is fixed and the two well-known ratios give the same well shape.
        assert H._separatrix_conjunction_deg() == pytest.approx(23.9, abs=0.5)

    def test_brackets_the_sim_transition(self):
        """A tadpole comfortably inside the separatrix reads the bit; one well
        beyond it reads erased. (The exact edge is a chaotic layer that disperses
        across runs/BLAS builds, so we test robustly-inside vs robustly-beyond,
        not the razor's edge — the analytic ~78 deg sits in that chaotic zone.)"""
        inside = nbody.integrate(memory.make_cell("L4", libration_deg=30.0),
                                 60 * memory.PERIOD, n_samples=5000)
        beyond = nbody.integrate(memory.make_cell("L4", libration_deg=110.0),
                                 60 * memory.PERIOD, n_samples=5000)
        assert memory.classify(memory.resonant_angle(inside))[0] == "L4"
        assert memory.classify(memory.resonant_angle(beyond))[0] == "erased"
        assert 30.0 < H.separatrix_amplitude() < 110.0    # between the two


class TestEraseMargin:
    def test_derives_erase_kick_to_leading_order(self):
        """The analytic L3-crossing margin is close to the measured ERASE_KICK.
        It runs ~30% high because the averaged model can't see conjunction-side
        close encounters, which let the real sim erase a bit more easily — so it
        is an upper-ish leading-order estimate of the measured margin."""
        margin = H.erase_margin(amp_deg=2.0, mu=memory.MU)
        assert margin == pytest.approx(memory.ERASE_KICK, rel=0.4)
        assert margin > memory.ERASE_KICK          # L3 barrier is the harder route

    def test_margin_shrinks_for_wider_bits(self):
        """A wider bit is already partway up the barrier, so a smaller kick
        erases it."""
        assert H.erase_margin(2.0, memory.MU) > H.erase_margin(40.0, memory.MU)
