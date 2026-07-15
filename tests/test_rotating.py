"""Contract for the rotating-frame restricted integrator (written first).

The rotating frame is the memory's natural coordinate system: the primaries
sit still, L4/L5 are fixed points, and the resonant angle is just atan2(y, x).
It must agree with the full inertial N-body integrator wherever they overlap.
This is the analysis/phase-portrait engine; the memory operations run in the
inertial nbody engine.
"""

import numpy as np
import pytest

from orbital import memory, nbody, rotating, theory


class TestAgreesWithFullNBody:
    def test_trajectory_matches_inertial_integrator(self):
        """Same L4 tadpole, two engines: the rotating-frame restricted run must
        reproduce the full 3-body inertial run (transformed) to high accuracy."""
        cell = memory.make_cell("L4", libration_deg=6.0)
        full = nbody.integrate(cell, 10 * memory.PERIOD, n_samples=500)
        ref = rotating.to_rotating_frame(full, 2)

        state0 = rotating.from_inertial(cell[2]["pos"], cell[2]["vel"], t=0.0)
        run = rotating.integrate(state0, 10 * memory.PERIOD, n_samples=500)
        assert np.max(np.linalg.norm(run["xy"] - ref, axis=1)) < 1e-6

    def test_to_rotating_frame_helper(self):
        """The shared helper equals the textbook transform, and a body on the
        corotation circle maps to a fixed point."""
        cell = memory.make_cell("L4", libration_deg=0.0)
        full = nbody.integrate(cell, 5 * memory.PERIOD, n_samples=200)
        rp = rotating.to_rotating_frame(full, 2)
        assert np.max(np.linalg.norm(rp - rp[0], axis=1)) < 1e-3
        # planet itself is pinned at (1-mu, 0) in this frame
        planet = rotating.to_rotating_frame(full, 1)
        assert np.allclose(planet, [1 - memory.MU, 0.0], atol=1e-9)

    def test_jacobi_conserved_at_constant_mu(self):
        state0 = rotating.lagrange_tadpole("L4", libration_deg=6.0)
        run = rotating.integrate(state0, 40 * memory.PERIOD, n_samples=2000)
        C = run["jacobi"]
        assert (C.max() - C.min()) < 1e-9

    def test_L4_is_a_fixed_point(self):
        pos = np.array([0.5 - memory.MU, np.sqrt(3) / 2])
        run = rotating.integrate(np.array([*pos, 0.0, 0.0]),
                                 20 * memory.PERIOD, n_samples=400)
        assert np.max(np.linalg.norm(run["xy"] - pos, axis=1)) < 1e-3


class TestPhiReadout:
    def test_phi_is_plain_polar_angle(self):
        """In the rotating frame the resonant angle is just atan2(y,x)."""
        state0 = rotating.lagrange_tadpole("L4", libration_deg=6.0)
        run = rotating.integrate(state0, 30 * memory.PERIOD, n_samples=1500)
        label, center, amp = memory.classify(run["phi"])
        assert label == "L4"
        assert abs(center - 60) < 6

    def test_blank_medium_circulates(self):
        """A co-orbital particle outside the co-orbital zone (da larger than
        the horseshoe width ~ mu^(1/3)) drifts through all longitudes — an
        unwritten (blank) medium reads 'erased' and truly circulates."""
        state0 = rotating.circular_coorbital(phi0_deg=90.0, da=0.06)
        run = rotating.integrate(state0, 60 * memory.PERIOD, n_samples=3000,
                                 mu=1e-5)
        assert memory.classify(run["phi"])[0] == "erased"
        unwrapped = np.unwrap(np.radians(run["phi"]))
        drift = unwrapped[-1] - unwrapped[0]
        assert abs(drift) > 2 * np.pi                       # full circulation
        assert drift == pytest.approx(-1.5 * 0.06 * 60 * memory.PERIOD,
                                      rel=0.25)             # matches -1.5 da t


class TestCallableMu:
    def test_constant_callable_matches_scalar(self):
        """mu may still be passed as a callable (kept for generality); a
        constant callable must reproduce the scalar path exactly."""
        state0 = rotating.lagrange_tadpole("L4", libration_deg=6.0)
        a = rotating.integrate(state0, 10 * memory.PERIOD, n_samples=400)
        b = rotating.integrate(state0, 10 * memory.PERIOD, n_samples=400,
                               mu=lambda t: memory.MU)
        assert np.allclose(a["xy"], b["xy"], atol=1e-9)
