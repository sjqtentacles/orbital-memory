"""Contract for bit cooling.

A pinch-written bit librates wide (~±70°, inherited from the horseshoe's
phase-space area). Cooling shrinks it with tangential station-keeping burns
that damp the co-orbital pendulum's momentum (the radial offset da = r − 1),
each far below the erase threshold. Cooling must preserve the stored value,
move the Jacobi constant toward the deep-band floor C_L4, and produce a bit
that survives in the full inertial N-body engine.

(The docstring of orbital/cool.py records the three failed cooling schemes —
including the tempting-but-backwards 'raise C_J' certificate — as physics
documentation; these tests pin the one that works.)
"""

import numpy as np
import pytest

from orbital import cool, memory, nbody, rotating, theory, write


@pytest.fixture(scope="module")
def wide():
    """An honest wide bit: freshly pinch-written '1'."""
    return write.write_bit("1")


@pytest.fixture(scope="module")
def cooled(wide):
    return cool.cool(wide["yf"], t0=wide["t"][-1])


def _read_tail(run, window_orbits=20):
    t = run["t"]
    idx = int(np.searchsorted(t, t[-1] - window_orbits * memory.PERIOD))
    return memory.classify(run["phi"][idx:])


class TestCooling:
    def test_wide_bit_starts_wide(self, wide):
        bit, center, amp = write.read(wide)
        assert bit == "1"
        assert amp > 55.0

    def test_cooling_shrinks_amplitude(self, cooled):
        label, center, amp = _read_tail(cooled)
        assert label == "L4"                 # same value, deeper well
        assert amp < cool.TARGET_AMP

    def test_cooling_needs_only_a_few_burns(self, cooled):
        assert 1 <= cooled["n_kicks"] <= 8

    def test_jacobi_gap_shrinks(self, cooled):
        """Cooling moves C_J toward the deep-band floor C_L4 — but no
        stronger claim is honest: C_J mixes the slow pendulum energy with
        fast epicyclic energy, so a deep-but-slightly-eccentric tadpole can
        sit BELOW C_L4 (ours does) while a kicked-erased orbit sits nearby.
        C_J stratifies the dynamics; it does not classify orbits. The
        operational tests (amplitude, readback, nbody hold) carry the
        correctness burden; this one just pins the direction of travel."""
        gap0, gap1 = cooled["jacobi_gap"]
        assert gap1 < 0.6 * gap0

    def test_burns_are_gentle(self, cooled):
        assert np.max(cooled["kick_sizes"]) <= cool.CAP
        assert cool.CAP < memory.ERASE_KICK / 4

    def test_cooled_bit_survives_full_nbody_hold(self, cooled):
        # the hold's clock must continue from the handoff time: jacobi_of
        # reconstructs the rotating frame from t, and the primaries' phase
        # is t mod 2*pi (the write tests dodge this because their t_end is
        # an exact orbit multiple; cooling ends at an arbitrary phase)
        t0 = float(cooled["t"][-1])
        bodies = rotating.to_inertial_bodies(cooled["yf"], t0, mu=memory.MU)
        held = nbody.integrate(bodies, t0 + 40 * memory.PERIOD,
                               n_samples=2000, t0=t0)
        assert memory.classify(memory.resonant_angle(held))[0] == "L4"
        C = theory.jacobi_of(held)
        assert (C.max() - C.min()) < 1e-9

    def test_cooling_is_deterministic(self, wide, cooled):
        again = cool.cool(wide["yf"], t0=wide["t"][-1])
        assert np.array_equal(again["xy"], cooled["xy"])
