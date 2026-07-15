"""Contract for the real Saturn perturber.

Saturn is why the Jupiter Trojan swarm is sculpted rather than pristine — but
the full erosion (secular resonances shaking Trojans loose) is a Gyr process,
far beyond any feasible integration here. What we can show, and assert, is the
ONSET: over a short run Saturn measurably pumps a Trojan's libration compared
to a Saturn-free control. The assertion is the relative invariant (Saturn
pumps more), never a platform-specific amplitude — chaotic trajectories
disperse across BLAS builds.
"""

import numpy as np
import pytest

from orbital import memory, nbody


def _peak_amplitude(with_saturn, periods=130):
    cell = memory.sun_jupiter_saturn_cell("L4", libration_deg=10.0,
                                          with_saturn=with_saturn)
    res = nbody.integrate(cell, periods * memory.PERIOD, n_samples=periods * 40)
    return memory.classify(memory.resonant_angle(res))


@pytest.mark.slow
class TestSaturnErodesTheTrojan:
    def test_saturn_pumps_the_libration(self):
        label_c, _, amp_control = _peak_amplitude(False)
        label_s, _, amp_saturn = _peak_amplitude(True)
        assert amp_saturn > 1.5 * amp_control      # the onset of erosion
        # over this short horizon the bit is pumped, not yet lost (escape is Gyr)
        assert label_c == "L4"

    def test_control_barely_moves(self):
        label, _, amp = _peak_amplitude(False)
        assert label == "L4"
        assert amp < 15.0                          # seeded 10 deg, nearly steady
