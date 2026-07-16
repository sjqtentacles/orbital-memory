"""Contract for the multi-bit register.

N co-orbital cells at spaced radii round-trip a bit-string; capacity is bounded
by mutual stability (too-close hosts go chaotic). These are the honest "more
than one bit" tests.
"""

import numpy as np
import pytest

from orbital import register


class TestRoundTrip:
    def test_nibble_round_trips(self):
        run = register.hold_register("1101")
        assert register.read_register(run) == "1101"

    @pytest.mark.slow
    @pytest.mark.parametrize("bits", ["0110", "1010", "0000", "1111"])
    def test_patterns_round_trip(self, bits):
        assert register.read_register(register.hold_register(bits)) == bits

    def test_deterministic(self):
        a = register.hold_register("1010")
        b = register.hold_register("1010")
        assert np.array_equal(a["traj"], b["traj"])


class TestCapacityIsStabilityBounded:
    @pytest.mark.slow
    def test_safe_spacing_has_low_crosstalk(self):
        """At the default spacing the neighbours induce only a few degrees of
        libration — the cells stay clean bits."""
        assert register.crosstalk("101", spacing=register.SPACING) < 20.0

    @pytest.mark.slow
    def test_too_close_goes_chaotic(self):
        """Pack the hosts too tightly and the register destabilizes — the
        capacity bound is real, not cosmetic."""
        assert register.crosstalk("101", spacing=1.5) > 30.0


class TestGeometry:
    def test_radii_are_geometric(self):
        rs = register.radii(4, spacing=1.8)
        assert rs[0] == 1.0
        assert rs[1] == pytest.approx(1.8)
        assert rs[3] == pytest.approx(1.8 ** 3)

    def test_register_has_star_plus_pair_per_bit(self):
        bodies = register.make_register("101")
        assert len(bodies) == 1 + 2 * 3            # star + (planet,trojan)*3
        assert bodies[0]["m"] == 1.0               # star
        assert bodies[2]["m"] == 0.0               # first trojan is massless

    def test_total_momentum_is_zeroed(self):
        bodies = register.make_register("1010")
        P = sum(b["m"] * np.array(b["vel"]) for b in bodies)
        assert np.allclose(P, 0.0, atol=1e-12)
