"""Contract for retention — how long a bit lasts, honestly.

The testable claims live in the escape-observable band near the separatrix
(deep bits are censored survivors of any feasible run, so 'deep outlives wide'
is only meaningful as: deep survives, wide escapes). Relative/platform-invariant
per the chaotic-dispersion discipline.
"""

import pytest

from orbital import memory, retention


class TestEscapeDetection:
    def test_deep_bit_is_a_censored_survivor(self):
        """A deep bit does not escape within a feasible horizon — its true
        lifetime is beyond what we integrate (Nekhoroshev-stable)."""
        assert retention.is_censored_survivor(20.0, max_orbits=300)

    @pytest.mark.slow
    def test_near_separatrix_bit_escapes(self):
        """A bit seeded close to the L3 separatrix escapes in a feasible run."""
        e = retention.escape_time(memory.make_cell("L4", libration_deg=82.0),
                                  max_orbits=300)
        assert e is not retention.CENSORED and e <= 300


class TestRetentionDatasheet:
    @pytest.mark.slow
    def test_escape_time_drops_toward_the_separatrix(self):
        """Within the escape band, a wider bit (closer to the separatrix) escapes
        sooner than a less-wide one — the steep retention falloff."""
        scan = retention.escape_scan([74.0, 82.0], max_orbits=400)
        by_amp = {r["amp"]: r for r in scan}
        assert not by_amp[74.0]["censored"] and not by_amp[82.0]["censored"]
        assert by_amp[82.0]["escape_orbits"] <= by_amp[74.0]["escape_orbits"]

    @pytest.mark.slow
    def test_deep_survives_while_wide_escapes(self):
        """The honest 'deep outlives wide': a deep bit is censored (survives the
        horizon) while a near-separatrix bit escapes within it."""
        scan = retention.escape_scan([25.0, 82.0], max_orbits=300)
        by_amp = {r["amp"]: r for r in scan}
        assert by_amp[25.0]["censored"]
        assert not by_amp[82.0]["censored"]
