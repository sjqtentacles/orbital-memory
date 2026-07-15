"""Erase, and the rewrite wall — a pinned negative result.

The erase pulse (mu shrunk back to blank) releases the stored bit into a
bounded, bit-less horseshoe: that is tested as working. What is ALSO tested
— deliberately, as physics — is that the released medium cannot be
re-written at this design point: release inflates the offset past the
recapture ceiling, and the slow-erase 'fix' ejects the moonlet entirely.
Write-once, erase-once; the rewrite channel is the roadmap's open problem.
"""

import numpy as np
import pytest

from orbital import memory, rewrite, write


@pytest.fixture(scope="module")
def erased():
    """The canonical post-erase state: blank -> write '1' -> store -> erase
    -> settle, one continuous run."""
    return rewrite.erased_prefix("1")


class TestErase:
    def test_erase_restores_a_blank(self, erased):
        label, center, amp = memory.classify(erased["phi"][-500:])
        assert label == "erased"             # no stored bit anymore

    def test_released_medium_is_bounded_not_ejected(self, erased):
        """The release lands in a horseshoe: no net circulation, offset in
        the co-orbital band (not flung away)."""
        tail = erased["phi"][-500:]
        u = np.unwrap(np.radians(tail))
        assert abs(u[-1] - u[0]) < 1.5 * np.pi
        da = rewrite.horseshoe_offset(erased)
        assert da < 0.15

    def test_release_inflates_the_medium(self, erased):
        """The measured invariant jump: captured at da = write.DA = 0.030,
        released at ~0.07 — past the recapture ceiling. The number behind
        the write-once conclusion."""
        da = rewrite.horseshoe_offset(erased)
        lo, hi = rewrite.RELEASE_DA
        assert lo < da < hi
        assert da > rewrite.RECAPTURE_CEILING

    def test_erase_is_deterministic(self, erased):
        again = rewrite.erased_prefix("1")
        assert np.array_equal(again["xy"], erased["xy"])


class TestRewriteWall:
    def test_recapture_ceiling_is_analytic(self):
        """The tadpole band can never outrun the Hill radius in this mu
        range: W/r_H = 2.35 mu^(1/6) < 1, so any medium released wider than
        ~RECAPTURE_CEILING has its would-be pinch inside Hill-scattering
        territory."""
        for mu in (3e-4, 1e-3, 3e-3):
            W = np.sqrt(8 * mu / 3)          # tadpole band full width
            r_H = (mu / 3) ** (1 / 3)        # Hill radius
            assert W / r_H < 1.0

    @pytest.mark.slow
    def test_no_second_write_delay_recaptures(self):
        """The recapture experiment that measured the wall: standard growth
        pulses at several delays from the settled post-erase state — every
        one reads 'erased'."""
        T = memory.PERIOD
        results = rewrite.scan_rewrite_delays([0.0, 8 * T, 16 * T])
        assert all(r["bit"] == "erased" for r in results)

    @pytest.mark.slow
    def test_slower_erase_ejects_instead_of_helping(self):
        """The tempting fix fails harder: at 3x the erase duration the
        moonlet lingers in the separatrix layer and is ejected from the
        co-orbital region entirely."""
        run = rewrite.erased_prefix("1",
                                    t_erase_duration=300 * memory.PERIOD)
        assert rewrite.horseshoe_offset(run) > 1.0  # gone


class TestCooledBitsAreLocked:
    @pytest.mark.slow
    def test_erase_pulse_does_not_release_a_cooled_bit(self):
        """Cooling write-protects: a cooled bit's phase-space area is below
        the small-mu island's, so the same erase pulse leaves it stored."""
        from orbital import cool
        wide = write.write_bit("1")
        cooled = cool.cool(wide["yf"], t0=wide["t"][-1])
        run = rewrite.erase(cooled["yf"], float(cooled["t"][-1]))
        label, _, _ = memory.classify(run["phi"][-500:])
        assert label == "L4"                 # still a bit, not released
