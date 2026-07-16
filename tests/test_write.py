"""Contract for WRITE by orbit insertion.

Writing is placing the body in the chosen island: insert '1' -> it librates at
L4, insert '0' -> L5, and the cost is a real insertion burn reported in m/s.
No planet is grown. This replaces the old mass-growth write entirely.
"""

import numpy as np
import pytest

from orbital import memory, nbody, units, write


class TestInsertionWritesTheBit:
    def test_write_one_lands_at_L4(self):
        assert write.read(write.write_bit("1"))[0] == "1"

    def test_write_zero_lands_at_L5(self):
        assert write.read(write.write_bit("0"))[0] == "0"

    def test_written_amplitude_matches_request(self):
        run = write.write_bit("1", amplitude_deg=6.0)
        _, _, amp = write.read(run)
        assert amp == pytest.approx(6.0, abs=2.5)

    def test_wider_request_gives_wider_bit(self):
        narrow = write.read(write.write_bit("1", amplitude_deg=6.0))[2]
        wide = write.read(write.write_bit("1", amplitude_deg=30.0))[2]
        assert wide > narrow + 15.0


class TestInsertionBurn:
    def test_dv_is_the_actually_applied_burn(self):
        """dv is the magnitude of the velocity change the write applies, not an
        accounting number: written velocity == arrival velocity + the burn."""
        cell, dv = write.insert("1")
        arrival = write.arrival_state("1")
        applied = np.linalg.norm(np.array(cell[2]["vel"])
                                 - np.array(arrival[2]["vel"]))
        assert applied == pytest.approx(dv, rel=1e-12)
        assert 100.0 < units.SUN_JUPITER.mps(dv) < 2000.0   # realistic L4 insertion

    def test_the_burn_is_what_writes_the_bit(self):
        """Without the insertion burn the body coasts in on the transfer and is
        NOT a stored bit (erased); applying the burn captures it into the island.
        The maneuver is load-bearing, not decorative."""
        arrival = write.arrival_state("1")
        run = nbody.integrate(arrival, 45 * memory.PERIOD, n_samples=2700)
        run["phi"] = memory.resonant_angle(run)
        assert write.read(run)[0] == "erased"          # no burn -> no bit
        assert write.read(write.write_bit("1"))[0] == "1"   # burn -> bit

    def test_bigger_transfer_costs_more(self):
        _, small = write.insert("1", approach_da=0.08)
        _, big = write.insert("1", approach_da=0.18)
        assert big > small


class TestBlankMedium:
    def test_blank_reads_erased(self):
        run = nbody.integrate(write.blank_cell(), 60 * memory.PERIOD, n_samples=3000)
        run["phi"] = memory.resonant_angle(run)
        assert write.read(run)[0] == "erased"


class TestWrittenBitIsHonest:
    def test_survives_full_nbody_hold_with_conserved_jacobi(self):
        from orbital import theory
        run = write.write_bit("1", periods=45)
        C = theory.jacobi_of(run)
        n = len(C) // 3
        assert abs(C[:n].mean() - C[-n:].mean()) < 1e-6     # no secular drift
        assert write.read(run)[0] == "1"

    def test_write_is_deterministic(self):
        a = write.write_bit("1")
        b = write.write_bit("1")
        assert np.array_equal(a["traj"], b["traj"])
