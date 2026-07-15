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
    def test_dv_is_finite_and_reported_in_mps(self):
        _, dv = write.insert("1")
        assert dv > 0.0
        mps = units.SUN_JUPITER.mps(dv)
        assert 50.0 < mps < 1000.0            # a plausible L4 station-insertion cost

    def test_bigger_transfer_costs_more(self):
        _, small = write.insert("1", approach_da=0.03)
        _, big = write.insert("1", approach_da=0.10)
        assert big > small

    def test_both_bits_cost_the_same(self):
        _, a = write.insert("1")
        _, b = write.insert("0")
        assert a == pytest.approx(b, rel=1e-9)


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
