"""The WRITE operation: same blank, same pulse — timing selects the bit.

The full memory cycle is now closed: blank -> write -> hold -> read (-> erase).
The last test hands a freshly written bit to the full inertial N-body
integrator and holds it there — write in the fast engine, verify in the
honest one.
"""

import numpy as np
import pytest

from orbital import memory, nbody, rotating, write


class TestBlankMedium:
    def test_blank_is_coorbital_but_bitless(self):
        run = rotating.integrate(write.blank(), 60 * memory.PERIOD,
                                 n_samples=1500, mu=write.MU0)
        label, _, _ = memory.classify(run["phi"])
        assert label == "erased"          # horseshoe: no stored bit
        # but genuinely co-orbital: no net circulation
        u = np.unwrap(np.radians(run["phi"]))
        assert abs(u[-1] - u[0]) < 1.5 * np.pi


class TestWrite:
    def test_write_one(self):
        run = write.write_bit("1")
        bit, center, amp = write.read(run)
        assert bit == "1"
        assert center > 0                 # leading side (L4)
        assert amp < 120.0                # a genuine tadpole, not a horseshoe

    def test_write_zero(self):
        run = write.write_bit("0")
        bit, center, amp = write.read(run)
        assert bit == "0"
        assert center < 0                 # trailing side (L5)
        assert amp < 120.0

    def test_same_hardware_only_timing_differs(self):
        """The write pulse is identical in shape; only its start time changes.
        This is the defining property of the mechanism."""
        assert write.WRITE_DELAY["1"] != write.WRITE_DELAY["0"]
        b1, b2 = write.blank(), write.blank()
        assert np.array_equal(b1, b2)

    def test_write_is_deterministic(self):
        a = write.write_bit("1")
        b = write.write_bit("1")
        assert np.array_equal(a["xy"], b["xy"])

    def test_read_window_is_in_orbits_not_samples(self):
        """read() must give the same verdict regardless of sampling density —
        the window is physical time (orbits), not an index count."""
        dense = write.write_bit("1", n_samples=3000)
        sparse = write.write_bit("1", n_samples=600)
        assert write.read(dense)[0] == "1"
        assert write.read(sparse)[0] == "1"

    def test_scan_delays_uses_production_path(self):
        """The datasheet scan must be the production write with a different
        delay — same blank, same pulse shape, same read. Canonical delays in
        the scan must reproduce the canonical bits."""
        delays = [write.WRITE_DELAY["1"], write.WRITE_DELAY["0"]]
        out = write.scan_delays(delays)
        assert out[0]["bit"] == "1" and out[1]["bit"] == "0"


class TestWriteThenHold:
    @pytest.mark.parametrize("bit,label", [("1", "L4"), ("0", "L5")])
    def test_written_bit_survives_in_full_nbody(self, bit, label):
        """Write with the fast rotating-frame engine, then hold the result in
        the full inertial 3-body integrator for 40 more orbits: the bit must
        persist across engines — no artifact of the restricted approximation."""
        run = write.write_bit(bit)
        bodies = write.written_cell_bodies(run)
        held = nbody.integrate(bodies, 40 * memory.PERIOD, n_samples=2000)
        phi = memory.resonant_angle(held)
        got, _, _ = memory.classify(phi)
        assert got == label
        # the moonlet's own invariant (energy_drift only sees the primaries)
        from orbital import theory
        C = theory.jacobi_of(held)
        assert (C.max() - C.min()) < 1e-9
