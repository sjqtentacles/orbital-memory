"""Contract for the bullet gate (written before orbital/gate.py exists).

The unification step: a flyby — the mechanism slingshot-computing uses for
logic — acts on a stored orbital bit. A massive bullet on a fast hyperbolic
pass is aimed to shave past the stored moonlet: bullet present, the bit
erases; bullet absent (or merely grazing at a safe distance), the bit
survives. Logic input, memory output.

Physics the contract encodes: the bullet is MASSIVE (a massless one exerts
no force), so it kicks the system barycenter — readout must use the
COM-corrected resonant angle; the bullet must leave (hyperbolic); the gate
must erase the BIT, not the cell (primaries' orbit intact); and the
close-encounter runs get a looser energy bound, asserted explicitly.
"""

import numpy as np
import pytest

from orbital import gate, memory, nbody


@pytest.fixture(scope="module")
def cell():
    """The deep 6-degree L4 bit: hardest to erase, cheapest to make."""
    return memory.make_cell("L4", libration_deg=6.0)


@pytest.fixture(scope="module")
def aimed(cell):
    psi, achieved = gate.aim(cell, miss=gate.MISS_ERASE)
    return {"psi": psi, "achieved": achieved}


class TestAiming:
    def test_aim_converges_to_requested_miss(self, aimed):
        assert abs(aimed["achieved"] - gate.MISS_ERASE) < 0.5 * gate.MISS_ERASE

    def test_bullet_is_hyperbolic_and_exits(self, cell, aimed):
        res = gate.fire(cell, aimed["psi"])
        r_b = np.linalg.norm(res["traj"][3], axis=1)
        assert r_b[-1] > gate.BULLET_R0          # left the system
        # outbound: radial distance increasing at the end
        assert r_b[-1] > r_b[-100]


class TestConditionalErase:
    def test_bullet_present_erases(self, cell, aimed):
        res = gate.fire(cell, aimed["psi"])
        label, _, _ = memory.classify(gate.resonant_angle_com(res))
        assert label == "erased"

    def test_bullet_absent_bit_survives(self, cell, aimed):
        """Same duration, same COM readout path — only the bullet differs."""
        res = nbody.integrate([dict(b) for b in cell], gate.T_VERDICT,
                              n_samples=2000)
        label, center, _ = memory.classify(gate.resonant_angle_com(res))
        assert label == "L4"
        assert abs(center - 60) < 10

    def test_graze_does_not_erase(self, cell):
        """A pass at MISS_SAFE must leave the bit intact — the locality datum
        dual-rail writing rests on."""
        psi, achieved = gate.aim(cell, miss=gate.MISS_SAFE)
        assert achieved > gate.MISS_ERASE * 3
        res = gate.fire(cell, psi)
        label, _, _ = memory.classify(gate.resonant_angle_com(res))
        assert label == "L4"

    def test_gate_erases_the_bit_not_the_cell(self, cell, aimed):
        """The primaries' orbit must ride out the shot. The bullet
        (m = 2e-4) legitimately wobbles the separation by ~2 G m_b / v ~
        a few 1e-3; the bound scales with BULLET_M and just needs to show
        the cell hardware stays a cell."""
        res = gate.fire(cell, aimed["psi"])
        sep = np.linalg.norm(res["traj"][0] - res["traj"][1], axis=1)
        assert np.max(np.abs(sep - 1.0)) < 25 * gate.BULLET_M

    def test_close_encounter_numerics(self, cell, aimed):
        res = gate.fire(cell, aimed["psi"])
        assert res["energy_drift"] < 1e-6       # looser: near-collision pass

    def test_gate_is_deterministic(self, cell, aimed):
        a = gate.fire(cell, aimed["psi"])
        b = gate.fire(cell, aimed["psi"])
        assert np.array_equal(a["traj"], b["traj"])


class TestComReadout:
    def test_com_correction_matters(self, cell, aimed):
        """The bullet drags the barycenter, so the origin-anchored readout
        and the COM-corrected one must visibly differ during/after the pass
        — the reason resonant_angle_com exists."""
        res = gate.fire(cell, aimed["psi"])
        raw = memory.resonant_angle(res)
        com = gate.resonant_angle_com(res)
        assert np.max(np.abs(raw - com)) > 0.5   # degrees
