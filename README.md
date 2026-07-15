# Orbital Memory

**A complete, nonvolatile memory cell made of gravity — write, hold, read, erase — that stores its bit the way Jupiter's Trojans and Saturn's co-orbital moons store theirs.**

Where its sibling project [slingshot-computing](https://github.com/sjqtentacles/slingshot-computing) does *logic* with transient gravitational flybys — and is fundamentally memoryless — this one is the other half of a computer: **storage**. A bit is stored as *which stable island a moonlet orbits in*. Nothing acts but `F = Gm₁m₂/r²` — the one non-gravitational term in the engine is an optional drag switch that exists *solely to prove that dissipation destroys this memory*.

<p align="center">
  <img src="docs/write.gif" width="520" alt="Rotating frame: a horseshoe moonlet sweeps the ring while the secondary grows, then pinches into a tadpole at L4 — the bit is written">
</p>

<p align="center"><em><b>Writing a bit by growing a planet.</b> The blank medium is a horseshoe moonlet; a slow mass-growth pulse pinches it into L4 (bit 1) or L5 (bit 0) — selected purely by pulse <b>timing</b>. This is the mechanism by which a growing Jupiter captured its Trojans.</em></p>

---

## The full memory cycle

| Operation | Mechanism | Physics |
|---|---|---|
| **HOLD** | tadpole libration around L4/L5 | topological protection (invariant island, KAM) |
| **WRITE** | grow the secondary's mass — the horseshoe pinch | adiabatic capture (Henrard, Neishtadt) |
| **READ** | which side the resonant angle librates on | the separatrix-crossing classifier |
| **ERASE** | any kick past the noise margin | separatrix crossing, `C_J` drops below the held value |

- librating around **L4** (60° ahead of the secondary) → reads **1**
- librating around **L5** (60° behind) → reads **0**
- horseshoe / circulation → **erased / blank**

**This is moon-scale hardware.** Saturn's moons **Telesto** and **Calypso** ride Tethys's L4/L5, **Helene** and **Polydeuces** ride Dione's, and **Janus & Epimetheus** live on the horseshoe orbits our blank medium uses. The dynamics depend only on the mass ratio `μ`, so the same code covers star+planet, planet+moon, and moon+moonlet.

## The states of the medium

<p align="center">
  <img src="docs/anatomy.png" width="560" alt="Rotating-frame orbit families: cyan tadpole at L4, pink tadpole at L5, amber horseshoe sweeping the ring, dashed grey circulation">
</p>

## The write: timing is the data

Same blank horseshoe, same growth pulse — only the **firing time** differs. The moonlet alternates sides of the ring as it runs its horseshoe; the pinch captures it into whichever island it is transiting when the growing tadpole band engulfs its orbit:

<p align="center">
  <img src="docs/timing.png" width="720" alt="Write timing diagram: pulse delay of 12 orbits writes 1, delay 30 (the canonical `WRITE_DELAY['0']`) writes 0, other delays remain erased guard bands">
</p>

Why it must be this way — three write mechanisms, two of which provably fail:

1. **Impulsive kicks cannot write.** A conservative kick moves the state along the energy surface; measured, every super-threshold kick scatters the moonlet out of the resonance. Kicks only erase.
2. **Dissipation cannot write.** L4/L5 are Coriolis-stabilized potential extrema — drag *destabilizes* them (verified in code). There is no attractor to relax into.
3. **A slow parameter change can.** Adiabatic capture needs no aim, only timing — robustness is the adiabatic theorem's gift. The write pulse is `μ(t)`: mass transferred from primary to secondary at fixed total, so the circular kinematics stay exact throughout.

The written bit then **survives an engine swap**: written in the fast rotating-frame integrator, each bit is handed to the full inertial N-body integrator and held for 40 more orbits (test-enforced).

## The energy landscape

<p align="center">
  <img src="docs/landscape.png" width="560" alt="The rotating-frame effective potential with zero-velocity curves, the five Lagrange points, and the two stored tadpole orbits">
</p>

The rigorous backbone is the **Jacobi constant** `C_J` — the one conserved quantity of the circular restricted three-body problem. A held bit sits at the exact triangular value `C_L4 = 3 − μ(1−μ)` (sim matches to 2e-4 and conserves `C_J` to 1e-10, both test-enforced; measured drift is ~1e-12); writing energy in (a kick) lowers `C_J` toward the separatrix; past it, the bit erases. The noise margin — kicks below `memory.ERASE_KICK` = **3.5% of orbital speed** — is a named constant, a statement about `C_J`, and tested at both sides of the threshold. And because the moonlet is massless, `C_J` — not the system energy, which only sees the massive bodies — is the correct accuracy metric for the cell; the tests are anchored to it.

<p align="center">
  <img src="docs/flipflop_2d.gif" width="440" alt="A particle librates at L4 holding bit 1, then a super-threshold kick sends it across the separatrix and it circulates away erased">
</p>

## 3D

`orbital/nbody.py` infers its dimension from the bodies, so the same integrator runs the flat cell and an **inclined** Trojan (`demos/flipflop_3d.py`) that holds its bit while bobbing ±0.16 through the orbital plane once per orbit — genuinely three-dimensional storage, with the full 3D Jacobi integral conserved to 1e-9 (test-enforced):

<p align="center">
  <img src="docs/orbital_3d.gif" width="520" alt="Rotating 3D camera: two inclined Trojan orbits coiling through the tilted orbital plane at the L4 and L5 corners">
</p>

## Run it

```bash
pip install -r requirements.txt

python -m demos.flipflop_demo    # HOLD + NOISE MARGIN (+ Jacobi readout)
python -m demos.write_demo       # WRITE: timing diagram + the write gif
python -m demos.landscape        # the energy-landscape & anatomy figures
python -m demos.flipflop_3d      # the inclined-Trojan 3D gif
python -m demos.make_gifs        # the 2D hold->erase gif
python -m pytest                 # 51-test suite
```

## Tests (TDD, physics-validated)

51 tests check the simulation against closed-form theory, not just against itself:

- **Kepler** — a moon on a circular orbit stays circular and obeys Kepler's third law.
- **Conservation** — energy `<1e-9` and momentum in 2D & 3D; barycenter pinned; the massless moonlet exerts no back-reaction; time-reversal retraces.
- **Lagrange theory** — L4/L5 exactly equilateral and true equilibria; measured libration period matches `2π/√(27/4·μ)`; stability across mass ratios.
- **Jacobi constant** — conserved along held *and* erased orbits; held bit sits at analytic `C_L4`; erasing kicks provably lower `C_J`.
- **Rotating engine** — matches the inertial integrator trajectory-for-trajectory; blank medium circulates at the theoretical drift rate.
- **WRITE** — same blank + same pulse, timing alone selects the bit; both bits write correctly, deterministically, and survive a 40-orbit hold after an engine swap to the full N-body integrator.
- **Memory** — holds 80 (and 300, slow-marked) orbits with no secular drift; sub-threshold kicks preserve, super-threshold erase; the separatrix-crossing reader classifies tadpole/horseshoe/circulation physically.

## Layout

```
orbital/    nbody.py (2D/3D inertial integrator) · rotating.py (co-rotating frame, μ(t))
            memory.py (cell, reader, kicks) · write.py (the pinch write) · theory.py (C_J, L-points)
demos/      flipflop_demo · write_demo · landscape · flipflop_3d · make_gifs
tests/      test_nbody · test_memory · test_theory · test_rotating · test_write   (51 tests)
docs/       all figures and GIFs above
```

## Physics & numerics

- Circular restricted three-body problem, `G = 1`, total mass 1, separation 1, mean motion `n = 1`. Cell mass ratio `μ = 0.003` (< 0.0385 Gascheau/Routh limit); blank medium at `μ₀ = 3e-4`.
- Adaptive high-order Runge–Kutta (scipy `DOP853`); energy drift ~`1e-11`, Jacobi drift ~`1e-12` on held cells. A symplectic integrator (e.g. REBOUND's WHFast) is the right upgrade for Gyr-scale retention claims.

## Honest caveats & prior art

- Every mechanism here is textbook celestial mechanics: triangular Lagrange stability (Gascheau/Routh), tadpole/horseshoe co-orbitals (Janus & Epimetheus), adiabatic resonance capture (Henrard; Neishtadt; Malhotra — Neptune capturing Pluto), zero-velocity curves and the Jacobi integral. What appears unclaimed is the *construction*: engineering these into a memory cell with a write pulse, a timing diagram, a noise margin, and a test suite. The novelty is the artifact, not the mechanism.
- Wide pinch-written tadpoles librate broadly (~±70°); they hold and read robustly, but deep-cooling a written bit (shrinking its libration) needs a non-adiabatic trick — an open problem, alongside a genuine erase-to-blank cycle (write is currently one-way: μ stays grown).
- Turing-completeness is not claimed. This is a memory element; pairing it with slingshot-computing's flyby gates (flyby = logic, orbit = storage) is the longer arc.

## Roadmap

- [x] A bit that holds: L4/L5 tadpole memory (80–300 orbits, no drift)
- [x] Noise margin: separatrix threshold ~3.5% of orbital speed, restated in `C_J`
- [x] 3D: dimension-agnostic integrator + inclined-Trojan bit
- [x] Finding: memory is topological, not dissipative (drag destabilizes L4/L5)
- [x] **WRITE: adiabatic horseshoe pinch — bit selected by pulse timing alone**
- [x] Rotating-frame engine with time-dependent μ (write pulses)
- [x] Jacobi-constant theory module; 51-test physics-validated suite
- [ ] Re-writable cell: erase back to blank (shrink μ) and write again
- [ ] Bit cooling: shrink a written tadpole's libration (non-adiabatic pulse shaping)
- [ ] Averaged 1-DOF Hamiltonian: closed-form noise margin & write windows
- [ ] A register: several cells at different radii; crosstalk vs spacing
- [ ] Symplectic integrator for astronomical-timescale retention
