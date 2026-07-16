# Orbital Memory

[![tests](https://github.com/sjqtentacles/orbital-memory/actions/workflows/tests.yml/badge.svg)](https://github.com/sjqtentacles/orbital-memory/actions/workflows/tests.yml)

**A nonvolatile memory cell made of gravity — written by orbit insertion, held by nothing but `F = Gm₁m₂/r²`, read by watching an angle, and rewritten by a real gravitational flyby — storing its bit the way Jupiter's ~13,000 Trojan asteroids store theirs.** No fake actuators: every operation maps to something that actually happens in space, and the simulation is validated against the real sky.

Where its sibling project [slingshot-computing](https://github.com/sjqtentacles/slingshot-computing) does *logic* with transient gravitational flybys — and is fundamentally memoryless — this one is the other half of a computer: **storage**. A bit is stored as *which stable island a body librates in* — L4 (leading) or L5 (trailing).

<p align="center">
  <img src="docs/validation.png" width="760" alt="Left: the resonant angle librating with a 148-year period in real years. Right: the tadpole around L4 in AU. Newtonian gravity at the real Sun-Jupiter ratio reproduces the observed Trojan libration.">
</p>

<p align="center"><em><b>Validated against the sky.</b> Run the cell at the real Sun–Jupiter mass ratio and convert the clock to years: the simulated bit librates with a <b>~148-year period</b> — the period actually observed for Jupiter's Trojans. Nothing is fit to that number; it falls out of Newtonian gravity at the real ratio.</em></p>

---

## Is this real? — the honest boundary

An earlier version of this project wrote and erased bits by *growing and shrinking the planet's mass*. That was the one fake part — planets don't do that on command — and it's gone. Every operation now maps to a real space mechanism:

| Operation | Realistic mechanism | Real-world analog |
|---|---|---|
| **WRITE** | orbit **insertion** — deliver the body, one insertion burn (Δv) drops it onto the tadpole | stationing a spacecraft at L4/L5; natural Trojan capture |
| **HOLD** | pure gravity, no forces added | Jupiter's Trojans, held for ~4 Gyr |
| **READ** | measure the libration angle | astrometry |
| **COOL** | **station-keeping** burns, a few m/s of Δv | how real co-orbital spacecraft would hold station |
| **REWRITE** | **erase by flyby, then re-insert** — both halves real | co-orbitals really are scattered and re-injected by encounters |
| **ERASE** | an aimed massive **flyby** scatters the bit out of resonance | gravitational scattering |

The simulation is honest Newtonian *n*-body throughout (scipy `DOP853`, no softening), validated against the observed Jupiter-Trojan libration period, Jupiter's real orbital speed (13.06 km/s), and its real year (11.87 yr). **The one remaining caveat is the obvious one: nobody is going to build a memory device out of asteroids.** It is a faithful simulation of real orbital dynamics and a physics/art project — not a product. The physics is real; the application is a conceit, and it says so.

## The full memory cycle

| Operation | Mechanism | Physics |
|---|---|---|
| **HOLD** | tadpole libration around L4/L5 | topological protection (invariant island, KAM) |
| **WRITE** | orbit insertion — one burn onto the tadpole | targeting / orbit insertion (Δv reported in m/s) |
| **COOL** | tangential station-keeping burns | co-orbital pendulum damping |
| **READ** | which side the resonant angle librates on | the separatrix-crossing classifier |
| **REWRITE** | erase-by-flyby, then re-insert | gravitational scattering + insertion |
| **ERASE** | an aimed massive flyby scatters the bit out | separatrix crossing / logic acting on memory |

- librating around **L4** (60° ahead of the secondary) → reads **1**
- librating around **L5** (60° behind) → reads **0**
- horseshoe / circulation → **erased / blank**

**This is moon-scale hardware too.** Saturn's moons **Telesto** and **Calypso** ride Tethys's L4/L5, **Helene** and **Polydeuces** ride Dione's, and **Janus & Epimetheus** live on the horseshoe orbits the blank medium uses. The dynamics depend only on the mass ratio `μ`, so the same code covers star+planet, planet+moon, and moon+moonlet.

## The real swarm

Real Jupiter carries two clouds of Trojans — the **Greeks** leading at L4, the **Trojans** trailing at L5. Seed a cloud of massless particles at the real ratio and advance them together, and the same two lobes appear:

<p align="center">
  <img src="docs/swarm.gif" width="520" alt="Two clouds of particles librating 60 degrees ahead of and behind Jupiter — the Greeks at L4 and Trojans at L5 — on the ring at 5.2 AU">
</p>

## WRITE: placing a bit

Writing is *delivering the body into the chosen island* — exactly how you would station a spacecraft at L4/L5. The body arrives on a co-orbital transfer (on which, coasting, it is **not** a bit — it reads erased) and a single **insertion burn** captures it onto the tadpole; the cost is a real, *actually-applied* Δv (~510 m/s for a deep L4 bit at the Sun–Jupiter scale):

<p align="center">
  <img src="docs/insert.gif" width="440" alt="A body coasts in on a transfer orbit, an insertion burn flashes, and it settles into a librating tadpole at L4 — bit 1 written by a real applied delta-v">
</p>

Insertion into L4 writes **1**, into L5 writes **0**. No planet is grown; the burn is a genuine velocity change (`dv = |v_written − v_arrival|`, test-enforced), and it is *load-bearing* — the same body without the burn coasts off the island and reads erased.

## The phase portrait — why a bit is a bit

The stored bit is a slow pendulum whose coordinate is the resonant angle `φ` and whose momentum is the radial offset `da = r − 1`. The two tadpole islands (around L4 and L5) are the two memory states; the **separatrix** between them is why a small perturbation can't flip the bit:

<p align="center">
  <img src="docs/phase_portrait.png" width="620" alt="Phase portrait: nested tadpole loops around L4 (bit 1) and L5 (bit 0), the horseshoe band, and circulation">
</p>

The same four families, drawn as orbits in the rotating frame — two tadpole bits, the horseshoe (erased), and free circulation (blank):

<p align="center">
  <img src="docs/anatomy.png" width="520" alt="Rotating-frame orbit families: cyan tadpole at L4, pink tadpole at L5, amber horseshoe sweeping the ring, dashed grey circulation">
</p>

## Cooling: deeper, harder bits

A freshly placed bit can librate wide. `orbital/cool.py` tightens it the way a real co-orbital spacecraft would — **tangential station-keeping burns**, a few m/s each (capped at `0.008`; the erase threshold is `0.035`). The tadpole is a slow pendulum whose momentum is `da = r − 1`; burns at mid-swing damp it, taking **±62° to under 30° in a couple of burns**:

<p align="center">
  <img src="docs/cool.gif" width="440" alt="A wide tadpole tightens around L4 under flash-marked station-keeping burns, ending as a deep bit that still reads 1">
</p>

Getting here took three wrong schemes, kept in the module docstring as physics documentation (retro-kicks eject through the L1 neck; 'raise C_J' is exactly backwards; blind prograde kicks fall out the bottom of the band). The surviving lesson: **C_J stratifies but does not classify** — a cooled, slightly eccentric deep tadpole sits *below* C_L4 while erased orbits sit nearby, so the operational tests (amplitude, readback, honest-engine hold) carry the correctness, not the Jacobi value. A deep bit is also a **harder** bit: farther from every separatrix, it takes a closer flyby to erase.

## The gate, and rewriting by flyby

The unification of the two projects: slingshot-computing's mechanism — an aimed flyby, its launch direction root-found on the full simulation — pointed at this project's stored bit. A massive bullet (`m = 2e-4`) on a fast hyperbolic pass shaves the moonlet at closest approach `0.002`:

<p align="center">
  <img src="docs/gate.gif" width="640" alt="Two synchronized lanes: with the bullet the moonlet is flung off its island and reads ERASED; without it the bit keeps librating">
</p>

**Bullet present → bit erased. Bullet absent → bit survives. A graze at 25× the distance → bit survives** (locality). The bullet's presence is a logic input; the stored bit is the register. Two findings run it: the **guiding center stores the bit, and only tangential impulse moves it** (a radial pass at miss 0.004 pumps a monster epicycle yet leaves the bit readable), and a massive intruder **drags the system barycenter**, so readout uses a COM-corrected resonant angle.

**Rewrite** is then real, and reliable: **erase the old bit with a flyby, then insert the new one** — `gate.rewrite_cycle("1", "0")` reads back `0`. Both halves are individually tested, so the cell is genuinely rewritable. The tempting one-shot — a single flyby that *flips* L4→L5 — turns out to be impossible: sweeping the pass depth, the bit's amplitude pumps but stays in its own island until, past a threshold, it erases outright; no depth lands it in the other island (a pinned finding). One conservative impulse can knock the guiding center out of a tadpole, but not settle it into the opposite one.

## Eccentric orbits: a breathing bit

Real planets are eccentric — Jupiter's is `e ≈ 0.0489`. In the elliptic restricted problem the primaries breathe in and out once per year (the equilateral point is a central configuration, so L4 is a fixed linear image of the star→planet vector, seeded exactly). The bit survives, its libration modulating at the orbital frequency while the separation breathes exactly `1 ± e` (test-enforced). There is no exact Jacobi integral here, so the tests assert the honest invariants — retention and bounded, orbit-locked breathing — not a conserved scalar.

## Saturn: the onset of erosion

Add real Saturn as a fourth body (real mass ratio, real 1.84× spacing) and it periodically tugs the Trojan — the mechanism that sculpts the real swarm. Full erosion is a **Gyr** process (secular resonances), far beyond a feasible run; what a short integration shows is its *onset*, and that is all the panels claim:

<p align="center">
  <img src="docs/erosion.gif" width="640" alt="Two lanes: with Saturn the Trojan's libration pumps wider over 1500 years; without Saturn the control bit sits nearly still">
</p>

With Saturn the libration amplitude pumps to ~2× the Saturn-free control over the run (a relative, platform-invariant assertion; the exact figure is chaotic).

## Where a bit survives

The triangular points are linearly stable only for `μ < 0.0385` (the Gascheau/Routh limit); beyond it L4/L5 come apart. Sweep `(μ, amplitude)` and integrate each cell to see where a bit lives — every real co-orbital sits far on the stable side:

<p align="center">
  <img src="docs/stability.png" width="640" alt="Stability map in mass ratio vs libration amplitude: a green stable region up to the Routh limit, red beyond, with Jupiter/Earth/Saturn co-orbitals marked deep in the stable zone">
</p>

## The energy landscape

<p align="center">
  <img src="docs/landscape.png" width="560" alt="The rotating-frame effective potential with zero-velocity curves, the five Lagrange points, and the two stored tadpole orbits">
</p>

The rigorous backbone (circular case) is the **Jacobi constant** `C_J`. A held bit sits at the exact triangular value `C_L4 = 3 − μ(1−μ)` (sim matches to 2e-4; no secular drift, test-enforced); a kick lowers `C_J` toward the separatrix; past it, the bit erases. The noise margin — kicks below `memory.ERASE_KICK` = **3.5% of orbital speed** — is a named constant, a statement about `C_J`, and tested on both sides. Because the moonlet is massless, `C_J` — not the system energy, which only sees the massive bodies — is the correct accuracy metric for the cell.

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

python -m demos.validation       # the money figure: sim vs observed 148-yr Trojan period
python -m demos.insert_demo      # WRITE: orbit insertion (Δv in m/s)
python -m demos.swarm            # the real Greek & Trojan clouds
python -m demos.phase_portrait   # tadpoles / horseshoe / circulation
python -m demos.stability_map    # where a bit survives (Routh limit + real objects)
python -m demos.gate_demo        # conditional erase by aimed flyby
python -m demos.cool_demo        # station-keeping burns tighten a bit
python -m demos.saturn_erosion   # Saturn pumping the Trojan (the onset)
python -m demos.landscape        # the energy-landscape & anatomy figures
python -m demos.flipflop_3d      # the inclined-Trojan 3D gif
python -m demos.make_gifs        # the 2D hold->erase gif
python -m pytest                 # 121-test suite
```

## Tests (TDD, physics-validated)

**121 tests** check the simulation against closed-form theory *and against real observations*, not just against itself:

- **Validation against the sky** — at the real Sun–Jupiter ratio the simulated libration period converts to **~148 years** (the observed Trojan value), matches linear theory, and the bit is stable with no secular drift; the Earth-Trojan case is documented as the honest harder case (linear theory underestimates 2010 TK7's large-amplitude libration).
- **Real units** — `orbital/units.py` reproduces Jupiter's 11.86-yr period and 13.06 km/s orbital speed; round-trip conversions exact.
- **Kepler & conservation** — circular stays circular, Kepler's third law; energy/momentum in 2D & 3D; barycenter pinned; massless moonlet exerts no back-reaction; time-reversal retraces.
- **Lagrange & Jacobi** — L4/L5 exactly equilateral; measured libration period matches `2π/√(27/4·μ)`; `C_J` conserved along held and erased orbits; held bit at analytic `C_L4`; erasing kicks provably lower `C_J`.
- **WRITE (insertion)** — insert '1'→L4, '0'→L5; achieved amplitude matches request; insertion Δv finite and in a plausible m/s range; deterministic; survives a 45-orbit hold with no secular drift.
- **REWRITE / GATE** — `rewrite_cycle` reads back the new bit both directions; aim converges; present erases / absent survives / 25× graze survives; COM-corrected readout provably differs; the single-flyby "flip" boundary is pinned (never reaches the other island).
- **COOL** — burns shrink a wide bit below target with the value preserved; bounded ≪ threshold; survives the honest engine; deterministic.
- **Eccentric** — reduces to circular at e=0; the bit survives real Jupiter eccentricity at L4 and L5; separation breathes exactly `1 ± e` at the orbital frequency.
- **Saturn perturber** (slow) — Saturn measurably pumps the Trojan vs a Saturn-free control (relative, platform-invariant).
- **Memory reader** — the separatrix-crossing classifier handles tadpole/horseshoe/circulation physically, wide tadpoles included; drag destabilization is a committed test.

## Layout

```
orbital/    nbody.py (2D/3D inertial) · rotating.py (co-rotating frame, analysis)
            memory.py (cell, reader, kicks, eccentric, Saturn) · units.py (real units)
            write.py (orbit insertion) · cool.py (station-keeping) · gate.py (flyby + rewrite)
            theory.py (C_J, Lagrange points)
demos/      validation · insert_demo · swarm · phase_portrait · stability_map
            gate_demo · cool_demo · saturn_erosion · landscape · flipflop_3d
            make_gifs · flipflop_demo · style.py (shared visual system)
tests/      121 tests: real-sky validation, physics invariants, every operation
docs/       all figures and GIFs above
```

## Physics & numerics

- Circular restricted three-body problem, canonical units (`G = 1`, total mass 1, separation 1, mean motion `n = 1`); a thin units layer rescales to real kg / AU / years / m·s⁻¹. Cell mass ratio `μ = 0.003` for the fast demos (< 0.0385 Gascheau/Routh limit); validation runs at the real `μ = 9.54e-4`.
- Adaptive high-order Runge–Kutta (scipy `DOP853`), no softening; Jacobi drift ~`1e-12` on held cells. A symplectic integrator (e.g. REBOUND's WHFast) is the right upgrade for true Gyr-scale retention claims.

## Honest caveats & prior art

- Every mechanism here is textbook celestial mechanics: triangular Lagrange stability (Gascheau/Routh), tadpole/horseshoe co-orbitals (Janus & Epimetheus), the elliptic-restricted equilateral solution, Jacobi integral and zero-velocity curves, orbit insertion, gravitational scattering. What is unclaimed is the *construction*: engineering these into a memory cell with an insertion write, a station-keeping cool, a flyby rewrite, a noise margin, and a real-units validation. The novelty is the artifact, not the mechanism.
- Nobody will build this out of asteroids — it is a simulation, validated against real dynamics, not a device. The write actuator (an insertion burn) and rewrite (an aimed flyby) are physically real but wildly impractical to deploy at solar-system scale. That's the whole conceit, stated plainly.
- Turing-completeness is not claimed. This is a memory element; pairing it with slingshot-computing's flyby gates (flyby = logic, orbit = storage) is the longer arc.

## Roadmap

- [x] A bit that holds: L4/L5 tadpole memory (80–300 orbits, no drift)
- [x] Noise margin: `ERASE_KICK` = 3.5% of orbital speed, tested both sides
- [x] 3D: dimension-agnostic integrator + inclined-Trojan bit (full 3D Jacobi integral)
- [x] Finding: memory is topological, not dissipative (drag destabilizes L4/L5 — tested)
- [x] **Real units + validation: the simulated bit librates at the observed ~148-yr Trojan period**
- [x] **WRITE by orbit insertion — a real Δv, no grown planet**
- [x] **COOL: honest station-keeping burns**
- [x] **REWRITE by real physics: erase-by-flyby then re-insert (single-flyby flip pinned impossible)**
- [x] **GATE: conditional erase by aimed flyby — logic acts on memory (the capstone)**
- [x] **Eccentric (elliptic restricted) orbits — a breathing bit survives real Jupiter eccentricity**
- [x] **Real Saturn perturber — the onset of swarm erosion**
- [x] CI; 121-test physics-and-observation-validated suite; findings pinned as tests
- [ ] Dual-rail write head: two moonlets, a routed bullet erases one — the surviving rail IS the written bit (locality datum already tested)
- [ ] Averaged 1-DOF Hamiltonian: closed-form noise margin & write windows
- [ ] A register: several cells at different radii; crosstalk vs spacing
- [ ] Symplectic integrator (REBOUND/WHFast) for astronomical-timescale retention
```
