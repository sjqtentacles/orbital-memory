# Orbital Memory

**A nonvolatile memory bit made of gravity — one that holds its value the way Jupiter's Trojan asteroids have held their positions for 4.5 billion years.**

Where its sibling project [slingshot-computing](https://github.com/sjqtentacles/slingshot-computing) does *logic* with transient gravitational flybys — and is fundamentally memoryless — this one does the other half of a computer: **storage**. A bit is stored as *which stable island a body orbits in*. It's persistent, self-protecting, and requires no forces but `F = Gm₁m₂/r²`.

<p align="center">
  <img src="docs/orbital_3d.gif" width="560" alt="A rotating 3D view: two inclined Trojan orbits coiling up and down through a tilted orbital plane, parked at the L4 and L5 corners">
</p>

<p align="center"><em>An inclined Trojan bit in 3D: the particle librates at a Lagrange corner (that's the bit) while bobbing through the orbital plane once per orbit.</em></p>

---

## The idea

In the circular restricted three-body problem — a heavy **star**, a **planet**, and a massless **test particle** — there are two triangular equilibrium points, **L4** (60° ahead of the planet) and **L5** (60° behind). A particle near one of them doesn't fall in or drift away; it *librates* in a stable "tadpole" orbit around it. The two islands are separated by a **separatrix**, so a small nudge can't move the particle from one to the other.

That's a bit:

- librating around **L4** → reads **1**
- librating around **L5** → reads **0**

<p align="center">
  <img src="docs/flipflop_2d.gif" width="440" alt="Rotating frame: a particle librates at L4 (holding bit = 1), then a kick pushes it across the separatrix and it circulates away (erased)">
</p>

<p align="center"><em>The noise margin, as motion: the bit holds, until a kick past the separatrix erases it.</em></p>

## What the demo shows

```
python -m demos.flipflop_demo
```

**1 · It holds.** Both states librate steadily around ±60° for 80 orbits at an energy drift of `3e-11`. This is the thing flyby logic fundamentally cannot do — retain state. Nature's version (Jupiter's Trojans) has held for ~4.5 Gyr.

**2 · It has a noise margin.** Kicks below **~3.5% of the orbital speed** leave the bit intact (just a wider libration); above it, the particle crosses the separatrix and is ejected — erased. A real, quantifiable memory margin.

<p align="center">
  <img src="docs/flipflop.png" width="760" alt="Left: rotating-frame tadpoles at L4 and L5. Right: the resonant angle holding flat around plus and minus 60 degrees for 80 orbits.">
</p>

## The physics finding

The memory is protected **topologically, not dissipatively.** L4/L5 are *maxima* of the effective potential — they're stable only because of the velocity-dependent Coriolis force. Add friction and you kill that stabilization: **drag destabilizes the Lagrange points** (verified in code). So there's no attractor pulling the bit back to a clean value; instead the bit is held by the *geometry of phase space* — an invariant island bounded by a separatrix, the KAM / Nekhoroshev way. Refresh isn't damping; it's the topology.

A consequence: **writing** (flipping 0↔1) is genuinely subtle — you can't just drag a particle into the island, because drag repels it. Real Trojan capture happens through slow orbital migration and adiabatic resonance capture (Henrard, Neishtadt). That's the honest frontier of this project, and the natural place for a proper Hamiltonian-theory pass.

## 3D

The dynamics generalize cleanly: `orbital/nbody.py` infers its dimension from the bodies, so the same integrator runs the flat bit and an **inclined** Trojan (`demos/flipflop_3d.py`), which holds the L4/L5 bit in-plane *and* oscillates ±0.16 through the orbital plane once per orbit — real out-of-plane motion, drift `1.8e-11`. The rotating-camera GIF at the top is that, not a projected 2D scene.

## Run it

```bash
pip install -r requirements.txt

python -m demos.flipflop_demo   # HOLD + NOISE MARGIN, writes out/flipflop.{png,json}
python -m demos.flipflop_3d     # inclined Trojan, writes docs/orbital_3d.gif
python -m demos.make_gifs       # the 2D hold->erase gif
```

## Layout

```
orbital/    nbody.py (2D/3D integrator, G=1) · memory.py (L4/L5 cell, readout, kicks)
demos/      flipflop_demo.py · flipflop_3d.py · make_gifs.py
docs/       the figures and GIFs above
```

## Physics & numerics

- Circular restricted three-body problem, `G = 1`, total mass 1, separation 1,
  mean motion `n = 1` (period `2π`). Planet mass fraction `μ = 0.003` (`< 0.0385`,
  so L4/L5 are linearly stable).
- Adaptive high-order Runge–Kutta (scipy `DOP853`, `rtol ~ 1e-11`); energy drift
  ~`1e-11` over the runs. A symplectic integrator (e.g. REBOUND's WHFast) is the
  right upgrade for retention claims over astronomical timescales.

## Honest caveats & prior art

- The stored bit rides on real, well-studied physics: triangular Lagrange
  stability (Gascheau/Routh), tadpole vs horseshoe libration, resonance capture
  (Malhotra; Neptune capturing Pluto into 3:2), KAM/Nekhoroshev stability of
  invariant tori. What appears unclaimed — as with the slingshot project — is
  engineering it into a *built memory artifact* with a read, a noise margin, and
  a route to a write. The novelty is the construction, not the mechanism.
- Turing-completeness is not claimed. This is a memory element; pairing it with
  slingshot-computing's gates (flyby = write head, orbit = storage) is the
  longer arc.

## Roadmap

- [x] A bit that holds: L4/L5 tadpole memory, read + retention (80 orbits, 3e-11)
- [x] Noise margin: the separatrix threshold (~3.5% of orbital speed)
- [x] 3D: dimension-agnostic integrator + inclined-Trojan visualization
- [x] Finding: the memory is topological, not dissipative (drag destabilizes L4/L5)
- [ ] Write: adiabatic resonance capture (migration-driven), not drag
- [ ] Theory pass: averaged Trojan Hamiltonian → closed-form noise margin & retention
- [ ] A multi-cell "byte"; the inclination bob as an analog sub-register
- [ ] Symplectic integrator for astronomical-timescale retention
