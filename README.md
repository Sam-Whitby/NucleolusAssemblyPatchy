# NucleolusAssemblyPatchy

VMMC lattice simulations for Samuel Whitby's PhD thesis (*"Towards a Model of Annealing in Spatial Gradients"*). Fork of [NucleolusAssembly](https://github.com/Sam-Whitby/NucleolusAssembly), replacing the isotropic Gō-model with an **orientation-dependent Kern–Frenkel patchy particle model** to prevent spurious inter-complex aggregation: every active patch on a finished complex faces an intra-complex partner, leaving the outer surface smooth.

Three simulation types:

| Binary | Geometry |
|---|---|
| `run_nucleolus` | Column (linear gradient along x, periodic y) |
| `run_condensate` | Circular (radial gradient, ring injection, open at r > R_c) |
| `run_box` | Periodic square box (no recycling; optional `--anneal` ramp) |

---

## Build

```bash
make
```

Requires C++14 (`g++`). Python 3 + `numpy` + `matplotlib` for visualisation.

---

## Target complex

A 16-particle complex of 4 polymers × 4 segments arranged on a 4×4 grid following the n=2 Moore curve. Backbone bonds connect consecutive segments within each polymer (d=1). The 4 polymers interleave to form 12 cross-polymer patch contacts in the assembled state.

```
y=3:  P0  P0  P0  P1        P0: (0,2)→(0,3)→(1,3)→(2,3)
y=2:  P0  P2  P1  P1        P1: (3,3)→(3,2)→(2,2)→(2,1)
y=1:  P2  P2  P1  P3        P2: (0,0)→(0,1)→(1,1)→(1,2)
y=0:  P2  P3  P3  P3        P3: (3,1)→(3,0)→(2,0)→(1,0)
```

---

## Interactions

### Backbone bonds (J_bb = 1000, gradient-independent)

Active at d=1 and d=√2 between consecutive same-polymer segments.

### Cross-polymer patch contacts (J = 8, gradient-scaled)

| Condition | d=1 |
|---|---|
| Patches aligned (both facing each other) | −J |
| Patches misaligned | 0 |
| d=√2 or beyond | 0 |

Same-type repulsions are absent; hard-core excluded volume (d<1) is sufficient.

### Patch slot convention (body frame)

Each particle has an orientation vector `(ox, oy)` ∈ {(1,0), (0,1), (-1,0), (0,-1)}.

| Slot | Body-frame direction |
|---|---|
| 0 | local-east (`ori` direction) |
| 1 | local-north (90° CCW) |
| 2 | local-west (opposite `ori`) |
| 3 | local-south (90° CW) |

A bond forms only when the separation vector from i to j maps to an active patch slot of i **and** the separation vector from j to i maps to an active patch slot of j.

### Reference complex energy

```
E_ref = 12 backbone bonds × (−1000) + 12 patch contacts × (−8) = −12 096
```

---

## VMMC moves

Each outer iteration performs N_particles VMMC move attempts. Three move types:

| Type | Probability | Description |
|---|---|---|
| Translation | 1 − φ_rot − φ_reorient | Rigid-body cluster translation to an adjacent lattice site |
| Cluster rotation | φ_rot | Rigid-body rotation of entire bonded cluster by 90°/180°/270° |
| In-place reorientation | φ_reorient | Single-particle Metropolis rotation: orientation changes, position fixed |

**Cluster rotation** recruits all bonded partners via VMMC link weights, preserving bond geometry. **In-place reorientation** skips partner recruitment — it is a standard single-particle move accepted by Metropolis, enabling isolated or dilute particles to explore all four patch orientations. Both move types correctly update the orientation vector.

Stokes hydrodynamic drag (`--stokes`) rejects moves with probability 1 − r_ref/r_eff (D ∝ 1/R).

---

## Column model (`run_nucleolus`)

Linear gradient γ(x) = x/L. Hard wall at x=0, open at x>L. Exited particles are recycled to x≈0. The gradient suppresses patch contacts near x=0 (denaturing) and restores them at x≥L (physiological).

### Usage

```
./run_nucleolus [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Outer iterations |
| `--snapshots N` | 1000 | Trajectory frames |
| `--length L` | 60 | Column length (lattice units) |
| `--width W` | 10 | Column width (periodic y) |
| `--gradient` | off | Enable linear gradient |
| `--stokes` | off | Stokes drag (D ∝ 1/R) |
| `--phi-rot φ` | 0.2 | Cluster rotation fraction |
| `--phi-reorient φ` | 0.2 | In-place reorientation fraction |
| `--output PREFIX` | `nucleolus` | Output file prefix |
| `--seed S` | 1 | RNG seed (0 = hardware random) |

### Example

```bash
./run_nucleolus \
    --steps 200000  --snapshots 1000 \
    --length 60     --width 10       \
    --gradient      --stokes         \
    --phi-rot 0.2   --phi-reorient 0.2 \
    --seed 2        --output my_column
```

### Visualisation

```bash
python3 visualize_nucleolus.py my_column_traj.txt --gradient-length 60 --width 10
python3 visualize_nucleolus.py my_column_traj.txt --output my_column.mp4 --fps 10
```

---

## Circular model (`run_condensate`)

Radial gradient γ(r) = γ₀ + (1−γ₀)·r/R_c. Hard wall at centre, open at r>R_c. Exited particles are returned to an injection ring at r=1 and queue for re-insertion. Three sequential phases: equilibration (g=1), denaturation (g=0), main run (gradient active).

### Usage

```
./run_condensate [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Main-phase outer iterations |
| `--snapshots N` | 1000 | Total frames across all phases |
| `--t-equil N` | 0 | Equilibration iterations (g=1) |
| `--t-denat N` | 0 | Denaturation iterations (g=0) |
| `--copies N` | 4 | Number of complex copies |
| `--radius R` | 60 | Condensate radius (lattice units) |
| `--gamma0 γ` | 0.0 | Minimum coupling at r=0 |
| `--gradient` | off | Enable radial gradient |
| `--stokes` | off | Stokes drag |
| `--coupling MODE` | `product` | `product` or `midpoint` coupling |
| `--phi-rot φ` | 0.2 | Cluster rotation fraction |
| `--phi-reorient φ` | 0.2 | In-place reorientation fraction |
| `--output PREFIX` | `condensate` | Output file prefix |
| `--seed S` | 1 | RNG seed (0 = hardware random) |

### Example

```bash
./run_condensate \
    --copies 4      --t-equil 5000   --t-denat 20000  \
    --steps 200000  --snapshots 1000 \
    --radius 60     --gamma0 0.0     --gradient        \
    --stokes        --coupling product                  \
    --phi-rot 0.2   --phi-reorient 0.2                 \
    --seed 1        --output my_circle
```

### Visualisation

```bash
python3 visualize_condensate.py my_circle_traj.txt
python3 visualize_condensate.py my_circle_traj.txt --output my_circle.mp4 --fps 10
```

---

## Output files

### `PREFIX_traj.txt`

```
<N_particles>
step=S energy=E [phase/geometry fields]
<id> <poly_type> <x> <y> <copy> <ox> <oy>
...
```

`ox`, `oy` — world-frame orientation unit vector (updated by rotation and reorientation moves).

### `PREFIX_stats.txt`

```
# step  energy  exitedParticles  exitedPerfect  acceptRatio  [phase]
```

---

## Periodic box (`run_box`)

All particles remain in the system at all times. Starts with `--copies` fully assembled target complexes. No spatial gradient — coupling is spatially uniform and set by the current gamma value.

Three phases (all optional):

| Phase | Flag | Gamma | Purpose |
|---|---|---|---|
| Equilibration | `--t-equil N` | 1 | Verify assembled complexes are stable at full coupling |
| Denaturation | `--t-denat N` | 0 | Break assembled complexes apart |
| Main | `--steps N` | 1 (or ramped if `--anneal`) | Assembly / equilibration |

With `--anneal`, gamma ramps linearly from 0 to 1 over the main-phase steps: at outer iteration i (0-indexed), gamma = i/(N-1). Without `--anneal`, gamma=1 throughout the main phase.

### Usage

```
./run_box [OPTIONS]
```

| Flag | Default | Description |
|---|---|---|
| `--steps N` | 10000 | Main-phase outer iterations |
| `--snapshots N` | 1000 | Trajectory frames across all phases |
| `--t-equil N` | 0 | Equilibration iterations at gamma=1 |
| `--t-denat N` | 0 | Denaturation iterations at gamma=0 |
| `--copies N` | 4 | Number of complex copies |
| `--box-size L` | 20 | Square box side length (lattice units) |
| `--anneal` | off | Ramp gamma 0→1 over main-phase steps |
| `--stokes` | off | Stokes drag (D ∝ 1/R) |
| `--phi-rot φ` | 0.2 | Cluster rotation fraction |
| `--phi-reorient φ` | 0.2 | In-place reorientation fraction |
| `--output PREFIX` | `box` | Output file prefix |
| `--seed S` | 1 | RNG seed (0 = hardware random) |

### Example

```bash
# Equil then anneal from denatured state
./run_box \
    --copies 4     --box-size 40      \
    --t-equil 5000 --t-denat 5000     \
    --steps 50000  --anneal --stokes  \
    --phi-rot 0.2  --phi-reorient 0.2 \
    --seed 1       --output my_box

# Fixed coupling (equilibrium at g=1)
./run_box \
    --copies 4     --box-size 30  \
    --steps 100000 --stokes       \
    --seed 1       --output equil_box
```

### Visualisation

```bash
python3 visualize_box.py my_box_traj.txt
python3 visualize_box.py my_box_traj.txt --output my_box.mp4 --fps 10
```

The left panel shows particles in the periodic square box.  The box background
colour is a uniform blue tint whose intensity tracks the current coupling
strength γ: dark blue (γ≈0, denaturing) → light blue (γ=1, fully coupled).
The top-right panel shows energy vs step; the bottom-right shows γ vs step
(useful for `--anneal` runs) with the VMMC acceptance ratio overlaid if a
stats file is present (auto-detected as `PREFIX_stats.txt`).

---

## Source layout

```
run_nucleolus.cpp          Column model driver
run_condensate.cpp         Circular model driver
run_box.cpp                Periodic box driver
visualize_nucleolus.py     Column model visualiser
visualize_condensate.py    Circular model visualiser
visualize_box.py           Periodic box visualiser
src/
  VMMC.h / VMMC.cpp        Core VMMC algorithm (cluster translation, rotation, reorientation)
  NucleolusModel.h / .cpp  Column geometry
  CondensateModel.h / .cpp Circular geometry
  Model.h / .cpp           Base class: energy, interactions, cell list
  StickySquare.h / .cpp    Kern-Frenkel patch energy on square lattice
  Box / CellList / Particle / Initialise / InputOutput
  MersenneTwister.h        MT19937 RNG
```

---

## Provenance

VMMC core: Hedges (2015), [vmmc.xyz](http://vmmc.xyz/). Extended for lattice models and Saturated Links by Holmes-Cerfon; patchy particle model and condensate geometry by Samuel Whitby.
