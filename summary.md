# NucleolusAssemblyPatchy — Implementation Summary

## Overview

This codebase is a fork of `NucleolusAssembly` (commit `af34d7e`) that replaces the existing
Gō-model interaction scheme with an **orientation-dependent, body-frame patchy particle model**.
The primary motivation is to prevent spurious aggregation between completed target complexes:
in the parent codebase, two fully assembled complexes lower their energy by 8 kT (one cross-type
native contact) when they sit adjacent, because the interaction matrix is agnostic to whether a
contact is intra-complex or inter-complex. Directional patches fix this by making each native
contact face an *inward-pointing surface patch* — in a completed complex, all patch-bearing faces
are occupied, so the outer surface is smooth and non-sticky.

---

## Repository Structure

```
NucleolusAssembly/          ← original model (github.com/Sam-Whitby/NucleolusAssembly)
NucleolusAssemblyPatchy/    ← this repo   (github.com/Sam-Whitby/NucleolusAssemblyPatchy)
```

Both are independent git repositories at the same filesystem level. Changes to one do not
affect the other.

---

## Physical Model: What Changes

### Parent model (NucleolusAssembly)

- Same-type (same polymer) particles repel at d=1 (+J), d=√2 (+J), d=2 (+εJ)
- Cross-type particles at native d=1 separation attract (-J), regardless of particle orientation
- Cross-type particles at native d=√2 separation attract (-εJ), regardless of orientation
- All interactions are isotropic in orientation

### Patchy model (this codebase)

| Interaction | Change |
|-------------|--------|
| Same-type repulsions | **Removed** — matrix entries set to zero |
| Cross-type native d=1 | **Kept, but gated by patch alignment** |
| Cross-type native d=√2 | **Removed** — set to zero (pure patch model) |
| Backbone bonds (~-1000) | **Unchanged** — robust to monomer rotation |
| Hard-core excluded volume | **Unchanged** |
| Saturated linking (phi_sl) | **Disabled** (set to 0) |
| Particle orientation | **Active** — each particle has an orientation vector that rotates with it during VMMC rotation moves |

A cross-type native d=1 bond forms **only** when both participating particles present their
active patch facing toward the other. Patches are defined in each particle's **body frame**
(local coordinates relative to its orientation vector) and rotate with the particle during
VMMC rotation moves. A completed 4×4 complex has every active patch occupied by its native
partner; exterior faces carry no patches, so two completed complexes cannot form attractive
inter-complex contacts regardless of their relative orientation.

### Physical motivation

This is the standard Kern–Frenkel patchy particle model applied to a square lattice
(Kern & Frenkel 2003, *J. Chem. Phys.* 118:9882). The original continuous-angle orientation
masking function `f = 1 iff ê_i·r̂_ij ≥ cos δ` reduces on a square lattice to the discrete
condition: the separation vector must map to an active local-frame slot on both particles
(half-opening angle δ = 45° implicitly, since on the square lattice there are exactly four
cardinal bond directions). Relative (body-frame) patches rather than global (world-frame) patches
are physically motivated because a molecular bonding surface is part of the molecule's structure
and tumbles with it in solution (Rovigatti, Russo & Romano 2018, *Eur. Phys. J. E* 41:59).

---

## Particle Orientations and the Patch Slot Convention

Each particle has a 2D unit orientation vector `(ox, oy)`, constrained to the four lattice
directions `{(1,0), (0,1), (-1,0), (0,-1)}` by VMMC's lattice-snapping after each rotation.

The four **local-frame patch slots** relative to orientation `(ox, oy)` are:

| Slot | Local direction | World direction when `ori=(1,0)` |
|------|----------------|----------------------------------|
| 0 | local-east (same as `ori`) | world-east  (+x) |
| 1 | local-north (90° CCW from `ori`) | world-north (+y) |
| 2 | local-west (opposite to `ori`) | world-west  (-x) |
| 3 | local-south (90° CW from `ori`) | world-south (-y) |

The mapping from a world-frame separation vector `(dx, dy)` to local slot is implemented in
`getPatchSlot()` (`src/StickySquare.cpp:122`):

```cpp
int lx = (int)round(dx * ori[0] + dy * ori[1]);   // local-x component
int ly = (int)round(-dx * ori[1] + dy * ori[0]);  // local-y component
```

A bond between particles 1 and 2 at separation `sep = pos1 - pos2` forms when:
- `patchSlots[id1][s1]` is true, where `s1 = getPatchSlot(ori1, -sep[0], -sep[1])` (direction *from* 1 *toward* 2)
- `patchSlots[id2][s2]` is true, where `s2 = getPatchSlot(ori2, +sep[0], +sep[1])` (direction from 2 toward 1)

---

## Patch Map

Patches are assigned to each of the 16 particle identities (local IDs 0–15) based on
which faces point toward **cross-type** (different polymer) native d=1 contact partners in
the target structure `TARGET_X/Y`.

All particles are initialised with orientation `(1,0)` (local frame = world frame at t=0),
so the patch slot for each particle is directly the world-frame direction to its native partner.

```
TARGET_X = { 0,0,1,2,  3,3,2,2,  0,0,1,1,  3,3,2,1 };
TARGET_Y = { 2,3,3,3,  3,2,2,1,  0,1,1,2,  1,0,0,0 };
```

Polymer membership: P0 = ids {0,1,2,3}, P1 = {4,5,6,7}, P2 = {8,9,10,11}, P3 = {12,13,14,15}.

| id | Pos | Polymer | Active slots | Partners (native d=1) |
|----|-----|---------|-------------|----------------------|
| 0  | (0,2) | P0 | 0(E), 3(S) | id11 via E, id9 via S |
| 1  | (0,3) | P0 | *none* | corner — no cross-type d=1 neighbours |
| 2  | (1,3) | P0 | 3(S) | id11 via S |
| 3  | (2,3) | P0 | 0(E), 3(S) | id4 via E, id6 via S |
| 4  | (3,3) | P1 | 2(W) | id3 via W |
| 5  | (3,2) | P1 | 3(S) | id12 via S |
| 6  | (2,2) | P1 | 1(N), 2(W) | id3 via N, id11 via W |
| 7  | (2,1) | P1 | 0(E), 2(W), 3(S) | id12 via E, id10 via W, id14 via S |
| 8  | (0,0) | P2 | 0(E) | id15 via E |
| 9  | (0,1) | P2 | 1(N) | id0 via N |
| 10 | (1,1) | P2 | 0(E), 3(S) | id7 via E, id15 via S |
| 11 | (1,2) | P2 | 0(E), 1(N), 2(W) | id6 via E, id2 via N, id0 via W |
| 12 | (3,1) | P3 | 1(N), 2(W) | id5 via N, id7 via W |
| 13 | (3,0) | P3 | *none* | corner — no cross-type d=1 neighbours |
| 14 | (2,0) | P3 | 1(N) | id7 via N |
| 15 | (1,0) | P3 | 1(N), 2(W) | id10 via N, id8 via W |

**Verification**: Every contact is symmetric — if id-A has a patch toward id-B, id-B has one
toward id-A. There are 12 patch-pair contacts, matching the 12 cross-type d=1 native contacts
in the reference complex. Particles id=1 and id=13 (the two outer-corner particles of P0 and P3
respectively) have no cross-type d=1 contacts in the target and carry no patches.

**Key property**: In the assembled complex (all particles at orientation (1,0)), every active
patch faces exactly one native partner, and every native partner provides the reciprocal patch.
The outer surface of the complex — all faces not pointing toward a native contact — carries no
patches and is therefore non-sticky with respect to other complexes.

---

## Interaction Matrix Changes

### `buildCouplingMatrices` (in `run_condensate.cpp` and `run_nucleolus.cpp`)

**Before (NucleolusAssembly):**
```cpp
if (sameType) {
    wD1[i][j]   = -J;         // repulsive (returns +J)
    wDsq2[i][j] = -J;         // repulsive (returns +J)
    wD2[i][j]   = -eps * J;   // repulsive (returns +eps*J)
} else {
    if (dsqd < 1.0 + TOL)   wD1[i][j]   = J;        // attractive (returns -J)
    else if (dsqd < 2.0+TOL) wDsq2[i][j] = eps * J;  // attractive (returns -eps*J)
}
```

**After (this codebase):**
```cpp
// Same-type repulsions removed: model is purely attractive at d=1 (patch-gated).
// d=sqrt(2) cross-type contacts also removed to avoid ungated inter-complex
// attraction (the patch system only gates cardinal-direction d=1 contacts).
if (!sameType && dsqd < 1.0 + TOL) {
    wD1[i][j] = J;   // attractive at d=1, gated by patch alignment
}
// All other entries remain zero.
```

This leaves `wDsq2`, `wD2`, `wDsq5` as all-zero matrices.

### Effect on backbone bond energy

Backbone pairs (same-type) previously accumulated `wD1[id1][id2] = -J` on top of `bbEnergy`:

```
old backbone energy = -(bbEnergy + (-J)) = -(1000 - 8) = -992
new backbone energy = -(bbEnergy + 0)   = -1000
```

Backbone bonds become marginally stronger (-1000 vs -992). This is intentional: the -J
same-type term was a coupling artefact of the old model, not a physical backbone stiffness.

---

## Reference Complex Energy

With the changes above, the reference energy of a perfectly assembled complex (at g=1) is:

```
E_ref = 12 × (-bbEnergy) + 12 × (-J)
      = 12 × (-1000)     + 12 × (-8)
      = -12,000 + -96
      = -12,096
```

The `referenceComplexEnergy()` function in both driver files must be updated:

```cpp
static double referenceComplexEnergy(double J, double bbEnergy)
{
    // 12 backbone bonds: -bbEnergy each (same-type wD1 is now 0)
    double E = (double)N_BB_PAIRS * -bbEnergy;
    // 12 cross-type d=1 native contacts: -J each (no d=sqrt(2) contacts in patchy model)
    for (int i = 0; i < N0; i++) {
        for (int j = i+1; j < N0; j++) {
            if (polyType(i) == polyType(j)) continue;  // skip same-type
            double dsqd = targetDistSqd(i, j);
            if (dsqd < 1.0 + 1e-6) E -= J;
            // d=sqrt(2) contacts removed from model; do not add -eps*J
        }
    }
    return E;
}
```

Note: the `eps` parameter is no longer needed in this function. The function signature can be
simplified to `referenceComplexEnergy(double J, double bbEnergy)`.

---

## Code Changes: File by File

### 1. `run_condensate.cpp`

#### a. `buildCouplingMatrices` (lines ~160–190)
Replace the matrix construction logic as shown in the section above.

#### b. Add `buildPatchSlots` function (after `buildCouplingMatrices`)

```cpp
static vector<array<bool,4>> buildPatchSlots()
{
    // For each particle identity, activate patches on faces pointing toward
    // cross-type (different polymer) native d=1 contact partners in TARGET.
    // All particles initialised with orientation (1,0), so local frame = world frame:
    //   slot 0 = world-east (+x), slot 1 = world-north (+y),
    //   slot 2 = world-west (-x), slot 3 = world-south (-y).
    vector<array<bool,4>> slots(N0);
    for (auto& s : slots) s.fill(false);

    for (int i = 0; i < N0; i++) {
        for (int j = 0; j < N0; j++) {
            if (i == j) continue;
            if (polyType(i) == polyType(j)) continue;  // same-type: no patch
            double dsqd = targetDistSqd(i, j);
            if (dsqd > 1.0 + 1e-6) continue;           // not a native d=1 contact
            // Direction from i to j in world frame
            int dx = TARGET_X[j] - TARGET_X[i];
            int dy = TARGET_Y[j] - TARGET_Y[i];
            // Map to slot (with ori=(1,0), local = world):
            // E=(1,0)→0, N=(0,1)→1, W=(-1,0)→2, S=(0,-1)→3
            int slot = -1;
            if      (dx ==  1 && dy ==  0) slot = 0;
            else if (dx ==  0 && dy ==  1) slot = 1;
            else if (dx == -1 && dy ==  0) slot = 2;
            else if (dx ==  0 && dy == -1) slot = 3;
            if (slot >= 0) slots[i][slot] = true;
        }
    }
    return slots;
}
```

#### c. Enable patches in `Interactions` setup (after `buildCouplingMatrices` call, before VMMC)

```cpp
auto patchSlots = buildPatchSlots();
interactions.patchesEnabled = true;
interactions.patchSlots = patchSlots;
```

#### d. Set `isIsotropic = false` for all particles (line ~642)

```cpp
// Change:
isIsotropic[i] = true;
// To:
isIsotropic[i] = false;   // orientation updated by VMMC rotation moves
```

#### e. Disable saturated linking

```cpp
// Change:
double phi_sl = 0.2;
// To:
double phi_sl = 0.0;      // disabled: no saturated-link moves in patchy model
```

Also remove or comment out the `--phi-sl` argument parser, or keep it as a no-op.

#### f. Update `referenceComplexEnergy`

Remove `eps` parameter; update the loop as shown in the section above.

#### g. Update the `referenceComplexEnergy` call in `main()`

```cpp
// Change:
const double refComplexEnergy = referenceComplexEnergy(J, eps, bbEnergy);
// To:
const double refComplexEnergy = referenceComplexEnergy(J, bbEnergy);
```

#### h. Remove `eps` from `Interactions` constructor call where only `wD1` is non-zero

Since `wDsq2`, `wD2`, `wDsq5` are now all-zero, they can still be passed (the
constructor accepts them) — just ensure they are populated correctly by `buildCouplingMatrices`.

### 2. `run_nucleolus.cpp`

Identical changes to items a–h above, adapted for the column model (`NucleolusModel` instead
of `CondensateModel`, `L_col` geometry, `checkAndReplace` instead of `handleExitsAndQueue`).

### 3. `src/StickySquare.cpp`

**No changes required.** The existing patch-gating logic at lines 218–231 is correct for this
model:

```cpp
if (!interactions.weakD1.empty()) {
    bool patchOk = true;
    if (interactions.patchesEnabled && !interactions.patchSlots.empty()) {
        int dx12 = (int)round(-sep[0]);
        int dy12 = (int)round(-sep[1]);
        int s1 = getPatchSlot(orientation1,  dx12,  dy12);
        int s2 = getPatchSlot(orientation2, -dx12, -dy12);
        patchOk = (s1 >= 0 && s2 >= 0 &&
                   interactions.patchSlots[id1][s1] &&
                   interactions.patchSlots[id2][s2]);
    }
    if (patchOk) energy += interactions.weakD1[id1][id2];
}
```

Since same-type `weakD1` entries are now 0, gating same-type pairs by patches is harmless (0 ×
anything = 0). All non-zero `weakD1` entries are cross-type native contacts that should be
patch-gated — the existing logic handles this correctly.

### 4. `src/CondensateModel.cpp` and `src/NucleolusModel.cpp`

**The `sameChain` guard should be reviewed.** Previously it prevented same-type repulsions
between intra-chain backbone-adjacent particles. Since same-type repulsions are removed from
`weakD1`, the guard no longer affects repulsions. However, it also prevents cross-type weak
coupling between particles in the same chain *and* the same polymer — which is vacuously the
case for backbone partners only (cross-type means different polymer, so sameChain cannot be
true for cross-type pairs). **No change is strictly required**, but a comment clarifying the
new role of `sameChain` would be helpful.

The gradient factor `g = couplingFactor(position1, position2)` still multiplies `weakD1`,
so near the condensate core (r ≈ 0, γ ≈ γ0), even aligned patches give reduced bonding —
the intended physics of the gradient model is preserved.

### 5. `src/StickySquare.h`

No changes needed. The `patchesEnabled` flag and `patchSlots` vector already exist in
`Interactions` (lines 67–68).

---

## Design Decisions and Rationale

### Why remove same-type repulsions?

Same-type repulsions were introduced to prevent same-polymer particles from stacking. In the
patchy model, the hard-core excluded volume (d < 1 → INF) already prevents overlap, and
patches ensure that attractive contacts are directional. Keeping same-type repulsions in a
purely-attractive patch model would create an asymmetric interaction (attraction via patches,
repulsion regardless of orientation) that lacks physical motivation and complicates the
analysis. The patch system makes each particle's sticky faces explicit; a non-sticky face is
simply neutral (hard-core only), not repulsive.

### Why remove d=√2 contacts?

Cross-type d=√2 contacts are **not** gateable by the current patch system — `getPatchSlot`
returns -1 for diagonal separations. If kept, they would remain attractive (−εJ = −4 kT) for
any pair of particles at d=√2 whose local IDs correspond to a native diagonal contact, regardless
of orientation. Two completed complexes can achieve d=√2 inter-complex proximity, forming
shallow attractive wells. Removing them makes the model purely patch-based at d=1 and
eliminates all thermodynamic driving force for aggregation.

### Why are backbone bonds not patch-gated?

Backbone bonds represent chain connectivity — the covalent-like sequential links along each
polymer chain. Chain integrity should be robust to monomer rotation, analogous to a glycosidic
bond allowing base rotation without breaking the sugar-phosphate backbone. The existing
direction-selective approach (bonds form only when partners are at d=1 in the correct
world-frame cardinal direction, stored in the `east`/`north` Triples arrays) is sufficient to
maintain polymer connectivity without orientation dependence.

### Why relative (body-frame) patches rather than global (world-frame)?

With global patches, every particle's sticky face always points in the same world direction
regardless of orientation. Two completed complexes would be sticky whenever their surfaces
face each other in those fixed world directions — the aggregation problem would be worse, not
better. Body-frame patches allow a complex to adopt any of the four lattice orientations
(0°, 90°, 180°, 270°) and still assemble correctly; a rotated complex has all its patches
rotated too, so they still align with native partners. This is the physically correct model:
a protein's binding surface is part of the protein and tumbles with it.

### Saturated linking (phi_sl)

The saturated-link (SL) algorithm in VMMC favours moves that change the saturation state of
strongly-bound clusters. With patches, the strongly-bound clusters are the assembled complexes
themselves. Disabling SL (phi_sl = 0) simplifies the VMMC setup and removes potential
interference between SL move selection and the patch-gating logic.

---

## Compilation

No new source files are needed. Build with the existing makefile:

```bash
cd NucleolusAssemblyPatchy
make run_condensate
make run_nucleolus
```

All changed code is in `run_condensate.cpp`, `run_nucleolus.cpp`, and the unchanged
`src/StickySquare.cpp` (which already supports patches).

---

## Verification

1. **Single assembled complex**: Run with `--copies 1 --steps 0 --snapshots 1`. Check that the
   reported system energy matches E_ref = −12,096. With patches enabled and all particles at
   orientation (1,0), all 12 patch pairs should be active.

2. **Two complexes, no aggregation**: Run with `--copies 2 --t-equil 0 --t-denat 0`. After
   equilibration, two fully-assembled complexes should diffuse independently with zero
   inter-complex attractive energy. Monitor `getSystemEnergy()` — it should equal
   2 × E_ref throughout, never dropping below 2 × E_ref + noise.

3. **Assembly from denatured state**: Run with `--copies 1 --t-denat 1000 --steps 5000`.
   Confirm that the complex reassembles after denaturation. Because orientation must align for
   patch-gated contacts, phi_rot should be non-zero (default 0.2 retained). If assembly is
   slow, increase phi_rot (see rotation_report.md for SED-motivated value of ≈ 0.43).

4. **Backbone integrity**: Confirm that no backbone bonds break during simulation (energy jumps
   of +992 or +1000 in the stats file). The existing detection code in run_condensate.cpp
   reports backbone violations; there should be none.

5. **Patch orientation test**: Place two particles (cross-type, native d=1 partners) at d=1
   in the correct world direction. With both at orientation (1,0): bond forms. Rotate one
   particle by 90°: bond does not form. This can be verified via a two-particle unit test.

---

## Key Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| J | 8.0 | Native contact energy |
| eps | 0.5 | Not used in patchy model (d=√2 contacts removed) |
| bbEnergy | 1000 | Backbone bond strength |
| E_ref | −12,096 | Reference perfect-complex energy |
| phi_sl | 0.0 | Saturated linking disabled |
| phi_rot | 0.2 | Rotation attempt probability (consider 0.43 for SED) |
| isIsotropic | false | All particles update orientation during rotation moves |

---

## References

1. Kern, N. & Frenkel, D. (2003). Fluid–fluid coexistence in colloidal systems with
   short-ranged strongly directional attraction. *J. Chem. Phys.* **118**, 9882.

2. Rovigatti, L., Russo, J. & Romano, F. (2018). How to simulate patchy particles.
   *Eur. Phys. J. E* **41**, 59.

3. Whitelam, S. & Geissler, P. L. (2007). Avoiding unphysical kinetics in Monte Carlo
   sampling. *J. Chem. Phys.* **127**, 154101.

4. Perrin, F. (1934). Mouvement brownien d'un ellipsoïde. *J. Phys. Radium* **5**, 497.

5. Doi, M. & Edwards, S. F. (1986). *Theory of Polymer Dynamics*. Oxford University Press.
