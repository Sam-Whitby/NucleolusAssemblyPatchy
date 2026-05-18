# Better Lattice Models for Ribosomal Assembly

*Design notes — S. Whitby, 2026*

---

## 1. The Problem with the Current Model

The existing 4×4 Moore-curve complex has four equal-length polymers that assemble with twelve cross-polymer patch contacts distributed roughly symmetrically. Assembly proceeds by any of the four polymers arriving in any order — the landscape has approximate four-fold symmetry, and the kinetic frustration that arises is *radial* (many partially assembled states with one or two polymers attached, no clear directionality to the pathway).

This does not resemble ribosome assembly. The real large and small ribosomal subunits are built around a single large rRNA scaffold that must fold before ribosomal proteins (r-proteins) can bind. Assembly is deeply *hierarchical*: a set of six or seven primary r-proteins binds directly to the folded rRNA; each primary r-protein creates the binding site for several secondary r-proteins; and so on in three or four tiers. The assembly pathway is therefore directed — there are dominant early events (RNA folding, primary protein binding) and late events (tertiary protein binding) — rather than the symmetric, interchangeable assembly of the current model.

Furthermore, the King et al. (2024) pH data introduces a physical asymmetry between contact types that the current model cannot exploit:

- **RNA-RNA contacts** (intra-RNA tertiary structure: base-stacking, triplex contacts, A-minor interactions) are **strengthened** at lower pH. This is because protonation of adenine (pKa ~3.5 but shifted to ~6 in structured contexts) and cytosine (pKa ~4.2, similarly shifted) stabilises hydrogen bonds and stacking interactions, while partial neutralisation of phosphate repulsion improves backbone packing. Base pairs and tertiary contacts in rRNA are consistently more stable in slightly acidic conditions (pH 6.0–6.5 vs. 7.2).

- **RNA-protein contacts** are less effective at the condensate interior, but through indirect mechanisms rather than direct charge-state changes. The dominant RNA-binding residues in r-proteins — Arg (pKa ≈ 12.5) and Lys (pKa ≈ 10.5) — account for ~75% of all r-protein:rRNA contacts and are fully protonated throughout the 6.5–7.2 range: a pH-independent contribution. RNA phosphodiester groups (pKa ~1–2) are also fully charged at pH 6.5, so the negative charge density on the rRNA backbone is unchanged. The two mechanisms that do reduce effective RNA-protein coupling at lower pH are: (i) **assembly factor impairment** — DDX21 and related helicases have pH optima of 6.8–7.5 (Rai et al. 2023) and are sub-optimal at pH 6.5, reducing the rate at which RNA kinetic traps are resolved and correct protein binding sites are presented; and (ii) **increased RNA secondary structure stability** at lower pH (first bullet) occludes protein binding sites that require partial unfolding of RNA for access. Both are kinetic effects that reduce the effective propensity for RNA-protein contacts to form in the interior.

This distinction is the physical basis for a **dual-gradient model**: the inside of the condensate (low pH) strengthens RNA-RNA contacts while simultaneously weakening RNA-protein contacts. Moving outward (toward physiological pH), RNA-RNA contacts relax to physiological values while RNA-protein contacts are allowed to form stably. The spatial gradient thus implements a temporally sequential assembly pathway: RNA folds first (interior), proteins bind second (exterior). This is the observed in vivo assembly sequence.

The models proposed here are designed to exploit this dual-gradient mechanism.

---

## 2. Central Design Principle: The Dual-Gradient

Define two distinct coupling parameters mapped to two types of native contact:

**γ₁(r)** — coupling for RNA-RNA contacts (intra-RNA tertiary structure):
```
γ₁(r) = 1 − (1−γ₀) × r/R_c     [inverted gradient, highest at centre]
```
This represents the pH effect: lower pH (interior) strengthens RNA-RNA contacts. γ₁ decreases from 1 at r=0 to γ₀ at r=R_c.

**γ₂(r)** — coupling for RNA-protein contacts (r-protein binding):
```
γ₂(r) = γ₀ + (1−γ₀) × r/R_c    [same as current model, highest at periphery]
```
This represents the pH effect: lower pH (interior) weakens RNA-protein contacts. γ₂ increases from γ₀ at r=0 to 1 at r=R_c.

**The key competition**: Near the centre (low pH), RNA-RNA contacts are near-maximal but RNA-protein contacts are suppressed. This promotes RNA-RNA structure formation and prevents premature lock-in by misplaced proteins. Near the periphery (physiological pH), both contact types are strong but RNA-RNA contacts are slightly weaker than the maximum, allowing the RNA to reorganise while proteins stabilise the correct fold. This is a spatial implementation of the annealing principle with a built-in hierarchy.

With a single γ₀ parameter, the entire two-channel gradient is specified. The cross-over radius where γ₁ = γ₂ = (1+γ₀)/2 is at r = R_c/2 — the midpoint of the condensate, roughly the FC/DFC to GC boundary.

In the code, this requires separating the weakD1 matrix into two sub-matrices (weakD1_intraRNA and weakD1_interRNA_prot) and applying different coupling functions to each:

```cpp
energy += couplingFactor_1(pos1, pos2) * weakD1_intra[id1][id2];   // RNA-RNA
energy += couplingFactor_2(pos1, pos2) * weakD1_inter[id1][id2];   // RNA-protein
```

This is a minor code change to CondensateModel.cpp.

---

## 3. Model Taxonomy

Four models are presented in increasing order of biological realism and computational cost. All share the dual-gradient architecture.

| Model | Dimension | Particles/complex | Chains | Description |
|---|---|---|---|---|
| **S** (Simple) | 2D | 16 | 1 | Single RNA chain, 4×4 grid, two contact classes by sequence separation |
| **A** (Asymmetric) | 2D | 22 | 5 | 1 RNA chain (16) + 4 protein dimers (8); explicit multi-chain hierarchy |
| **B** (Better) | 3D | 27 | 1 | Single RNA chain, 3×3×3 cube, cubic-lattice Hilbert path |
| **C** (Complete) | 3D | 39 | 7 | 1 RNA chain (27) + 6 protein dimers (12); full dual-gradient hierarchy |

Model **A** and Model **C** are the recommended choices for publication (A for initial demonstration, C for the full result). Model **S** is a drop-in replacement for the current code. Model **B** provides the 3D single-chain intermediate.

---

## 4. Model S: 2D Single-Chain Refactoring

### 4.1 Motivation

The most immediate improvement to the current code with zero new code beyond re-labelling: treat the 16-particle Moore-curve complex as a *single RNA chain* rather than four equal polymers, and classify contacts by sequence separation rather than by polymer identity.

### 4.2 Construction

Keep the existing 4×4 = 16-particle target complex with the same Moore-curve TARGET_X/Y coordinates. Instead of four polymers of four segments each, reinterpret as:

**One polymer of 16 segments**, ordered by the Moore-curve path:
```
Segment order: 0→1→2→3→4→5→6→7→8→9→10→11→12→13→14→15
```
Using the existing path ordering within the Moore curve, consecutive segments are already at d=1, so backbone bonds are unchanged. The only change is the identification of contact *types*:

**Local contacts** (|i−j| ≤ 4): segments close in the chain — represent RNA secondary structure and local base-stacking. Apply γ₁ (inverted, pH-enhanced).

**Global contacts** (|i−j| > 8): segments far apart in the chain that come together in the native 3D fold — represent long-range tertiary contacts mediated by r-proteins (implicit protein model). Apply γ₂ (pH-suppressed).

**Intermediate contacts** (4 < |i−j| ≤ 8): middle-range contacts; apply γ₁ or γ₂ depending on the specific geometry (tertiary structure type). A simple default is to use γ₂ for all non-local contacts.

### 4.3 Implicit Protein Model

This is the "neat" option: proteins are represented implicitly as the long-range bridging contacts in the RNA chain. The physical picture is that a pair of RNA segments at positions i and j (|i−j| > 8), whose native contact requires an r-protein to bridge them, has its effective coupling scaled by γ₂ (pH-suppressed). The protein is never explicitly in the simulation — it is encoded in the energy function. This reduces the parameter count to:

| Parameter | Meaning |
|---|---|
| J_bb | Backbone bond strength (~1000) |
| J₁ | Local (RNA-RNA) contact strength |
| J₂ | Global (RNA-protein-RNA bridge) contact strength |
| γ₀ | Coupling floor for both channels |
| cutoff k | Sequence separation threshold (default k=4 local, k=8 global) |

Five parameters. With J₁ = J₂ = J (same energy scale, just different γ), this reduces to four parameters total.

### 4.4 Patch Assignments

Patches are assigned exactly as in the current model: each segment has active patch slots toward its native d=1 cross-type (cross-position) partners. The patch system is unchanged — the only change is which of the two γ values multiplies the interaction energy.

### 4.5 Assembly Pathway

The inverted γ₁ gradient means that in the condensate interior, local RNA-RNA contacts form readily, producing a partially folded RNA chain that may contain locally correct but globally misoriented secondary structure. As the chain moves outward:
- Local contacts (γ₁) weaken slightly — misfolded local structures can correct
- Global contacts (γ₂) strengthen — long-range bridging (protein binding) is now allowed

The pathway is: *local RNA structure in the interior → global tertiary contacts at the periphery*. This is the observed ribosome assembly sequence.

---

## 5. Model A: 2D Multi-Chain with Explicit Proteins

### 5.1 Geometry

**RNA chain** (chain 0, segments 0–15): the existing 4×4 Moore-curve arrangement on the 4×4 grid. This is the rRNA scaffold.

**Protein chains** (chains 1–4, 4 dimeric proteins, segments 16–23): each protein occupies 2 adjacent lattice sites adjacent to the outer surface of the RNA chain. Four protein binding sites are chosen at four geometrically distinct locations on the RNA chain perimeter.

The protein binding sites are chosen to require the RNA chain to be in its native folded state before the protein can bind — i.e., the protein's patch faces a site on the RNA that is only exposed in the correctly folded structure. This creates genuine hierarchy: proteins cannot bind until the RNA has folded.

### 5.2 Target Complex

The assembled complex occupies a 4×4 grid (RNA) plus 4 adjacent protein pairs:

```
Protein 1         Protein 2
    [P1a P1b]     [P2a P2b]
[P0a][RNA 4×4 grid][RNA   ]
[P0b][RNA         ][RNA   ]
    [P3a P3b]     [P4a P4b]
```

(Exact positions depend on which RNA surface sites are chosen for protein binding; see below.)

**Protein binding sites**: choose 4 RNA segments on the outer perimeter whose d=1 sites face outward from the assembled complex. These are the protein docking sites. The protein dimer (two segments, backbone bond) has one segment docking on the RNA surface site and one segment extending outward.

### 5.3 Contact Classification

| Contact type | Pairs | γ | Physical meaning |
|---|---|---|---|
| RNA backbone | consecutive (i, i+1) | gradient-free | covalent backbone (always strong) |
| RNA-RNA local | |i-j| ≤ 4, d=1 native | γ₁ (inverted) | base-stacking, secondary structure |
| RNA-RNA global | |i-j| > 8, d=1 native | γ₂ (normal) | cross-domain tertiary contacts |
| RNA-protein | RNA segment ↔ protein segment, d=1 native | γ₂ (normal) | r-protein binding |
| Protein backbone | consecutive protein segments | gradient-free | protein chain connectivity |
| Protein-protein | between proteins at d=1 in native structure | γ₂ (normal) | r-protein–r-protein contacts |

### 5.4 Why This Is Better

1. **Explicit hierarchy**: The RNA chain must adopt a folded state before the four protein binding sites become accessible. Simulations starting from a denatured state require the RNA to fold first — primary assembly event — before proteins can stably bind.

2. **The dual-gradient does real work**: In the condensate interior (low pH), γ₁ is high (RNA-RNA contacts are strong) but γ₂ is low (protein binding is suppressed). The RNA chain will fold — possibly with wrong local contacts — but proteins cannot lock in the misfolded state. As the particle moves outward, protein binding becomes possible and selectively stabilises the correctly folded RNA. This is a genuine selectivity mechanism that the symmetric 4-polymer model cannot produce.

3. **Realistic frustrated intermediates**: Misassembled states include: (a) RNA locally folded but with wrong global topology; (b) one or two proteins bound to a partially folded RNA, locking in a metastable state; (c) two RNA chains tangled around each other with wrong inter-chain contacts. These are all qualitatively realistic failure modes.

4. **Generalisation**: By changing the number and size of protein chains, the same framework describes any RNA-protein complex: snRNPs (1 RNA + 7 proteins), SRP (1 RNA + 6 proteins), RNase P, telomerase.

### 5.5 Backbone Flexibility

**RNA chain backbone** (J_bb = 20–50 kBT): significantly lower than the current 1000, allowing rare but possible backbone fluctuations. This is motivated by the fact that RNA chains in the cell undergo conformational changes (kink-turn motifs, A-minor interactions) that require local backbone reorientation. A backbone of 20 kBT means backbone bonds break at equilibrium with probability ~e⁻²⁰ ≈ 10⁻⁹ per step — essentially never broken, but the softer restoring force allows the chain to adopt non-extended configurations more naturally.

**Protein chain backbone** (J_bb_prot = 50 kBT): proteins are slightly stiffer than the RNA (protein secondary structure is more rigid than RNA backbone conformational isomers).

Note: on the current lattice, "backbone bond" controls which VMMC cluster-building recognises chain connectivity. The value of J_bb just needs to be much greater than J to ensure chains are never broken at equilibrium. J_bb = 20 × J is sufficient. Using J_bb = 1000 when J = 8 is a factor of 125×, which is excessive and may slow exploration of conformational space by making energy barriers to chain rearrangement very high. J_bb = 20 × J = 160 (if J = 8) is physically cleaner and easier to justify.

---

## 6. Model B: 3D Single-Chain on a Cubic Lattice

### 6.1 Motivation

Moving to 3D is the most important physical improvement. In 3D, each particle has 6 faces (±x, ±y, ±z) and 24 distinct orientations. The configurational entropy of the complex is dramatically larger, as is the space of possible misassembled states. Kinetic frustration is qualitatively different in 3D: wrong contacts are harder to escape because there are more neighbouring sites to which a particle could be accidentally bonded, and collective motions of the chain are more constrained.

### 6.2 The 3×3×3 Hilbert Path

A single RNA chain of 27 segments traces a Hamiltonian path through the 27 sites of a 3×3×3 cubic grid. The path is chosen to maximise locality (consecutive segments are d=1 apart; as many as possible of the d=1 contacts are between non-consecutive segments that are distant in the chain — these will be the tertiary contacts).

One specific path that works well is a layer-by-layer boustrophedon:

**Layer z=0** (9 sites, indexed 0–8):
```
(0,0,0)→(1,0,0)→(2,0,0)→(2,1,0)→(1,1,0)→(0,1,0)→(0,2,0)→(1,2,0)→(2,2,0)
```

**Layer z=1** (9 sites, indexed 9–17, entered from (2,2,0)):
```
(2,2,1)→(1,2,1)→(0,2,1)→(0,1,1)→(1,1,1)→(2,1,1)→(2,0,1)→(1,0,1)→(0,0,1)
```

**Layer z=2** (9 sites, indexed 18–26, entered from (0,0,1)):
```
(0,0,2)→(1,0,2)→(2,0,2)→(2,1,2)→(1,1,2)→(0,1,2)→(0,2,2)→(1,2,2)→(2,2,2)
```

This path has 26 backbone bonds and produces 27 particles with many non-backbone d=1 contacts between segments in different layers — these are the tertiary contacts.

### 6.3 Contact Classification in 3D

Identify all d=1 non-backbone pairs in the native structure (segments at d=1 but not consecutive in the chain). Classify by sequence separation:

- |i−j| ≤ 6 (within same layer or adjacent layers, close in chain): **local RNA-RNA**, apply γ₁
- |i−j| > 12 (far apart in chain, crossing layers): **global bridging (implicit protein)**, apply γ₂
- 6 < |i−j| ≤ 12: **intermediate**, apply γ₂

### 6.4 Counting Contacts

For the 3×3×3 Hilbert path above, the number of d=1 non-backbone contacts in the native structure can be computed. Each layer (9 sites in a 3×3 grid) has 12 d=1 pairs; along the chain, consecutive segments account for 8 of these. The remaining 4 per layer are non-backbone contacts. Between layers, there are 9 pairs of sites at d=1 in z-direction; minus the 1 inter-layer backbone bond = 8 inter-layer non-backbone contacts. Total non-backbone d=1 contacts: approximately 3×4 + 2×8 = 28 contacts.

For a 27-particle complex with 26 backbone bonds and ~28 non-backbone contacts: this is a complex with a rich energy landscape. At J = 8, the reference energy would be:
```
E_ref = 26 × (−J_bb) + 28 × (−J) = 26×(−20) + 28×(−8) = −520 − 224 = −744 kBT
```
(using J_bb = 20 and J = 8 as recommended). This is a deeply stable complex with many intermediate states.

### 6.5 Patch Assignments in 3D

The patch system extends naturally to 3D. Each particle has 6 body-frame patch slots (±x, ±y, ±z in local frame). The assignment proceeds identically to the 2D case: for each non-backbone d=1 contact (i, j) in the native structure, particle i receives a patch in the direction of j (in i's local frame), and particle j receives a patch toward i.

The 3D orientation state space has 24 elements (the rotation group of the cube): all 24 proper rotations map the 6 cardinal directions to one another. VMMC rotation moves in 3D propose rotations from this set. An efficient implementation stores the orientation as an integer from 0–23 indexing the 24 rotations, with a precomputed table for composition of rotations.

### 6.6 VMMC in 3D

The key change for 3D VMMC:
1. Translation moves: 6 cardinal directions (±x, ±y, ±z) instead of 4
2. Rotation moves: cluster rotation by 90°, 180°, or 270° around any of 3 axes (9 non-identity rotations instead of 3)
3. In-place reorientation: single-particle orientation drawn from 24 possible orientations instead of 4
4. Stokes drag: D ∝ 1/R in 3D; for a sphere of N particles, R ~ N^(1/3), so the drag factor scales as N^(1/3) rather than N^(1/2)

The φ_rot parameter now controls the rate of the 9 rotation types. The SED-motivated value for φ_rot in 3D is:

D_rot / D_trans = k_BT/(8πηr³) / (k_BT/6πηr) = 3/(4r²)

For r ≈ 1 lattice unit: D_rot ≈ 0.75 × D_trans, so φ_rot ≈ 0.43 (same order as the 2D SED-motivated value).

---

## 7. Model C: 3D Multi-Chain — The Full Model

### 7.1 Overview

This is the recommended model for a publication aimed at eLife or PLOS Computational Biology. It combines the 3D cubic lattice with explicit RNA and protein chains, the dual-gradient architecture, and a hierarchical assembly structure that closely mirrors the known ribosome assembly pathway.

### 7.2 Complex Structure

**RNA chain** (1 chain, 27 segments): traces the 3×3×3 Hilbert path as in Model B. This is the rRNA analog.

**Protein chains** (6 dimers, 12 segments): each dimer occupies 2 adjacent lattice sites on the *surface* of the 3×3×3 cube. Six protein binding sites are placed at six distinct locations on the RNA surface, arranged to require that specific RNA tertiary contacts have formed before the protein can dock.

**Total complex**: 27 + 12 = 39 particles.

The 3×3×3 cube has a surface of 26 sites, of which 6 face outward in distinct cardinal directions. Proteins are placed on these 6 outer faces, each protein dimer extending one step outward from the cube surface. The outer particle of each protein dimer has no patches (non-sticky outward face — preventing inter-complex protein-protein contacts).

### 7.3 Assembly Hierarchy (Three Tiers)

**Tier 1 — RNA folding** (local RNA-RNA contacts, γ₁, favoured at low pH/interior):
RNA chain forms its secondary structure (local contacts, |i−j| ≤ 6). This is required before any protein can bind, because protein binding sites are at RNA sites that are only accessible when the local RNA structure is correct. These contacts form in the interior of the condensate.

**Tier 2 — Primary protein binding** (RNA-protein contacts for 3 "primary" proteins, γ₂, favoured at periphery):
Three of the six proteins have binding sites that require only Tier 1 to be complete. As the particle moves toward the condensate boundary, γ₂ increases and these proteins can dock. Their binding stabilises the RNA local structure.

**Tier 3 — Cross-domain contacts + secondary protein binding** (global RNA-RNA contacts + remaining 3 proteins, γ₂, also favoured at periphery):
With the primary proteins bound, the RNA chain's global cross-domain contacts (|i−j| > 12) are stabilised. The remaining three proteins have binding sites that require both Tier 1 and Tier 2 to be complete — they are the "tertiary binders." They dock last, completing the complex.

This three-tier assembly pathway is directly analogous to the Nomura assembly map of the bacterial 30S ribosomal subunit (Nomura 1973; Held et al. 1974) and is emergent from the physics of the model rather than hard-coded.

### 7.4 Why the Three Tiers Emerge Naturally

- **Tier 1 occurs first** because γ₁ is large in the interior (local RNA-RNA contacts favoured) while γ₂ is small (protein binding suppressed). RNA folding happens without protein interference.
- **Tier 2 occurs second** because as γ₂ increases with radial position, the primary proteins (lowest binding cost, closest to RNA surface, fewest misassembly options) dock first.
- **Tier 3 occurs last** because tertiary protein binding sites are sterically inaccessible (blocked by RNA cross-domain contacts that haven't yet formed) or thermodynamically unfavourable (require additional correct neighbours) until Tiers 1 and 2 are complete.

No explicit time ordering or cooperative kinetics rule is needed — the hierarchy emerges from the spatial gradient and the geometry of the complex.

---

## 8. Comparison of Models

### 8.1 Kinetic Frustration

| Model | Frustrated intermediates | Dominant trap type |
|---|---|---|
| Current (4-polymer) | Partial n=1,2,3 polymer attachments | Any polymer in wrong order |
| S (single chain) | Local folding traps; wrong-direction threading | Local structure forming before global |
| A (2D, RNA+protein) | RNA trapped in wrong fold by early protein binding | Protein locks misfolded RNA |
| B (3D, single chain) | 3D knots, layer inversion, wrong-layer contacts | Wrong inter-layer contacts |
| C (3D, RNA+protein) | All of B plus protein-mediated RNA traps | Primary protein blocking tertiary binding |

Model C has the richest frustration landscape: wrong-order protein binding (a protein from Tier 3 binds before Tier 2 is complete) locks the RNA in a non-native inter-domain topology that cannot be corrected without first removing the protein. This is directly analogous to the known ribosome assembly failure mode: early binding of late-stage r-proteins traps misfolded rRNA intermediates.

### 8.2 Size Hierarchy

In the current model, all 4 polymers are equal in length (4 segments). In the real ribosome, the rRNA is ~100-fold longer by mass than any individual r-protein. In Model C:

- RNA chain: 27 segments (dominant component, the scaffold)
- Protein dimers: 2 segments each (6 proteins, total 12 segments)
- Mass ratio of RNA to one protein: 27:2 = 13.5 (more realistic than 1:1)
- Stokes drag: RNA chain (when folded, size ~3 lattice units diameter) has D ∝ 1/3; protein dimer (size ~1–2 units) has D ∝ 1/1–1/2. This size difference produces the sorting effect: free RNA chains diffuse 3-fold slower than free protein dimers, mimicking the nucleolar dynamics where rRNA is gel-like (slow) while proteins are mobile (FRAP t₁/₂ ~ seconds).

### 8.3 Parameter Comparison

| Parameter | Current model | Model A | Model C |
|---|---|---|---|
| J (weak contacts) | 1 (J = J₁ = J₂) | 2 (J₁, J₂) | 2 (J₁, J₂) |
| J_bb | 1 (= 1000) | 2 (RNA vs. protein) | 2 (RNA vs. protein) |
| γ₀ | 1 | 1 | 1 |
| Gradient type | 1 | 1 (dual-channel) | 1 (dual-channel) |
| Chain lengths | 1 (all equal = 4) | 2 (RNA = 16, protein = 2) | 2 (RNA = 27, protein = 2) |
| **Total** | **5** | **7** | **7** |

The addition of the dual-gradient requires only one extra free parameter (J₁ separately from J₂). This is very lean for the biological richness gained.

---

## 9. The pH Mapping in Detail

With the dual-gradient, the coupling of each contact type to the physical pH gradient (King et al. 2024) is explicit:

**RNA-RNA local contacts** (secondary structure):
```
ΔG_RNA-RNA(pH) ≈ −J₁ × [1 − (1−γ₀) × r/R_c]
```
At r=0: ΔG = −J₁ (full coupling). At r=R_c: ΔG = −J₁ × γ₀ (reduced). 

The pH corresponding to r/R_c: pH(r) = 6.5 + 0.7 × r/R_c. The RNA-RNA coupling at pH 6.5 is J₁ and at pH 7.2 is J₁ × γ₀. A γ₀ = 0.7 for RNA-RNA contacts means they are 70% as strong at physiological pH as at nucleolar interior pH — physically motivated by the ~0.43 kcal/mol per pH unit stabilisation of RNA duplexes measured at high ionic strength.

**RNA-protein contacts** (r-protein binding):
```
ΔG_RNA-prot(pH) ≈ −J₂ × [γ₀ + (1−γ₀) × r/R_c]
```
At r=0: ΔG = −J₂ × γ₀ (suppressed). At r=R_c: ΔG = −J₂ (full coupling).

The RNA-protein coupling at pH 6.5 is J₂ × γ₀ and at pH 7.2 is J₂. The direct pH effect on r-protein:rRNA binding thermodynamics through histidine charge-state changes is small (~2–5% of total binding energy; see gradients.md §1.2) and cannot alone justify γ₀ = 0.3–0.4. The route to this range is through assembly factor impairment: DDX21 activity is sub-optimal at pH 6.5, reducing the effective rate at which correct RNA-protein contacts form and kinetic traps are resolved. The 15-fold change in S4-rRNA Kd over 2–8 mM Mg²⁺ (PMC4604426) calibrates the energy scale (supporting J₂ ≈ 5–10 kBT) but describes Mg²⁺ dependence, not the pH gradient. γ₀ = 0.3–0.4 is treated as a phenomenological parameter consistent with the interior being substantially less assembly-competent than the periphery, with DDX21 pH sensitivity as the clearest direct evidence.

**Important**: For a given physical system, γ₀ should be the same parameter for both channels (it represents the minimum coupling at the condensate centre). The difference between the two channels is just the *direction* of the gradient:

- γ₁(r) = 1 − (1−γ₀) × r/R_c
- γ₂(r) = γ₀ + (1−γ₀) × r/R_c

Both equal (1+γ₀)/2 at r = R_c/2. Both are specified by one parameter γ₀. The physical interpretation of γ₀ is the ratio of coupling strength at the condensate periphery to the condensate centre for RNA-protein contacts (or, equivalently, the ratio of RNA-RNA coupling at the periphery to the centre). With γ₀ = 0.4, RNA-protein contacts are 40% as strong at the centre as at the edge; RNA-RNA contacts are 100% at the centre and 40% at the edge.

---

## 10. Space-Filling Curves: Explicit Coordinates

### 10.1 2D: 5×5 Peano-Hilbert Path (25 segments, alternative to 4×4)

For a longer single-chain 2D model:
```
y=4: 24 23 22 21 20
y=3: 15 16 17 18 19
y=2: 14 13 12 11 10
y=1:  5  6  7  8  9
y=0:  4  3  2  1  0
```
Path: (0,0)→(1,0)→(2,0)→(3,0)→(4,0)→(4,1)→(3,1)→(2,1)→(1,1)→(0,1)→(0,2)→(1,2)→...

This boustrophedon (back-and-forth) 5×5 path has 24 backbone bonds and ~30 non-backbone d=1 native contacts. With contact classification by sequence separation, it gives a well-defined local/global structure.

**Domain structure** for 5-segment "domains":
- Domain 1: segments 0–4 (bottom row)
- Domain 2: segments 5–9 (second row)
- Domain 3: segments 10–14 (middle row)
- Domain 4: segments 15–19 (fourth row)
- Domain 5: segments 20–24 (top row)

Inter-domain contacts (between non-adjacent domains, |i−j| > 10): these bridge the five domains together and are the long-range contacts that proteins would mediate implicitly.

### 10.2 3D: 3×3×3 Layer Path (27 segments)

Explicitly:

**Layer z=0** (segments 0–8):
```
(0,0,0) (1,0,0) (2,0,0)
(2,1,0) (1,1,0) (0,1,0)
(0,2,0) (1,2,0) (2,2,0)
```
Segment indices: 0→1→2→3→4→5→6→7→8 (boustrophedon in x, increasing y)

**Layer z=1** (segments 9–17, entered via (2,2,0)→(2,2,1)):
```
(2,2,1) (1,2,1) (0,2,1)
(0,1,1) (1,1,1) (2,1,1)
(2,0,1) (1,0,1) (0,0,1)
```
Segment indices: 9→10→11→12→13→14→15→16→17 (reverse boustrophedon)

**Layer z=2** (segments 18–26, entered via (0,0,1)→(0,0,2)):
```
(0,0,2) (1,0,2) (2,0,2)
(2,1,2) (1,1,2) (0,1,2)
(0,2,2) (1,2,2) (2,2,2)
```
Segment indices: 18→19→20→21→22→23→24→25→26

**Backbone bonds** (26 bonds):
Each consecutive pair is exactly d=1 apart. The layer-transition bonds are (8→9) at (2,2,0)→(2,2,1) and (17→18) at (0,0,1)→(0,0,2).

**Non-backbone d=1 contacts** in this native structure: there are contacts between z-adjacent layers (sites at same x,y but z differing by 1) that are non-consecutive in the chain. These are the tertiary contacts.

Count of inter-layer non-backbone contacts:
- Between layers 0 and 1 (z=0 to z=1): sites (x,y,0) and (x,y,1) for all x,y ∈ {0,1,2}: 9 pairs. Minus the backbone bond (2,2,0)→(2,2,1): 8 non-backbone inter-layer contacts between layers 0 and 1.
- Between layers 1 and 2: similarly 8 non-backbone contacts.

Within-layer non-backbone contacts: within each 3×3 boustrophedon layer, the 12 cardinal d=1 pairs minus 8 backbone bonds = 4 non-backbone contacts per layer.

**Total non-backbone d=1 contacts**: 3×4 (within layers) + 2×8 (between layers) = 12 + 16 = **28 native patch contacts**.

Reference energy at γ = 1:
```
E_ref = 26 × (−J_bb) + 28 × (−J)
```
With J_bb = 40, J = 8: E_ref = −1040 − 224 = −1264 kBT.

### 10.3 Contact Table for 3D Model (sample)

Non-backbone d=1 contacts by sequence separation for the 3×3×3 layer path:

| Pair (i,j) | |i−j| | Contact type | γ |
|---|---|---|---|
| (0,5) | 5 | within layer 0 | local RNA-RNA | γ₁ |
| (2,9) | 7 | layer 0→1 (corner) | local RNA-RNA | γ₁ |
| (0,18) | 18 | layer 0→2 | global (protein-bridged) | γ₂ |
| (6,24) | 18 | layer 0→2 | global (protein-bridged) | γ₂ |
| (3,12) | 9 | layer 0→1 | intermediate | γ₂ |

The cutoff between γ₁ and γ₂ is set at |i−j| = 9 (one full layer separation). Local contacts (same layer or adjacent-layer backbone neighbours, |i−j| ≤ 9) use γ₁; global contacts (cross-layer bridging, |i−j| > 9) use γ₂.

---

## 11. Implementation Roadmap

### 11.1 What Stays the Same

- VMMC algorithm (Whitelam–Geissler, symmetrised version)
- Kern–Frenkel patch model (patch slot convention)
- Condensate circular geometry and open-boundary recycling
- exitQuality metric and stats output

### 11.2 Model S (Drop-in, 2D) — Changes Required

1. **In run_condensate.cpp**: Change TARGET_X/Y to the single-chain path ordering; remove polyType() function (all same "polymer"); add contactType() function that classifies by |i−j|.
2. **In buildCouplingMatrices()**: Fill weakD1_intra[i][j] = J₁ for local contacts, weakD1_inter[i][j] = J₂ for global contacts.
3. **In CondensateModel.cpp**: `computePairEnergy()` applies `couplingFactor_1()` to intra-entries and `couplingFactor_2()` to inter-entries. Add `hasGradient_inverted` flag for γ₁.
4. **In gamma_r()**: Add a flag to return the inverted gradient.

Estimated implementation time: 1–2 days.

### 11.3 Model A (2D Multi-chain) — Additional Changes

5. Add explicit protein particles with their own polymer type, backbone bond strength, and patch assignments.
6. Modify `polyType()` to return RNA vs. protein.
7. Adjust `referenceComplexEnergy()` to account for RNA-protein contacts.

Estimated additional implementation time: 2–3 days.

### 11.4 Model B/C (3D) — Major Changes

8. **Box.h/cpp**: Extend to 3D box with 6-direction periodic boundary conditions.
9. **Particle.h**: Orientation now stored as integer 0–23 (the 24 cubic rotations) or as a rotation matrix; body-frame patch slot computation in 3D.
10. **VMMC.h/cpp**: Rotation proposals now choose from 9 non-identity rotations (90°, 180°, 270° around each of 3 axes); translation proposals choose from 6 cardinal directions.
11. **StickySquare.cpp**: Patch slot computation for 6-faced particle: `getPatchSlot()` returns one of 6 values (±x, ±y, ±z in local frame).
12. **CondensateModel.cpp**: Condensate becomes a sphere (3D), with `gamma_r(r) = gamma_r(sqrt(dx²+dy²+dz²))`. Exit condition: r > R_c in 3D.
13. **Stokes drag**: 3D Stokes drag: D ∝ 1/R ∝ N^(-1/3) for a cluster of N particles.

Estimated implementation time: 2–4 weeks (non-trivial; 3D VMMC requires careful implementation and testing of the orientation algebra).

---

## 12. Predicted New Results

### 12.1 The Dual-Gradient Outperforms All Single-Gradient Controls

With the dual-gradient (γ₁ inverted, γ₂ normal), four conditions can be compared at fixed J₁, J₂, and γ₀:

| Condition | γ₁ | γ₂ | Predicted exitQuality |
|---|---|---|---|
| Dual gradient | inverted | normal | Best |
| γ₂ only (no inversion) | uniform 1 | normal | Intermediate |
| γ₁ only (RNA-RNA gradient only) | inverted | uniform γ₀ | Poor |
| No gradient | uniform 1 | uniform γ₀ | Worst |
| Inverted (both wrong) | normal | inverted | Worst |

This is a clean 5-condition test of the core hypothesis. The prediction that the *dual* gradient outperforms either single gradient is strong and falsifiable. It directly demonstrates that the physical distinction between RNA-RNA and RNA-protein contacts — and their opposite pH sensitivities — is necessary for optimal assembly.

### 12.2 Hierarchical Assembly Timeline

Measure the mean radial position at which each contact tier first forms (on successful assemblies). The prediction is:
- Tier 1 (RNA-RNA local): forms at r/R_c ≈ 0.2–0.4 (interior, where γ₁ is high)
- Tier 2 (primary proteins): r/R_c ≈ 0.5–0.7 (mid-condensate, as γ₂ rises)
- Tier 3 (tertiary proteins + global RNA contacts): r/R_c ≈ 0.7–0.9 (periphery)

This would be a beautiful figure: a stacked radial probability distribution of when each tier's contacts form, showing the spatial separation of the assembly stages. It directly corresponds to the Quinodoz et al. (2025) sequencing data showing temporal separation of rRNA processing events.

### 12.3 Phase Diagram in (γ₀, J₂/J₁) Space

The ratio J₂/J₁ controls the relative importance of RNA-RNA versus RNA-protein contacts. When J₂/J₁ >> 1, protein binding dominates; when J₂/J₁ << 1, RNA folding dominates. The assembly phase diagram in (γ₀, J₂/J₁) space would show:
- A region of successful assembly (the gradient correctly sequences RNA folding then protein binding)
- A region of "protein-first trap" (J₂ too large, proteins bind before RNA folds)
- A region of "RNA trap" (J₁ too large, RNA misfolding traps prevent protein binding)
- A region of "no assembly" (both J₁ and J₂ too small)

The physical prediction is that the ribosome operates in the central "successful assembly" region, with the Mg²⁺/pH conditions setting the effective (J₂/J₁, γ₀) values.

---

## 13. Summary Recommendation

**For the first paper** (achievable in weeks to months with existing 2D code):

Implement **Model S** or **Model A** by modifying the existing run_condensate binary. The dual-gradient is the key new physics — it can be demonstrated in 2D and the biological interpretation is already clear from the pH gradient data. Run the 5-condition comparison (§12.1), produce the assembly phase diagram (§12.3), and include the radial contact-formation analysis (§12.2). This is a complete paper.

**For the full publication** (3D, larger complex, full hierarchy):

Implement **Model C** (3D, RNA + explicit proteins, 39 particles/complex). This is a more significant code project (2–4 weeks for 3D VMMC) but produces qualitatively new physics: genuine 3D frustration, richer intermediate states, and a direct correspondence with the three-tier ribosome assembly hierarchy. A paper based on Model C would be more likely to be accepted at eLife.

The single most important new physics in any of these models is the **inverted gradient for RNA-RNA contacts**. Even implementing this in the current 4-polymer model (treating same-polymer contacts as RNA-RNA and cross-polymer contacts as RNA-protein) would be a meaningful improvement and a direct test of the pH hypothesis, achievable in a day.
