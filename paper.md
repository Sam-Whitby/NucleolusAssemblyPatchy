# Towards a Paper: Spatial Chemical Gradients Drive Emergent Annealing in a Model of Ribosome Assembly Within the Nucleolus

*A critical research roadmap and assessment — S. Whitby, 2026*

---

## 1. Executive Summary

The central idea is this: the nucleolus is not merely a container for ribosome assembly; its radial chemical gradients (of Mg²⁺, pH, and ionic strength) function as a spatial annealing schedule that is traversed automatically by the outward advection of assembling particles. A ribosomal subunit particle born near the transcription centre experiences denaturing conditions that prevent premature lock-in of misassembled contacts; as it drifts outward, coupling strengthens gradually and the correct native complex can nucleate. This is thermodynamically passive annealing, requiring no chaperone energy expenditure, emerging purely from the geometry of the condensate.

This is a genuinely interesting, non-obvious idea. The minimal model in this codebase — a 2D circular condensate with a Kern–Frenkel patchy Go̊-model, open exit boundary, and a radial gradient of the coupling strength γ — is the right kind of model to demonstrate it. If the simulations are done correctly and thoroughly, and the biological connection is made quantitatively, this constitutes a publishable claim.

**The current state of the work, however, is insufficient for submission.** The comparison you describe between "working" and "non-working" parameter sets conflates three variables simultaneously; the statistical basis is single-seed runs; and the model's connection to real nucleolar chemistry is qualitative at best. This document sets out what you need to do.

---

## 2. The Central Claim and What Makes It Novel

### 2.1 The Claim

A spatially organised radial gradient of interaction strength inside a circular condensate, combined with outward advective drift, enables self-assembly of kinetically frustrated multi-polymer complexes into their correct target structure with a significant yield. The same assembly at uniform high coupling (no gradient) fails, producing trapped misassembled intermediates. The same assembly at uniform low coupling (strong denaturing) also fails, because bonds cannot form. The gradient, traversed spatially, automatically implements a near-optimal annealing schedule.

### 2.2 What Is Novel

This claim connects four things that have not previously been synthesised:

1. **Patchy particle assembly** — well-studied in colloidal systems (Romano & Sciortino, Wilber et al.), but not yet applied to ribosomal-geometry condensates.

2. **Spatial gradients as temporal annealing schedules** — Trubiano & Whitelam (2022, *PNAS*) showed theoretically that an optimal control schedule can maximise yield for a colloidal chain; the nucleolus implements this schedule *spatially* rather than requiring external temporal control.

3. **The nucleolus as an annealing machine** — the biological literature on ribosome reconstitution (Nierhaus 1990; Guth-Metzler et al. 2023) notes that a Mg²⁺ ramp is required for efficient assembly *in vitro*, and that assembly *in vivo* involves radial transit through distinct chemical phases, but no simulation has modelled this as an annealing phenomenon.

4. **Emergent advection from Stokes drag** — the size-dependent drag biases diffusion of larger clusters toward the periphery, meaning that the gradient and the advective biasing are coupled: assembled complexes exit preferentially, while free monomers remain in the interior and can re-fold.

None of these has been combined before in a single model. The novelty is real.

### 2.3 A Critical Warning About Your Current Comparison

The two command lines you show are **not a valid comparison**:

```
Case A: --gamma0 0.4  --J  8   --gradient   (gradient + annealing)
Case B: --gamma0 0.9  --J 80   --gradient   (gradient + strong coupling, no annealing?)
```

You have changed **both J and gamma0 simultaneously**. Case B has J=80 (ten times stronger binding) and gamma0=0.9 (almost no gradient — the inner boundary condition is already near physiological). This is not "the same advective current but no annealing": it is a completely different physical regime with far stronger binding and almost no gradient extent. A frustrated-assembly failure in Case B could simply be because J=80 creates deep misassembly traps, not because annealing is absent.

To demonstrate the core claim, you need a controlled comparison:
- Same J (e.g. J=8)
- Same R_c, same gamma0
- Case A: `--gradient` enabled → γ(r) = γ₀ + (1−γ₀)r/R_c
- Case B: no `--gradient` flag → γ(r) = 1 everywhere
- Case C: no `--gradient`, but with a box simulation using `--anneal` ramp → pure temporal annealing

This three-way comparison at fixed J would be the most important figure in any paper.

---

## 3. Critical Assessment of the Model

### 3.1 Strengths

- **Detailed balance is exact** — VMMC guarantees this by construction; it is the single most important virtue of this approach over coarse-grained MD alternatives.
- **Patchy Kern–Frenkel model prevents spurious aggregation** — this is the key improvement over the original NucleolusAssembly code and is physically well-motivated.
- **Open boundary with recycling** — physically realistic representation of the nucleolar condensate boundary.
- **Three-phase initialisation** (equil→denature→run) — a sensible way to ensure the system starts from a genuinely disordered state.
- **exitQuality metric** — using the fraction of exiting mass that is in perfect complexes is a good observable because it disentangles assembly quality from simple kinetic effects (how many particles exit at all).

### 3.2 Fundamental Model Weaknesses

#### 3.2.1 Two Dimensions

The ribosome is a three-dimensional object. On a 2D square lattice, particles have only four possible orientations. This severely restricts the configurational entropy of misassembled states — the frustration that the model generates is qualitatively different from 3D frustration. The 2D lattice underestimates the number of wrong-contacts states that are accessible, which means the real assembly problem is harder and the annealing effect more important.

**Remedy**: Extend to a 3D cubic lattice, where each particle has 6 faces and a richer set of misassembly intermediates is accessible. This would also allow a sphere (rather than circle) condensate geometry, more faithful to the real nucleolus.

#### 3.2.2 Lattice Geometry Limits Orientational Entropy

With only 4 discrete orientations, the entropy cost of correctly orienting a particle is ln(4) ≈ 1.4 k_BT, which is small. In a continuous-angle model the orientational entropy is larger, and the effect of the gradient on exploring orientational space would be more pronounced.

**Remedy**: Use a continuous-space VMMC (the original Whitelam–Geissler formulation, off-lattice) with Kern–Frenkel patches defined by angular width δ. The code exists at vmmc.xyz and the patch model has already been validated by Rovigatti et al. (2018).

#### 3.2.3 Extreme Backbone Stiffness (J_bb = 1000)

The backbone bonds at strength 1000 k_BT are essentially indestructible and mean that each polymer chain is completely rigid and never breaks apart. Real RNA and protein chains have finite flexibility and their own folding dynamics. The high J_bb prevents the simulation from exploring states where chain connectivity competes with wrong-contact formation.

**Remedy**: Keep J_bb high enough to maintain chain connectivity (J_bb ≈ 10–20 times J is sufficient to ensure chains stay together thermodynamically) but not so high that they never transiently stretch. At J_bb = 20 for J = 8, chains would rarely break but would explore more realistic chain conformations.

#### 3.2.4 Toy Complex Size and Geometry

A 4×4 = 16 particle complex with 4 identical-length polymers arranged on a Moore curve is elegant as a model system, but is far from any real ribosomal sub-unit. The large ribosomal subunit (60S in eukaryotes) contains 47 proteins and three rRNA molecules (~4500 nt total). The small subunit (40S) contains 33 proteins and one rRNA molecule (~1900 nt). The key structural features are:
- Large asymmetry in component sizes (rRNA is much larger than any single r-protein)
- Hierarchical assembly: some proteins bind early (primary binders), others require pre-assembled intermediates (secondary and tertiary binders)
- Specific binding geometries that are not well captured by the Moore curve

**Remedy**: Use an n=3 Moore curve (36 particles) or switch to a hierarchical assembly model where a large "RNA core" particle acts as the scaffold and smaller "protein" particles bind to it sequentially. This would be more biology-faithful and allow you to demonstrate that the gradient assists the sequential, hierarchical nature of assembly.

#### 3.2.5 No Explicit Size Hierarchy

In the current model, all 16 particles are identical in size. The Stokes drag (`--stokes`) depends on the cluster size, which is realistic for assembled complexes, but the individual particles are all point-like on the lattice. In the real nucleolus, ribosomal proteins are 10–100 times smaller by mass than the rRNA molecule they bind to. Larger particles diffuse more slowly, spend more time in the interior of the condensate, and encounter the full gradient profile. Smaller particles might equilibrate faster.

**Remedy**: Assign different effective diffusion coefficients to each polymer type (or segment), reflecting size differences. Even in the 2D lattice model this could be implemented by modulating the Stokes factor per particle type.

#### 3.2.6 Static Gradient

The gradient in the model is fixed: γ(r) = γ₀ + (1−γ₀)r/R_c. The real nucleolus has dynamic, fluctuating gradients because the chemical concentrations (Mg²⁺, pH) are maintained by active enzymatic processes that fluctuate in time. Temporal fluctuations in the gradient could actually enhance assembly by providing additional annealing noise.

**Remedy**: Add a time-dependent fluctuation to γ(r, t) = γ(r) + η(t) where η(t) is a Gaussian noise with amplitude σ and correlation time τ. Sweep over σ and τ to find the regime where fluctuations help vs. hinder assembly.

#### 3.2.7 Gradient Shape Is Linear — Is It Optimal?

The choice γ(r) = γ₀ + (1−γ₀)r/R_c is linear. There is no theoretical reason why the optimal gradient should be linear. The Sivak–Crooks framework for optimal control of thermodynamic transitions (and the Trubiano–Whitelam 2022 paper) suggests that the optimal schedule minimises dissipation, which in general gives a non-linear schedule shaped by the curvature of the free energy landscape.

**Remedy**: After establishing the baseline result with a linear gradient, sweep over different functional forms: concave (γ ~ (r/R_c)^α with α < 1), convex (α > 1), sigmoidal (γ ~ tanh(r−r₀)). This would be an impressive result: showing that the nucleolus has approximately the optimal gradient profile.

---

## 4. Simulation Roadmap

### Phase 1: Establish the Core Claim (Essential)

**Objective**: Demonstrate, with proper statistics, that a radial gradient enables assembly that a uniform coupling cannot.

**Simulations needed**:

```
# Fixed parameters for all:
--copies 8  --t-equil 500  --t-denat 1000
--steps 200000  --snapshots 200
--radius 60  --J 8
--stokes  --coupling midpoint
--phi-rot 0.2  --phi-reorient 0.2

# Condition A: Gradient (your working case)
--gamma0 0.4  --gradient  [seeds 1–20]

# Condition B: Full coupling, no gradient (uniform γ=1)
[no --gradient]  [seeds 1–20]

# Condition C: Low coupling everywhere (uniform γ=0.4)
# Cannot directly set this in condensate; use box model equivalent
# Or: --gamma0 0.4 but no --gradient (so gradient overrides to g=1 with no gradient flag)
# Actually: run_box --anneal --t-denat 1000 --steps 200000 (temporal annealing)
```

**Metrics to report**:
- Mean and standard deviation of `exitQuality` over 20 seeds
- Mean `exitedPerfect` / (total possible perfect exits) — assembly efficiency
- Fraction of time spent in each assembly state (free monomers, partial complexes, full complexes)
- Cumulative `exitedPerfect` vs time curve — assembly kinetics

**Expected result**: Condition A significantly outperforms B and C on all metrics.

### Phase 2: Parameter Sensitivity (Important)

**γ₀ scan**: At fixed J=8, vary γ₀ ∈ {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9} — 10 conditions × 20 seeds. Plot exitQuality vs γ₀. This produces the "optimal γ₀" figure and shows that too much denaturing (γ₀→0) is as bad as too little (γ₀→1), producing a peak.

**J scan**: At fixed γ₀=0.4, vary J ∈ {4, 6, 8, 10, 12, 16}. This shows the sensitivity to the absolute binding energy scale.

**Radius scan**: Vary R_c ∈ {30, 45, 60, 90, 120} — does a larger condensate (longer time in gradient) improve yield? This connects to the known 10 nm/s drift speed measured in pulse-chase experiments.

**Gradient shape**: Implement non-linear γ(r) profiles (see §3.2.7). Run with α=0.5, 1.0, 2.0, and sigmoid.

### Phase 3: Mechanistic Analysis (Important for Understanding)

**Trajectory analysis**: For successful assemblies (exitQuality > 0), trace the assembly pathway in position-space and energy-space. At what radial position does each native contact first form? This would produce a beautiful figure showing assembly progressing as particles move outward.

**Frustration measurement**: Track the number of non-native patch contacts (wrong-type particles with patches touching) as a function of radial position. If the gradient works correctly, wrong contacts should be enriched near the interior and native contacts at the periphery.

**Advection analysis**: Track mean radial position of: (i) free monomers, (ii) partial complexes of size 2–15, (iii) full 16-particle complexes. The Stokes drag should produce a sorting effect where larger clusters are biased outward, which is a key mechanistic element.

**Comparison with uniform temporal annealing**: Run `run_box --anneal` with the same J and copying the best γ₀ schedule as a linear temporal ramp. This comparison establishes that spatial annealing is not merely equivalent to temporal annealing.

### Phase 4: Physical Parameter Connection (Critical for Publication)

Quantitatively connect γ to real chemical parameters — see §6 of this document.

### Phase 5: Model Extensions (Desirable but Not Essential for First Paper)

- 3D extension
- Larger complex
- Gradient fluctuations

---

## 5. Relevant Literature

### 5.1 Ribosome Biology and the Nucleolus

**Warner (1999)** [*Trends Biochem. Sci.*] — The economics of ribosome biosynthesis in yeast. Establishes the scale of ribosome production (2000 per minute in yeast) and its metabolic cost, motivating efficient assembly mechanisms.

**Woolford & Baserga (2013)** [*Genetics*] — Ribosome biogenesis in the yeast *Saccharomyces cerevisiae*. The definitive review of the 200+ assembly factors involved. The assembly is hierarchical: ~20 "primary" ribosomal proteins bind first to the pre-rRNA, defining the scaffold for later binders. This hierarchy is exactly the kind of sequential process that an annealing gradient would facilitate.

**Nierhaus (1990)** [*Biochemistry*] — Reconstitution of ribosomes from isolated components demonstrates that *in vitro* assembly requires specific protocols: a two-temperature step (40°C then 50°C for the large subunit) and a Mg²⁺ ramp. This is the key biological fact motivating the entire thesis.

**Guth-Metzler et al. (2023)** [*PNAS*, "Goldilocks Mg²⁺"] — Quantifies how the Mg²⁺ concentration must be within a narrow "Goldilocks" range for ribosome assembly: too little and RNA cannot fold, too much and aggregation occurs. This defines the physical basis for the γ₀ parameter: γ₀ should be above the minimum threshold for RNA folding.

**Riback et al. (2023)** [*Mol. Cell*] — Pulse-chase and fluorescence tracking shows assembling ribosomal particles move radially outward through the nucleolus at ~10 nm/s. This is the advective drift that the Stokes mechanism models. The paper also characterises the viscoelastic properties of each nucleolar sub-compartment.

**Lafontaine et al. (2021)** [*Nature Rev. Cell Biol.*] — The most recent comprehensive review of nucleolus biology. Establishes the three-phase structure (FC, DFC, GC) with distinct chemical compositions, confirming that the radial position of a ribosomal particle corresponds to its maturation state.

**Detres et al. (2025)** [*preprint or recent*] — pH variation in the nucleolus affects RNA-protein interactions. This is direct evidence for the kind of chemical gradient the model proposes.

### 5.2 Chemical Gradients in Condensates

**King et al. (2024)** [*Cell*] — Demonstrates that chemical gradients can exist within single-phase condensates and across phase boundaries. Establishes the physical feasibility of intra-condensate gradients.

**Kilgore et al. (2024)** and **Ausserwoger et al. (2024)** — Recent experimental demonstrations of composition gradients within condensates; directly relevant as evidence that the gradient in the model is realistic.

**Yewdall et al. (2022)** — ATP-Mg²⁺ coupling within condensates suggests a mechanism by which the nucleolus could maintain a Mg²⁺ gradient through enzymatic activity near the transcription centres.

### 5.3 Patchy Particle Assembly

**Kern & Frenkel (2003)** [*J. Chem. Phys.*] — The original patchy particle paper. Defines the Kern–Frenkel potential used in this model. Essential to cite.

**Romano & Sciortino (2011)** [*Nature Materials*] — Self-assembly of patchy particles into designer structures. Demonstrates the power of patchy particle models for controlled self-assembly.

**Whitelam et al. (2009)** [*Biophys. J.*] — Patchy particle models for protein and RNA assembly show that directionality of interactions dramatically improves assembly yields.

**Rovigatti, Russo & Romano (2018)** [*Eur. Phys. J. E*] — "How to simulate patchy particles." Provides the technical details for implementing patchy particle models; the body-frame patch convention used in this code follows their prescription.

**Wilber et al. (2007)** [*J. Chem. Phys.*] — Self-assembly of patchy particles shows that assembly yield has a non-monotonic dependence on temperature (coupling strength), with a maximum at intermediate temperature — this is the key qualitative result that motivates the annealing approach.

**Jacobs et al. (2015)** [*Soft Matter*] — Rational design of patchy particle assembly: shows that for complex target structures, patches must be specifically placed to avoid kinetic traps. The Moore-curve complex design follows this principle.

**Sartori & Leibler (2020)** — Self-assembly under dissipation. Discusses how active processes can improve yield beyond passive self-assembly. Your work is the passive limit of this problem.

### 5.4 Annealing and Optimal Control of Assembly

**Kirkpatrick, Gelatt & Vecchi (1983)** [*Science*] — Simulated annealing. The conceptual ancestor of the whole approach.

**Hajek (1988)** — Cooling schedules for optimal annealing: the logarithmic schedule guarantees convergence to the ground state given infinite time. Relevant for understanding why slow gradients should work better.

**Trubiano & Whitelam (2022)** [*PNAS*] — Optimal bond-strength schedule for assembly of a colloidal chain. The most directly relevant prior work. Shows that optimal control theory can find the annealing schedule that maximises yield in finite time. **Your paper should explicitly position itself relative to this work**: you are showing that the nucleolus implements (approximately) such an optimal schedule *spatially*, without external control.

**Barzegar et al. (2024)** — Optimal annealing schedules from thermodynamic geometry. Extends the Sivak–Crooks framework to assembly. Discusses how the curvature of the free energy landscape determines the optimal γ(t) profile — relevant for §3.2.7.

**Sivak & Crooks (2012)** [*PRL*] — Thermodynamic geometry and optimal protocols. Provides the formal framework for "minimum dissipation" schedules.

### 5.5 VMMC and Lattice Models

**Whitelam & Geissler (2007)** [*J. Chem. Phys.*] — Original VMMC paper. Demonstrates that VMMC avoids unphysical kinetics in cluster self-assembly systems. Essential citation.

**Whitelam & Geissler (2009)** [*Phys. Rev. E*] — The symmetrised VMMC algorithm used in this code. Proves detailed balance for the cluster translation moves.

**Holmes-Cerfon (2016)** — Enumerating rigidity for colloidal clusters. Relevant to understanding which cluster geometries are kinetically accessible.

**Takada (2019)** [*Biophys. Rev.*] — Go̊ model review. Places the interaction scheme in the context of protein folding models.

### 5.6 RNA Folding Kinetics (Relevant for Biology)

**Pan & Sosnick (1997)** — RNA folding and the folding kinetics of the hammerhead ribozyme demonstrate that RNA misfolds readily under physiological conditions.

**Treiber & Williamson (1999)** [*Curr. Opin. Struct. Biol.*] — Exposing the kinetic traps in RNA folding. Shows that RNA tertiary folding is kinetically frustrated for the same reasons as protein complex assembly.

**Adilakshmi et al. (2005)** — Protein-independent folding pathway for the ribosomal RNA: even without proteins, the rRNA has kinetically frustrated folding, reinforcing the need for an annealing-like process.

---

## 6. Physical Chemistry: Connecting γ to Real Chemical Gradients

### 6.1 The γ Parameter in the Model

γ(r) scales all weak patch interactions: E_patch = −γ(r) × J. When γ = 0, no patch contacts form; when γ = 1, full physiological binding. The gradient γ(r) = γ₀ + (1−γ₀)r/R_c represents a smooth transition from a partially denaturing core (γ₀) to fully physiological conditions at the boundary.

### 6.2 Magnesium Ion Concentration

Mg²⁺ is critical for ribosome assembly: it coordinates RNA phosphate groups, neutralises backbone repulsion, and enables RNA tertiary structure formation. The Nierhaus reconstitution protocols require 4–20 mM free [Mg²⁺], while cytoplasmic free [Mg²⁺] is only ~0.5–1 mM, establishing that the condensate environment must enrich Mg²⁺ to support assembly.

**Physical mechanism**: Mg²⁺ coordinates with the negatively charged phosphate groups of RNA, neutralising electrostatic repulsion and allowing RNA strands to approach for hydrogen-bond formation. The Goldilocks effect (Guth-Metzler et al. 2023) means both too-low and too-high free [Mg²⁺] impair assembly: optimum ~3 mM for tRNA-sized domains, ~10 mM for larger rRNA domains.

**Critical limitation — total versus free [Mg²⁺]**: All available measurements of Mg²⁺ enrichment in condensates (ICP-MS, NMR, partition coefficients) report **total** [Mg²⁺], not **free** [Mg²⁺]. In the RNA-dense FC/DFC, the majority of total Mg²⁺ is sequestered onto RNA phosphate groups, and free [Mg²⁺] — the fraction available for assembly — could be substantially lower than the total and potentially lower than in the GC. No direct measurement of free [Mg²⁺] inside the nucleolus at sub-compartmental spatial resolution exists.

**Gradient direction and magnitude**: The direction of any free [Mg²⁺] gradient within the nucleolus is uncertain. Three scenarios are plausible: (i) Donnan enrichment dominates — free [Mg²⁺] higher at the RNA-dense centre and decreasing outward; (ii) RNA sequestration dominates — free [Mg²⁺] lower at the centre despite higher total; (iii) both centre and periphery are sub-optimal (centre bound to RNA, periphery at cytoplasmic ~0.5–1 mM which is below the Goldilocks peak). The cytoplasmic boundary condition means γ = 1 at r = R_c cannot be anchored to Mg²⁺ reaching optimal assembly concentrations there.

**Connection to model**: Mg²⁺ literature is useful for calibrating the energy scale J (typical r-protein:rRNA contacts are 2–8 kcal/mol; the S4-rRNA Kd shifts 15-fold over 2–8 mM Mg²⁺, setting ΔΔG ≈ 2.7 kBT as a reference scale). However, Mg²⁺ concentration cannot currently be used to derive γ₀ or the gradient shape, and any claim of a Mg²⁺-based spatial gradient within the nucleolus would be speculative.

### 6.3 pH

The pH gradient from pH 6.5 (FC interior) to pH 7.2 (nucleoplasm), directly measured by King et al. (2024), is the best-characterised chemical gradient in the nucleolus. Its effect on assembly must be decomposed by contact type:

**r-protein binding thermodynamics (direct effect, small)**: The dominant RNA-binding residues in r-proteins are Arg (pKa ≈ 12.5) and Lys (pKa ≈ 10.5), which account for ~75% of r-protein:rRNA contacts and are fully protonated throughout 6.5–7.2. Histidine (pKa ≈ 6.0) accounts for ~5–10% of contacts; the correct Henderson-Hasselbalch values are ~6% protonated at pH 7.2 and ~24% at pH 6.5 with pKa = 6.0 (~3-fold change with protein-context pKa ≈ 6.5). Where His+ forms direct salt bridges with RNA phosphates, lower pH actually strengthens binding. The net thermodynamic pH effect on r-protein binding is ~2–5% across the nucleolar range — too small and too ambiguous in sign to drive the γ gradient alone.

**RNA secondary structure stability**: RNA duplexes are modestly more stable at lower pH (~0.43 kcal/mol per pH unit; partial protonation of A and C bases). RNA-RNA contacts are **strengthened** in the acidic interior. In the dual-gradient model (complexes.md), this justifies γ₁ higher at the centre (inverted gradient for RNA-RNA contacts). For the single-γ model, this effect opposes the γ < 1 requirement for the interior.

**Assembly factor activity (primary mechanism for γ₂ < 1 at centre)**: DDX21, an essential DEAD-box helicase that resolves RNA secondary structure kinetic traps during assembly, has a pH optimum of 6.8–7.5 and is sub-optimal at pH 6.5 (Rai et al. 2023). This directly reduces the effective propensity for correct RNA-protein contacts to form in the interior — the most experimentally defensible mechanism linking the pH gradient to γ₂ < 1.

**Functional form**: Using the measured pH profile and DDX21 activity as the reference:

γ_eff(r) ≈ 1 / (1 + exp(−k(pH(r) − pH_ref)))

with pH(r) = 6.5 + 0.7×r/R_c, pH_ref ≈ 7.0 (assembly optimum), k ≈ 2–3. This gives γ_eff(0) ≈ 0.35–0.4 and γ_eff(R_c) ≈ 0.85–0.9, consistent with the default γ₀ = 0.4 and confirming the approximately linear approximation is reasonable.

### 6.4 Salt / Ionic Strength

**Physical mechanism**: Higher ionic strength screens electrostatic interactions following Debye-Hückel theory. For RNA-protein contacts with a net electrostatic contribution:

ΔG_contact(I) ≈ ΔG₀ − N_charges × k_BT × κa

where κ = √(2I/ε₀ε_r k_BT) is the inverse Debye length and a is the contact radius.

**Gradient direction**: Monovalent salt concentration in the nucleolus is not well characterised, but it is likely higher near the osmotically active transcription centres. Higher salt near the centre would weaken contacts (lower γ), consistent with the gradient direction in the model.

### 6.5 Summary Table

| Chemical | Effect on assembly | Gradient direction (for γ) | Status |
|----------|-------------------|---------------------------|--------|
| **pH** | RNA-RNA contacts strengthened at low pH (inverted for γ₁); assembly factor DDX21 impaired at low pH (supports γ₂ < 1 at centre); r-protein Arg/Lys binding pH-independent | Confirmed increasing outward for effective protein-binding propensity (King et al. 2024) | Directly measured; best-supported gradient |
| **Mg²⁺** | Essential for RNA folding and r-protein binding; free [Mg²⁺] gradient within nucleolus unknown | **Unknown**: total enrichment in condensate established; free [Mg²⁺] may be depleted at RNA-dense centre due to sequestration | No direct free [Mg²⁺] measurement at sub-nucleolar resolution |
| **Salt / ionic strength** | Higher ionic strength screens electrostatic r-protein contacts | Uncertain; monovalent ions weakly excluded from condensates | No direct measurement |

γ₀ is treated as a phenomenological parameter (range 0.3–0.5) rather than a quantity derivable from a specific chemical mechanism. The pH gradient (through DDX21 activity) provides the clearest biological motivation for γ₀ < 1.

### 6.6 Quantitative Connection: γ to kBT

With J = 8 k_BT and γ₀ = 0.4, the minimum coupling energy is 0.4 × 8 = 3.2 k_BT. At room temperature (k_BT = 0.593 kcal/mol), this is ~1.9 kcal/mol per contact. A typical r-protein:rRNA interface contact measured by ITC is ~2–8 kcal/mol. J = 8 k_BT sits at the lower end of this range, appropriate for secondary/tertiary binders — the contacts most likely to be disrupted during annealing. The γ₀ = 0.4 minimum (1.9 kcal/mol) represents a plausible partially denaturing regime where contacts can still form but with reduced stability, consistent with rRNA having some secondary structure in the FC while failing to achieve the fully assembled state there. The parameter range γ₀ = 0.3–0.5 should be explored systematically (as recommended in §4) rather than pegged to a specific chemical calculation.

---

## 7. VMMC: Strengths, Limitations, and Appropriateness

### 7.1 Strengths

**Detailed balance by construction**: This is the most important property. VMMC guarantees that the equilibrium distribution is the correct Boltzmann distribution for the defined Hamiltonian. Any observable computed from long VMMC trajectories converges to its true equilibrium value — something that cannot be guaranteed for many popular coarse-grained MD approaches.

**Cluster moves**: VMMC naturally handles the diffusion of multi-particle clusters. A 16-particle assembled complex diffuses as a rigid body under VMMC cluster translation, which is physically correct (diffusion of a rigid body with diffusion constant D ∝ 1/R under Stokes drag). Simple single-particle Monte Carlo cannot do this.

**Efficiency for assembly**: VMMC was originally designed for systems where assembly intermediates appear and disappear, and it handles the transition between single-particle and cluster states without special-casing.

**Computational efficiency**: For a 2D lattice model with N = 128 particles (8 copies × 16) in a circle of radius 60 (area ≈ 11,300 lattice sites), the system is very sparse. VMMC with cell lists scales as O(N) per iteration, making 200,000 iterations computationally fast (~minutes on a single core).

### 7.2 Limitations

**Lattice discretisation**: The square lattice restricts translational moves to 4 cardinal directions (or 8 if diagonal moves are included) and orientational states to 4 discrete values. This underestimates the orientational entropy of the real system, which has continuous orientations. The discrete patch model correctly gates bonds in 4 directions but cannot represent bonds at intermediate angles.

**No genuine time mapping**: VMMC iterations are not physical time steps. The relationship between VMMC iterations and real physical time is complex and depends on system-specific acceptance rates, cluster sizes, and move frequencies. This means you cannot directly compare to the 10 nm/s drift measured in pulse-chase experiments without additional calibration.

To calibrate: measure the mean squared displacement of a free monomer under VMMC and compare to the expected D for a ribosomal protein (~1 μm²/s from FRAP). This would allow a conversion of VMMC iterations to physical time and permit comparison with experimental timescales.

**Rotation fraction φ_rot as free parameter**: The probability of rotation moves (φ_rot = 0.2) is a free parameter with no direct physical meaning. It affects the rate at which misassembled orientational states are explored. In the SED (Stokes-Einstein-Debye) picture, rotational diffusion of a sphere of radius r has D_rot = k_BT/(8πηr³), while translational diffusion has D_trans = k_BT/(6πηr). For a monomer of radius r ≈ 5 nm, D_rot/D_trans ≈ 1/(4r²) ≈ 10⁻² nm⁻² — meaning rotational relaxation is much faster than translational diffusion. The value φ_rot ≈ 0.43 mentioned in the rotation_report.md is thus more physically justified than the default 0.2; running at 0.43 should improve assembly.

**VMMC cluster rotation is not fully realistic**: The cluster rotation move rotates all bonded particles by 90° around a chosen centre. For a 16-particle complex, rotation by exactly 90° is a very specific move and may not represent the orientational diffusion of a real compact body. For large clusters, rotational moves become increasingly rare and may lead to orientational traps where the entire cluster is slightly misaligned.

**No memory of hydrodynamic interactions between particles**: VMMC treats each particle's diffusion independently, not accounting for hydrodynamic correlations between nearby particles (which would create correlated diffusion in real systems). For dilute systems this is a reasonable approximation.

### 7.3 Verdict

VMMC is the right algorithm for this system. The key alternatives (simple Metropolis, MD, Brownian dynamics) are all worse: Metropolis cannot handle clusters, MD does not enforce detailed balance in coarse-grained settings, and BD at the coarse-graining level of this model would violate detailed balance when clusters form. The known issues (lattice discretisation, φ_rot calibration, time mapping) are all manageable limitations that should be discussed but do not invalidate the conclusions.

---

## 8. Model Improvements for a More Realistic Ribosomal Model

In roughly increasing order of implementation complexity:

### 8.1 Short-term (Feasible in Weeks)

1. **More independent seeds** (at least N=20 per condition) — already discussed in §4.
2. **φ_rot = 0.43** — use the SED-motivated value throughout.
3. **Non-linear gradient shapes** — modify `gamma_r()` in CondensateModel.cpp to support power-law and sigmoidal profiles.
4. **Gradient fluctuations** — add time-dependent noise to γ with tunable amplitude and correlation time.
5. **Asymmetric complex** — modify TARGET_X/Y to create an asymmetric complex (not all polymers identical in length), testing whether the gradient helps more with asymmetric assembly.

### 8.2 Medium-term (Feasible in Months)

6. **3D lattice model** — extend to cubic lattice with 6 faces per particle. This is a significant code refactor but would eliminate the "2D is unrealistic" criticism.
7. **Larger complex (n=3 Moore curve, 36 particles)** — a larger complex has more misassembly intermediates and more potential for kinetic frustration, making the annealing effect more dramatic.
8. **Heterogeneous J values** — assign different J values to different contact types, reflecting the known variation in ribosome assembly factor binding energetics.
9. **Scaffold particle** — introduce a large central "rRNA" particle (occupying a 4×4 or 6×6 square) that acts as the assembly scaffold, with protein-like particles binding to its periphery. This is much more biologically faithful.
10. **Calibrated time mapping** — measure D_monomer under VMMC and convert iterations to physical time.

### 8.3 Long-term (Major Extension)

11. **Off-lattice continuous VMMC** — use the original Whitelam–Geissler formulation with continuous positions and Kern–Frenkel patches at user-defined angular width δ. This removes all lattice artefacts and allows genuine orientational entropy.
12. **RNA flexibility** — allow the backbone chain to flex (modify backbone bonds from rigid to harmonic), producing more realistic polymer dynamics.
13. **Active advection** — add an explicit drift velocity field v(r) pointing radially outward, in addition to the Stokes drag. This would directly model the cytoplasmic streaming and/or rDNA transcription pressure that drives particle outflow.

---

## 9. New Ideas for Impressive Publications

### 9.1 The "Assembly Phase Diagram" Paper

**Idea**: Produce a 2D phase diagram of exitQuality as a function of (γ₀, J). At fixed R_c and N_copies, sweep over a grid of γ₀ ∈ [0, 1] × J ∈ [2, 20] with 20 seeds per point. The result will show an "assembly funnel" in parameter space: a region of (γ₀, J) where good assembly occurs, surrounded by regions of either no assembly (J too small, contacts don't form) or frustrated assembly (J too large, kinetic traps dominate, or γ₀ too high, no annealing).

**Why impressive**: This is a clean, falsifiable quantitative prediction. The claim is that the nucleolus operates within the assembly funnel, and the Mg²⁺ gradient controls which (γ₀, J) point is accessed. Real Mg²⁺ concentrations map to J values through the Hill function, and the optimal Mg²⁺ range from Guth-Metzler et al. (2023) should fall within the assembly funnel. This would constitute quantitative, testable biophysics.

### 9.2 Spatial vs. Temporal Annealing: Is the Nucleolus Optimal?

**Idea**: Compare three conditions, all at the same total "annealing budget" (same range of γ):
- (A) Spatial gradient traversed by diffusing+advecting particles (current model)
- (B) Temporal linear ramp γ(t) = γ₀ + (1−γ₀)t/T in a box (run_box --anneal)
- (C) Optimal temporal schedule from Trubiano & Whitelam (2022) applied to the same complex

The question: does the spatial gradient in the condensate geometry outperform a simple temporal ramp? If yes, is it because the condensate geometry automatically sizes the time spent in each coupling regime to the kinetic bottleneck? And how close does it come to the theoretically optimal schedule?

**Why impressive**: This directly addresses whether the nucleolus has been evolutionarily optimised, and provides a theoretical bound on how much better a perfectly designed gradient could be. If the spatial gradient approximately matches the optimal temporal schedule, that is a remarkable quantitative result suggesting the nucleolus geometry implements thermodynamic optimality.

### 9.3 Advective Sorting: Large Complexes Exit First

**Idea**: Explicitly demonstrate the Stokes-drag-mediated sorting effect. Track the size distribution of particles at the exit boundary as a function of time. The prediction is that assembled 16-particle complexes (D ∝ 1/4) exit later but with higher quality than free monomers (D ∝ 1) which exit early but with zero quality. There is an intermediate regime of partial complexes (size 2–12) that exit with intermediate quality.

Produce a scatter plot of exit time vs. exitQuality for individual exit events. This would show a clear trend: late-exiting particles are more fully assembled, which is a direct model prediction of the Stokes mechanism.

**Why impressive**: This connects directly to the pulse-chase experiments of Riback et al. (2023) where different maturation states exit at different times. A quantitative match to the observed size-sorting would be a strong experimental validation of the model.

### 9.4 Information Stored in the Gradient

**Idea**: Quantify how much "information" the spatial gradient encodes about the assembly pathway, using the framework of Trubiano & Whitelam. The gradient specifies a trajectory through (r, γ) space; this trajectory defines a path in the free energy landscape. For a given complex geometry and coupling matrix, how much Fisher information is available to the assembling particle about which path to take to the correctly assembled state?

This is a more theoretical paper connecting information theory to condensate biology — closer to the Sivak–Crooks framework.

### 9.5 Gradient-Mediated Error Correction

**Idea**: After a complex assembles in the periphery (γ ≈ 1) and then circulates back to the centre (γ ≈ γ₀), would a misassembled complex be more likely to be corrected than a correctly assembled one? If misassembled contacts have higher energy than native contacts (they involve non-complementary patches), they would melt more easily in the low-γ interior.

This is a model of **gradient-mediated error correction**: the condensate geometry not only assists assembly but also selectively destroys errors. This would be a very striking biological implication — the nucleolus as a quality control machine, not just an assembly machine.

To test: initialise some copies as misassembled complexes (wrong patches in contact) and measure their exit quality versus time compared to correctly assembled complexes. Do misassembled complexes get corrected more often when they visit the interior?

---

## 10. Honest Conclusion: Publishability Assessment

### 10.1 Is the Core Idea Novel and Important?

Yes, with caveats. The idea that spatial gradients in the nucleolus function as annealing schedules for ribosome assembly is not, to my knowledge, explicitly modelled in the literature. The experimental papers (Riback et al. 2023; Lafontaine et al. 2021; Detres et al. 2025; Guth-Metzler et al. 2023) contain all the ingredients — radial transit, chemical gradients, need for Mg²⁺ ramp — but none has connected them to the annealing concept from statistical mechanics or demonstrated the connection with a simulation.

However, the idea is *close* to things that have been published:
- Trubiano & Whitelam (2022) demonstrated optimal temporal annealing for colloidal assembly
- The concept of the nucleolus as a "production line" (Lafontaine et al. 2021) is established
- Patchy particle assembly with temperature-dependent yield is well-known (Wilber et al. 2007)

The novelty is in the synthesis: connecting these three threads in a single, minimal model that captures the geometry of the nucleolus. That synthesis is real and publication-worthy.

### 10.2 Where Could This Be Published?

**Realistic targets, roughly in order of reach**:

1. **eLife** (biophysics/computational biology) — This is the right level if you have proper statistics, a quantitative connection to biology, and the spatial vs. temporal comparison. eLife rewards conceptual novelty combined with clean, well-controlled simulations. The patchy particle model, circular condensate geometry, and exitQuality metric together constitute a complete study if done properly.

2. **PLOS Computational Biology** — A more certain home. Strong computational biology methodology with biological implications. Less risk of rejection on grounds of "insufficient experimental validation."

3. **Physical Review E** — Good home if the paper is framed as a statistical physics result (optimal control, annealing theory) with biology as the motivation. Harder to get the biology message out.

4. **Journal of Physical Chemistry B** — Strong option for the chemistry-focused angle (Mg²⁺ gradient maps to γ, etc.). Biophysics and soft matter community.

5. **Soft Matter** (RSC) — If framed as a patchy particle / condensate paper, with annealing as the key result.

6. **PNAS** — Possible if you can demonstrate experimental validation (collaboration with a lab measuring nucleolar assembly quality as a function of Mg²⁺ gradients), but very unlikely with simulations alone.

**Not realistic**: Nature, Science, or Cell-family journals. These require either breakthrough experimental validation or theoretical frameworks that overhaul a field. A well-executed minimal simulation, however elegant, is not at that level.

### 10.3 What Are the Fatal Flaws?

There are no **fatal** flaws in the scientific logic, but there are **critical insufficiencies** in the current state of the work:

1. **Insufficient statistics** — The results described (single seed=2 runs) are illustrations, not data. Referees will demand 20+ independent runs per condition with error bars.

2. **Flawed comparison** — The specific comparison described in the preamble (J=8,γ₀=0.4 vs J=80,γ₀=0.9) conflates three variables. This must be fixed before submission or it will be fatal in peer review.

3. **No quantitative connection to biology** — The claim that "annealing can emerge from gradients in the ribosome" requires at least a qualitative mapping from γ to a real chemical parameter (Mg²⁺ is the most defensible). Without this, the model is a mathematical curiosity, not biophysics.

4. **2D limitation** — This is a genuine weakness that must be acknowledged and justified. The justification is computational tractability and the precedent of 2D Go̊ models in the literature; referees will accept this if stated honestly rather than hidden.

5. **VMMC time calibration** — Without a conversion from iterations to physical time, claims about "assembly timescales" relative to the 10-minute residence time measured in vivo are speculative.

### 10.4 What Would Make This a Strong Paper?

A strong paper would contain:

- A systematic parameter scan (exitQuality vs γ₀ at fixed J, with ~20 seeds per point) showing an optimal γ₀
- The three-way comparison: gradient vs. uniform γ=1 vs. uniform γ=γ₀ (or temporal ramp)
- Spatial statistics showing where in the condensate assembly occurs and where traps form
- A quantitative mapping of γ to Mg²⁺ concentration using the Hill function and Guth-Metzler data
- A discussion of how the nucleolus's observed [Mg²⁺] range maps to the assembly funnel
- A comparison with the Trubiano-Whitelam optimal schedule
- Honest discussion of 2D limitations and a simple check (run a 3D version for one parameter set, or cite dimensionality studies)

With these elements, a paper in *eLife* or *PLOS Computational Biology* is achievable. Without them, it remains a proof-of-concept demonstration that belongs in a thesis, not a journal.

### 10.5 A Personal Note

The PhD work underlying this is strong. The observation that a spatial chemical gradient can function as an annealing schedule is physically insightful and deserves to be published. The patchy particle extension is a genuine improvement over the original model. The circular condensate geometry is elegant. What is needed is not more clever ideas but more systematic, statistically rigorous execution of the simulations that already exist. Get the statistics right, control the comparisons properly, and connect γ to a real chemical parameter — and you have a paper worth publishing.

---

## Appendix: Recommended First Simulations (Run These Next)

```bash
# Control: no gradient, full coupling (uniform γ=1)
for seed in $(seq 1 20); do
  ./run_condensate \
    --copies 8  --t-equil 500  --t-denat 1000 \
    --steps 200000  --snapshots 200 \
    --radius 60  --J 8  --gamma0 0.0 \
    --stokes  --coupling midpoint \
    --phi-rot 0.43  --phi-reorient 0.2 \
    --seed $seed  --output ctrl_uniform_$seed &
done
wait

# Treatment: radial gradient, gamma0=0.4
for seed in $(seq 1 20); do
  ./run_condensate \
    --copies 8  --t-equil 500  --t-denat 1000 \
    --steps 200000  --snapshots 200 \
    --radius 60  --J 8  --gamma0 0.4  --gradient \
    --stokes  --coupling midpoint \
    --phi-rot 0.43  --phi-reorient 0.2 \
    --seed $seed  --output grad_g04_$seed &
done
wait

# Control: no gradient, low coupling (uniform γ=0.4 — not directly available
# in condensate model; approximate with very large radius so gradient is flat)
# Better: use run_box with fixed gamma=0.4 to show that low uniform coupling also fails.
```

Then: compute mean ± std of `exitQuality` from the stats files for the final 50,000 steps of each run. Report the three conditions side by side. This is the Figure 1 of the paper.
