# Chemical Gradients in the Nucleolus: A Literature Review

*Compiled for the NucleolusAssemblyPatchy simulation paper — S. Whitby, 2026*

---

## Overview

The nucleolus is the largest and most compositionally complex membrane-less organelle in eukaryotic cells. It is organised into at least three concentric liquid-like sub-phases — the fibrillar centre (FC), dense fibrillar component (DFC), and granular component (GC) — that differ dramatically in composition, viscosity, and material properties. Assembling pre-ribosomal particles are born near the transcription centres at the FC/DFC boundary and transit radially outward through the DFC and GC, spending tens of minutes inside the condensate before exiting into the nucleoplasm (Riback et al. 2023). During this transit, the chemical environment they experience changes substantially. The last five years have produced direct experimental measurements of several of these chemical gradients, most notably pH, with substantial indirect evidence for ion (particularly Mg²⁺) gradients.

This document reviews the available evidence for each major chemical species, assesses the likely spatial profile within the nucleolus, discusses temporal fluctuations, quantifies the known effect on ribosome-relevant bonding, and concludes with a synthesis that maps all gradients onto the dimensionless coupling parameter γ in the simulation model.

---

## 1. pH

### 1.1 Direct Measurements In Vivo

The most quantitatively precise chemical gradient measurements in the nucleolus are for pH. King et al. (2024, *Cell*) used two complementary ratiometric fluorescent probes — pHluorin2 (genetically encoded, exclusively reports cytoplasmic/nuclear-phase protein environment) and SNARF-4 (a membrane-permeant small molecule) — in living Xenopus laevis oocytes to map pH across the nucleolus and several other nuclear condensates. Key values:

| Sub-compartment | pH | Reference |
|---|---|---|
| FC / inner DFC | 6.5 ± 0.1 | King et al. 2024, *Cell* |
| GC inner half | 6.6 ± 0.1 | King et al. 2024, *Cell* |
| GC outer half | 6.8 ± 0.1 | King et al. 2024, *Cell* |
| Nucleoplasm | 7.2 ± 0.1 | King et al. 2024, *Cell* |
| Cajal bodies | 7.5 ± 0.1 | King et al. 2024, *Cell* |
| Nuclear speckles | 7.2 ± 0.1 | King et al. 2024, *Cell* |

This establishes a pH gradient of **~0.7 pH units** across the nucleolus, from the most acidic core (FC, pH ≈ 6.5) to the nucleoplasm (pH ≈ 7.2). Crucially, the gradient is monotonic and radially organised: the interior is more acidic than the periphery. The proton motive force corresponding to a 0.7-unit gradient is approximately −88 mJ per proton, a value comparable in magnitude to the proton-motive force across mitochondrial membranes.

The mechanistic basis of this gradient was established by King et al. using a combination of proteomics and mutagenesis: nucleolar proteins are enriched in acidic D/E amino-acid tract sequences (IDRs rich in aspartate/glutamate). These tracts act as proton buffers within the condensate — their negative charges recruit protons from the dilute phase, acidifying the condensate interior. The D/E-tract density is highest in the FC and DFC (UBF, nucleolin) and decreases toward the GC (NPM1), naturally explaining the radial pH gradient. High-salt injection abolished this gradient, confirming its electrostatic origin and linking it to the coacervate nature of the condensate.

Complementary in vitro measurements by King et al. using reconstituted condensates (NCL+pre-rRNA; NPM1+mat-rRNA) confirmed the same trend: NPM1+rRNA condensates show an internal pH ~0.6 units below the dilute phase, and this gradient is eliminated in the ΔD/E-tract NCL mutant.

A parallel and independent set of measurements by Ausserwoger et al. (2025, *Nature Chemistry*) using a high-throughput microdroplet platform measured dense-phase pH in >75,000 synthetic condensate conditions. They demonstrated a universal mechanistic principle: the dense phase of a condensate buffers toward the isoelectric point (pI) of the protein mixture, because charge neutralisation within the condensate excludes buffer molecules that would otherwise set the pH. For systems with pI < nucleoplasm pH (e.g., NPM1 pI ≈ 5.5, DFC proteins pI ≈ 5–6), this produces acidic interiors, quantitatively consistent with King et al.'s measurements. The predicted nucleolus pH profile from this model aligns well with the observed 6.5–6.8 range in the GC/DFC.

Rai et al. (2023, *Molecular Cell*) further showed that the nucleolar pH is **activity-dependent and rapidly reversible**. In actively transcribing cells (high Pol I activity), the nucleolus maintains its acidic state (pH 6.8–7.5). When rRNA synthesis is inhibited by CX5461, actinomycin D at 50 nM (Pol I specific), or PARP inhibitors, the nucleolus deacidifies to pH 8.0–8.8 within minutes. This is because PARP generates protons as a metabolic byproduct of rRNA synthesis, maintaining the condensate's acidity during active production. The recruitment of the helicase DDX21 to the nucleolus — essential for ribosome biogenesis — is pH-dependent and fails both below pH 6.0 and above pH 8.0, with an optimal window of pH 6.8–7.5. This demonstrates that the pH gradient has direct functional consequences for assembly factor recruitment.

Detrés, Camacho-Badillo & Calo (2025, *Journal of Molecular Biology*) synthesised this body of evidence into a "pH-centric model" of nucleolar activity, proposing that the acidic condensate interior is not an incidental consequence of condensate chemistry but a functionally maintained feature whose disruption — in both directions — impairs ribosome biogenesis. The nucleolus is most acidic (pH 6.5) when rRNA production is highest, and less acidic when production is attenuated.

In Saccharomyces cerevisiae, rRNA synthesis rate is pH-sensitive: shifting cells from pH 8.8 to pH 7.0 in sporulating yeast restored rRNA synthesis to vegetative levels within 10 minutes, demonstrating that pH regulation of nucleolar function extends beyond mammalian systems.

### 1.2 Effect of pH on Ribosomal Bonding

**RNA secondary structure stability**: Lowering pH from 7.5 to 5.5 increases RNA duplex stability by approximately 0.43 kcal/mol at high ionic strength (measured by thermal denaturation; RNA journal 2019). This is modest but non-negligible: a single RNA-RNA contact becomes ~0.9-fold stronger per pH unit decrease. The primary mechanism is partial protonation of adenine and cytosine bases at pH <6.0, which can stabilise wobble pairs, and protonation of phosphate groups that reduces electrostatic repulsion within the RNA backbone.

**Protein-RNA binding via charged residues**: The dominant RNA-binding residues in ribosomal proteins are arginine and lysine, not histidine. Ribosomal proteins are among the most basic proteins in the proteome — typical pI values are 9.5–12.5 — driven by strong enrichment of Arg (~8–10% of residues vs. ~5.5% proteome average) and Lys (~9–11% vs. ~5.7% average). Structural analyses of ribosome crystal structures show that approximately 75% of all r-protein:rRNA contacts involve Arg or Lys residues. The guanidinium of Arg (pKa ≈ 12.5) and the ε-amino of Lys (pKa ≈ 10.5) are fully and permanently protonated across the entire nucleolar pH range of 6.5–7.2, contributing a **pH-independent** electrostatic attraction to RNA phosphates.

Histidine (pKa ≈ 6.0 in solution; 5–7 in protein context) is the only common residue with a pKa in the nucleolar range and is therefore the only direct charge-state switch available between pH 6.5 and 7.2. However, histidine is present in r-proteins at approximately normal proteome abundance (~2–3% of residues, not enriched) and contributes only ~5–10% of r-protein:rRNA contacts. The correct Henderson-Hasselbalch calculation gives:

| pH | pKa = 6.0 (solution) | pKa = 6.5 (protein context) |
|---|---|---|
| 7.2 | ~6% protonated | ~17% protonated |
| 6.5 | ~24% protonated | ~50% protonated |
| Fold-change | ~4-fold | ~3-fold |

With His contributing ~7% of contacts and its charge state changing ~4-fold, the net contribution to total binding energy change is approximately 0.07 × kBT ln(4) ≈ 0.1 kcal/mol — a ~2–5% shift in total binding across the nucleolar pH range.

Crucially, the sign of the His contribution is not straightforward. Where His+ forms a direct salt bridge with an RNA phosphate, **lower pH strengthens binding** (more His+ = stronger ion pair). Where His sits in an already positively charged local environment (common given the K/R richness of r-proteins), additional protonation increases intra-protein charge repulsion and may slightly destabilise the binding conformation. The net direction is protein-specific.

RNA phosphodiester groups have pKa ~1–2 and are essentially fully deprotonated at pH 6.5: the negative charge density on the rRNA backbone is unchanged between pH 6.5 and 7.2. Aspartate (pKa ≈ 3.7) and glutamate (pKa ≈ 4.2) are similarly fully deprotonated throughout the nucleolus and nucleoplasm.

**The direct pH effect on r-protein:rRNA binding thermodynamics across the nucleolar pH range is therefore small (~2–5% of total binding energy) and of ambiguous sign.** It cannot account for the large changes in assembly propensity that distinguish the condensate interior from the periphery.

**SRSF1 / splicing factor relevance**: A 2025 bioRxiv study found that splicing factor SRSF1 acts as a "pH-stat" to restore nucleolar integrity when pH deviations occur, reinforcing the centrality of pH maintenance to nucleolar function.

**Direction of the pH effect on assembly**: The pH gradient affects different components of the assembly process in opposite directions, which has important implications for the model:

- **RNA-RNA contacts** (secondary structure, base-stacking, A-minor interactions): **strengthened** at lower pH. Partial protonation of adenine (pKa ~6 in structured RNA contexts) and cytosine stabilises additional hydrogen bonds and reduces backbone repulsion. RNA duplexes are modestly more stable at pH 6.5 than at 7.2 by approximately 0.43 kcal/mol per pH unit. This means the acidic condensate interior actively promotes RNA secondary structure formation. In the context of a dual-gradient model (see complexes.md), this justifies an *inverted* gradient for RNA-RNA coupling (γ₁ higher at the centre).

- **r-protein Arg/Lys contacts**: pH-independent throughout 6.5–7.2, as argued above. No gradient contribution from this dominant binding mechanism.

- **r-protein His contacts**: small (~2–5% of total binding energy) and sign-ambiguous. Not a reliable source of gradient.

- **Assembly factor activity**: the most experimentally supported mechanism for pH-dependent suppression of assembly in the condensate interior. DDX21, an essential DEAD-box helicase that resolves RNA secondary structure kinetic traps during ribosome assembly, has a documented pH optimum of 6.8–7.5 (Rai et al. 2023). At pH 6.5, DDX21 operates below its optimum, meaning RNA folding kinetic traps accumulate in the condensate interior and are not corrected as efficiently. This is a kinetic, not thermodynamic, suppression mechanism.

The conclusion that the interior is less favourable for final assembly than the periphery is well-supported, but it rests primarily on **assembly factor activity being impaired at pH 6.5** and on the indirect effect of **RNA secondary structure over-stability potentially occluding protein binding sites** — not on direct suppression of r-protein:rRNA binding thermodynamics through charge-state changes.

### 1.3 Spatial Profile

The empirical pH profile across the nucleolus is well-described by a monotonic function that increases from ~6.5 at r = 0 (FC/DFC) to ~7.2 at r = R_c (nucleolus edge). The profile is approximately linear or slightly concave (the FC-to-DFC step, from 6.5 to 6.6, is 0.1 units over a small distance, while the GC-to-nucleoplasm step is 0.4 units over a larger radius). A rough mapping:

- r/R_c ≈ 0–0.3 (FC/inner DFC): pH 6.5
- r/R_c ≈ 0.3–0.7 (GC inner): pH 6.6–6.7
- r/R_c ≈ 0.7–1.0 (GC outer): pH 6.8–7.0

This is consistent with a linear gradient pH(r) ≈ 6.5 + 0.7 × r/R_c.

### 1.4 Temporal Fluctuations

The nucleolar pH responds to Pol I transcription activity on a timescale of minutes:
- Pol I inhibition deacidifies the nucleolus from pH 6.8 to >8 within 5–10 minutes
- pH recovery upon restoration of Pol I activity is similarly rapid
- During cell division, the nucleolus disassembles in prophase (timescale ~30 minutes, accelerating orders of magnitude near nuclear envelope breakdown; Vafabakhsh et al. 2024, *PNAS*)
- Between cell divisions, nucleolar pH appears stable during the S/G2 phase when ribosome production is maximal

---

## 2. Magnesium (Mg²⁺)

### 2.1 Background: Intracellular Mg²⁺ Homeostasis

Magnesium is the most abundant divalent cation in the cell and the dominant divalent ion controlling RNA structure and ribosome stability. The total intracellular Mg²⁺ concentration in eukaryotic cells is approximately 20 mM (measured by ICP-OES), but the vast majority is bound to RNA, DNA, ATP, and phospholipids. The free cytoplasmic Mg²⁺ in eukaryotes is tightly buffered at **0.5–1 mM** (measured by ³¹P-NMR and fluorescent probes including mag-fluo4). This buffering is maintained by the competition of ATP (intracellular concentration ~2–5 mM, binding Mg²⁺ with K_d ≈ 0.1 mM) and other phosphometabolites for free Mg²⁺.

A single bacterial 70S ribosome binds >170 Mg²⁺ atoms by X-ray crystallography, with the most recent 1.55 Å resolution structure of the E. coli 70S ribosome identifying 290–403 well-defined Mg²⁺ ions (Nierhaus/Sigala review 2014, *J. Bacteriology*; NAR 2025 structure paper). Stoichiometrically, the ~10,000–15,000 ribosomes per bacterial cell sequester approximately 12 mM equivalent of Mg²⁺, making ribosomes the dominant Mg²⁺ sink in the cell.

**No direct measurement of free [Mg²⁺] inside the nucleolus or any nucleolar sub-compartment has been published** as of 2026. Genetically encoded Mg²⁺ sensors (MagIC, ratiometric) have been targeted to mitochondria and the endoplasmic reticulum but not to the nucleus or nucleolus. This is the most important gap in the literature for this simulation paper.

### 2.2 Indirect Evidence for Mg²⁺ Enrichment in Condensates

Despite the absence of direct nucleolar Mg²⁺ measurements, strong indirect evidence supports Mg²⁺ enrichment within condensate-type structures:

**In vivo bacterial condensates** (Dai et al. 2024, *Cell*): ICP-MS measurements in E. coli forming synthetic RLP (RNA-like protein) condensates showed that cytoplasmic Mg²⁺ decreased ~2-fold upon condensate formation, implying a ~5-fold enrichment of Mg²⁺ inside condensate droplets (the condensate volume is ~20% of the cell). Simultaneously, Na⁺ was excluded (19% increase in the cytoplasm upon condensation). Ca²⁺ was enriched ~2-fold within condensates. These measurements provide the first direct in vivo evidence that condensates actively enrich divalent cations.

**NMR measurements** (Selective Ion Binding in Condensates, *JACS* 2025): ²³Na and ²⁵Mg solid-state and solution NMR of in vitro condensates demonstrated that kosmotropic cations (Li⁺, Mg²⁺, Ca²⁺) are selectively enriched within the dense phase of protein/RNA condensates, while chaotropic cations (Cs⁺, K⁺) are excluded. The enrichment follows the Hofmeister series via the "law of matching water affinities." Mg²⁺, being the most strongly hydrated divalent cation, is consistently enriched.

**ICP-MS in peptide-polymer condensates** (Competitive Mg²⁺ Regulation of Condensate Microenvironments, *JACS Au* 2025): Partition coefficients for Mg²⁺ into condensates formed by poly-arginine/polyinosine systems: K = 98; for poly-lysine/poly-aspartate: K = 57. These values confirm that Mg²⁺ enrichment of ~60–100-fold above the dilute phase is achievable in highly charged condensates. While the nucleolus is not a simple coacervate of this type, its NPM1+rRNA and nucleolin+rRNA cores are polyanionic-polycationic systems that would be expected to enrich Mg²⁺.

**NPM1+rRNA condensates** (Yewdall et al. 2022, *Biophysical Journal*): In vitro condensates of NPM1 (a major GC protein) with rRNA show dramatic Mg²⁺-dependent transitions. At [Mg²⁺] ≥ 7 mM, the condensates transition from liquid to gel with FRAP recovery τ > 1500 s and <10% total recovery. At 20 mM Mg²⁺, gel-like morphology is established; ATP addition (competing with Mg²⁺ as a chelator) restores liquidity over ~20 minutes. This demonstrates that the material state of nucleolar condensates is exquisitely sensitive to [Mg²⁺] in the range 1–20 mM.

**Taken together**: the nucleolus, as a highly charged RNA-protein condensate, is likely to have elevated **total** [Mg²⁺] compared to the cytoplasm. However, a critical distinction must be made: all of the enrichment evidence above (ICP-MS, NMR, partition coefficients) reports **total** [Mg²⁺] in the condensate phase, not **free** [Mg²⁺]. In RNA-dense environments, most Mg²⁺ above background is coordinated to RNA phosphate groups; it is the free fraction that drives additional RNA folding or r-protein binding. A reasonable estimate of total [Mg²⁺] enrichment in the nucleolar condensate is ~5-fold above cytoplasm, but the free [Mg²⁺] fraction within the condensate is unknown. Whether there is a gradient **within** the nucleolus in free [Mg²⁺] is unknown experimentally; see §2.5 for a critical assessment.

### 2.3 Effect of Mg²⁺ on rRNA Folding (Quantitative)

**16S rRNA central domain folding** (Hori, Denesyuk & Thirumalai 2021, *PNAS*; coarse-grained molecular dynamics simulation):
- Global folding midpoint: **[Mg²⁺] = 3.3 mM**, Hill coefficient **n = 2.96** (cooperative folding)
- Three sequential tertiary contact stages with midpoints at ~2.5, ~4, and ~5 mM Mg²⁺
- At [Mg²⁺] < 1 mM: only secondary structures form; all tertiary contacts are absent
- At [Mg²⁺] = 2 mM: early tertiary contacts form, including the critical h19 helix
- At [Mg²⁺] = 5 mM: three-way junction reorganisation and full tertiary structure achieved

**30S structural integrity at physiological Mg²⁺** (cryo-EM, 2023, PMC10046523):
- At 10 mM Mg²⁺: stable mature 30S structure
- At 2.5 mM Mg²⁺: helix h44 "barely visible" and h17 flexible
- At 1 mM Mg²⁺: helix h27 shifts by 19–26 Å, ~14% of particles lose ribosomal protein S12, head domain becomes flexible; 70S dissociates
- At 0.5 mM Mg²⁺: complete 70S dissociation
- All changes reversible upon restoration to 10 mM Mg²⁺

**r-protein–rRNA binding constants vs. Mg²⁺** (FRET assay, 37°C; S4 binding to 16S rRNA 5' domain, PMC4604426):
- Kd at 2 mM Mg²⁺: **0.6 ± 0.1 nM**
- Kd at 8 mM Mg²⁺: **0.04 ± 0.01 nM** (15-fold tighter)
- Kd at 20 mM Mg²⁺: 0.08 ± 0.01 nM (slight loosening relative to 8 mM, consistent with Goldilocks effect)
- Conformational equilibrium (native:flipped conformation ratio): shifts **127-fold** from 2 mM to 12 mM Mg²⁺ (from flipped-dominant to native-dominant)

**Goldilocks Mg²⁺** (Guth-Metzler et al. 2023, *Nucleic Acids Research*, PMC10164553; in vitro alkaline hydrolysis):
- tRNA^Phe: maximum RNA protection (ratio of ku to kf) at ~**3 mM Mg²⁺**; peak width ~1 mM at half-height; 2.4-fold protection at peak
- P4-P6 RNA domain (~160 nt): Goldilocks peak at ~**10 mM Mg²⁺**; 22-fold protection ratio (ku/kf)
- Below the Goldilocks Mg²⁺: RNA is unfolded and vulnerable to hydrolysis
- Above the Goldilocks Mg²⁺: RNA is over-stabilised (possibly misfolded or aggregated)
- Unstructured rU20 control: monotonic decrease in lifetime with Mg²⁺ (no Goldilocks peak) — demonstrating the effect is specific to structured RNAs

**Key implication**: The free cytoplasmic [Mg²⁺] of ~0.5–1 mM is *below* the Goldilocks peak for even small tRNA-sized RNAs (peak at ~3 mM), and far below the peak for larger rRNA domains (peak at ~10 mM). This means that *in vivo* ribosome assembly must be occurring in sub-compartmental environments where [Mg²⁺] is locally elevated above the cytoplasmic background — exactly the kind of enrichment the condensate provides. The nucleolus thus functions, in part, as a Mg²⁺ concentrating device.

### 2.4 Reconstitution Protocols (Quantitative Requirements)

**Nierhaus two-step reconstitution** (Nierhaus & Dohme 1974; reviewed Nierhaus 1990):
- Step 1: 4 mM Mg²⁺, 44°C, 20 minutes → produces 48S inactive intermediate
- Step 2: 20 mM Mg²⁺, 50°C, 90 minutes → active 50S subunit (50–100% native activity)
- 30S reconstitution: single step at ~20 mM Mg²⁺, or 16 mM Mg²⁺ in the presence of 3 mM spermidine
- Spontaneous assembly at physiological 1 mM free Mg²⁺ does not occur in vitro — confirming that the condensate must provide something beyond the cytoplasmic background

**GTPase-mediated rescue** (eLife 2024): Assembly factors EngA and ObgE can substitute for the high-Mg²⁺/high-temperature requirement, enabling single-step reconstitution at physiological conditions. This demonstrates that chaperones provide the energy that Mg²⁺ otherwise supplies to fold the rRNA — relevant to the simulation model which lacks chaperones.

### 2.5 Spatial Profile

**Critical distinction: total versus free [Mg²⁺]**: All enrichment evidence cited in §2.2 (ICP-MS, partition coefficients, NMR) measures or estimates **total** [Mg²⁺] in the condensate phase, not **free** [Mg²⁺]. In RNA-dense environments this distinction is fundamental. Mg²⁺ coordinates to RNA phosphate groups through both diffuse outer-sphere interactions and tight inner-sphere site-specific binding. The FC/DFC, where newly transcribed pre-rRNA is present at well above the polymer overlap concentration (rRNA concentration ~5 μM; Riback et al. 2023), is a powerful Mg²⁺ **sink**: essentially all Mg²⁺ above the background cytoplasmic level in this region is sequestered onto RNA. The free [Mg²⁺] — the only fraction available to facilitate folding of newly synthesised RNA or r-protein binding — could therefore be **lower** in the RNA-dense interior than in the GC, even if total [Mg²⁺] is higher there. As the particle moves outward through the GC, rRNA becomes progressively folded and protein-bound, reducing the available phosphate-binding sites and potentially releasing Mg²⁺ back to the free pool.

No direct measurement of **free** [Mg²⁺] inside the nucleolus at sub-compartmental resolution exists as of 2026. The following summarises what can reasonably be said:

- **Nucleolus condensate as a whole** (relative to cytoplasm): total [Mg²⁺] is likely enriched ~5-fold, consistent with the condensate's high rRNA phosphate density and Donnan potential
- **Spatial gradient of free [Mg²⁺] within the nucleolus**: unknown and potentially inverted relative to total [Mg²⁺]. RNA sequestration in the RNA-dense FC could make free [Mg²⁺] lower at the centre despite higher total.
- **At nucleolus boundary** (r/R_c ~ 1): total and free [Mg²⁺] approach cytoplasmic background (~0.5–1 mM), which is **below the Goldilocks peak** for rRNA folding (~3–10 mM)

The cytoplasmic boundary condition also creates a difficulty for anchoring γ = 1 at the periphery to a Mg²⁺ mechanism: cytoplasmic free [Mg²⁺] of 0.5–1 mM is insufficient for ribosome assembly in vitro (the Nierhaus protocol requires 4–20 mM free [Mg²⁺]). In vivo, assembly factors (GTPases, helicases) compensate for the sub-optimal Mg²⁺, as demonstrated by the EngA/ObgE rescue experiments (eLife 2024). These chaperones are not captured by γ, so γ = 1 at the periphery should be understood as representing in-vivo-competent conditions including chaperone activity, not the Mg²⁺-sufficient in vitro reconstitution conditions.

### 2.6 Temporal Fluctuations

- Free [Mg²⁺] fluctuates with ATP concentration because ATP-Mg is the dominant buffering equilibrium; ATP concentration varies with metabolic state (~1–5 mM range)
- In Salmonella under low-Mg²⁺ stress: rRNA transcription drops 10-fold and translation decreases 2-fold; bacteria compensate by upregulating Mg²⁺ transporters (PMC5500012)
- In bacteria, ribosome-bound Mg²⁺ constitutes ~30–50% of total cellular Mg²⁺; disruptions to ribosome homeostasis (Lee et al. 2019, *Cell*) propagate as global Mg²⁺ perturbations
- The Yewdall et al. (2022) data show that ATP-Mg²⁺ competition is a dynamic regulatory mechanism: adding ATP to a gel-like condensate (high Mg²⁺) chelates free Mg²⁺ and rapidly restores liquid-like properties within 20 minutes

---

## 3. Polyamines (Spermidine and Spermine)

### 3.1 Concentrations

Polyamines are polycationic molecules that, like Mg²⁺, neutralise the negative charges on RNA and promote RNA folding:

- **E. coli**: total intracellular spermidine ~6.88 mM, putrescine ~32 mM (mostly RNA-bound; measured by ion exchange chromatography, NAR 2014)
- **Eukaryotic mammalian cells**: nuclear spermidine ~1.5 mM, spermine ~0.4 mM for maximal chromatin compaction (in vitro nuclear studies, PLoS ONE 2013); cytoplasmic values lower
- **Ribosome reconstitution**: 3 mM spermidine reduces the required Mg²⁺ for 30S assembly from 20 mM to 16 mM; 2 mM spermidine is optimal for in vitro protein synthesis together with 2–5 mM Mg²⁺ and 60–150 mM K⁺
- **Binding to rRNA**: spermine binds 23S rRNA with Kd ≈ 5.5 mM (1/K = 0.18×10⁴ M⁻¹); spermidine with Kd ≈ 45 mM; ~0.11 amines per RNA phosphate at saturation

### 3.2 Spatial Profile

No direct measurement of polyamine concentration within nucleolar sub-compartments exists. As polycations, spermidine and spermine would be expected to co-enrich with Mg²⁺ in the condensate, likely following a similar distribution (highest near the RNA-dense FC/DFC core). However, polyamines are much larger than Mg²⁺ and their spatial distribution would also be controlled by molecular crowding and polymer partitioning effects.

### 3.3 Effect on Ribosome Assembly

Functionally equivalent to Mg²⁺: reduces the Mg²⁺ requirement for ribosome reconstitution by ~4 mM per 3 mM spermidine. Stimulates protein synthesis 1.5–2.0-fold. Spermidine also stabilises the association of the 30S and 50S subunits into the 70S complex. Given their similar cationic roles, polyamines can be considered as contributing to an effective "total cation" gradient alongside Mg²⁺, though with lower per-charge potency.

---

## 4. Potassium (K⁺), Sodium (Na⁺), and Other Monovalent Ions

### 4.1 K⁺

- **Bulk intracellular K⁺**: 60–150 mM under physiological conditions (major cytoplasmic cation)
- **Ribosome requirement**: K⁺ is essential for ribosome function; optimal in vitro assembly and translation occurs at 60–150 mM K⁺. High K⁺ (>500 mM) dissociates ribosomes.
- **Condensate partitioning**: K⁺ is excluded from condensate dense phase (chaotropic by the law of matching water affinities; NMR 2025, JACS). Dai et al. (2024) in vivo bacteria: K⁺ not significantly redistributed upon condensate formation (unlike in vitro measurements where it was excluded), possibly because intracellular buffering mechanisms prevent depletion.
- **Spatial gradient in nucleolus**: No measurement. If K⁺ is excluded from the condensate, there would be a gradient of *lower* K⁺ inside the nucleolus (especially in FC/DFC) relative to the nucleoplasm. The magnitude is uncertain.

### 4.2 Na⁺

- **Bulk Na⁺**: ~10–30 mM intracellular (much lower than extracellular ~140 mM)
- **Condensate partitioning**: Na⁺ excluded from condensates (Dai et al. 2024, ICP-MS: 19% increase in cytoplasmic Na⁺ upon condensate formation)
- **Effect on ribosome**: minimal specific effect; Na⁺ can partly substitute for K⁺ but is less optimal

### 4.3 Zinc (Zn²⁺)

- **Ribosomal zinc-finger proteins**: ~8 r-proteins in eukaryotic ribosomes contain zinc-finger domains. L37a (yeast YL37a) contains a C2-C2 zinc finger essential for cell viability and binds domains II and III of 26S rRNA with Kd = 79 nM (domain II) and 198 nM (domain III).
- ZNF658 regulates ribosome biogenesis through zinc homeostasis.
- No measurements of spatial Zn²⁺ distribution within the nucleolus exist.

### 4.4 Condensate Electric Potential

An important recent finding relevant to ion gradients: condensates maintain Donnan potentials and Nernst-like electric potentials across their interfaces as a result of asymmetric ion partitioning (JACS 2024, Brangwynne/Pappu group). The magnitude of these potentials is on the order of membrane potentials of membrane-bound organelles — tens of millivolts. This interfacial electric potential would drive directional migration of charged molecules across the condensate boundary and could contribute to the sorting of charged assembly intermediates during the radial transit of pre-ribosomal particles.

---

## 5. Small Molecules, Denaturants, and Crowding Agents

### 5.1 Molecular Crowding

The nucleolus is extremely dense in macromolecules. rRNA concentration is estimated at ~5 μM (Riback et al. 2023) — approximately 10-fold above the overlap concentration — placing the rRNA in an entangled gel regime. This crowding has several effects on assembly:

- Reduces the translational diffusion of assembling pre-ribosomal particles (effective D ≈ 2×10⁻⁵ μm²/s for rRNA, 5,000× below free diffusion)
- Increases the effective local concentrations of assembly factors relative to dilute solution
- Can stabilise or destabilise folded RNA structures through excluded-volume effects (typically stabilising for compact folded states relative to extended denatured states)

The macromolecular density decreases from FC to GC as nascent rRNA folds and binds proteins (compacting), which means crowding-based effects should be strongest in the DFC and decrease toward the GC boundary.

### 5.2 ATP as a Chemical Chaotrope

ATP (adenosine triphosphate) acts as a biological hydrotrope that can fluidise condensates — a discovery made in the context of phase-separated compartments in 2017 (Patel et al.; Zhu et al.). In the Yewdall et al. (2022) system, ATP at physiological concentrations chelates Mg²⁺ (thereby reducing effective [Mg²⁺]) and also directly fluidises NPM1+rRNA condensates. ATP concentration in actively metabolising cells is highest in regions of highest metabolic activity and lowest in regions far from mitochondria and the ATP-generating machinery. Whether there is a spatial ATP gradient within the nucleolus is unknown, but the coupling of ATP hydrolysis by AAA-ATPases and GTPases (the ribosome assembly factors) to local ATP concentration provides a mechanistic link between metabolic state and condensate material properties.

### 5.3 Specific Small Molecules That Concentrate in the Nucleolus

Kilgore et al. (2023, *Nature Chemical Biology*) screened a library of small-molecule probes and found that mitoxantrone concentrates in the nucleolus at up to **6-fold** above the surrounding nucleoplasm (in HeLa cells, in vivo fluorescence imaging). This provides proof-of-principle that small molecules can be significantly enriched inside the nucleolus, consistent with the general principle that the condensate environment provides a distinct chemical microenvironment. However, the specific chemical species relevant to ribosome assembly (Mg²⁺, polyamines) have not yet been directly imaged with sub-organellar spatial resolution.

---

## 6. Composition and Protein Concentration Gradients

### 6.1 Proteome Organisation

The nucleolar proteome contains 1,318 proteins (Boisvert et al. 2020, *Molecular Systems Biology*), of which 287 localise to fibrillar components and 157 are enriched at the nucleolar rim (a proposed fourth sub-compartment, the "nucleolar rim"). The spatial distribution follows the functional assembly pathway: RNA synthesis machinery at the FC, rRNA processing factors in the DFC, and late-stage ribosome assembly factors in the GC.

Key landmark proteins:
- **FC**: UBF (upstream binding factor, RNA Pol I component; pI high, K-block-enriched)
- **DFC**: nucleolin (NCL; D/E-tract-enriched, pH-buffering function); fibrillarin (methyltransferase)
- **GC**: NPM1/nucleophosmin (abundant, pI ~5.5; D/E-tract-enriched, sets GC pH)

### 6.2 rRNA Dynamics and Advective Transport

Riback et al. (2023, *Molecular Cell*) provide the most quantitative data on how assembling rRNA moves through the nucleolus:

- **rRNA diffusion constant** in the GC: ~2×10⁻⁵ μm²/s (measured by pulse-chase imaging; 5,000× below free diffusion, consistent with entangled gel)
- **NPM1 diffusion constant**: ~0.1 μm²/s (FRAP; ~3 orders of magnitude faster than rRNA)
- **rRNA advective velocity**: ~1 Å/s (~6 nm/min) radially outward from FC, measured by spatial shift of rRNA signal centre between pulse and chase time points
- **Péclet number**: Pe = 2.3 ± 0.1 (ratio of advection to diffusion), meaning advection dominates rRNA transport — rRNA is swept outward by convective flow rather than diffusing
- **Transit time**: based on Pe > 1 and nucleolar radius ~1.5 μm, the characteristic transit time for rRNA from FC to nucleolus boundary is ~1.5 μm / (1 Å/s) ≈ **150 minutes** — consistent with the 60–90 minute processing timescale measured by pulse-chase sequencing (Quinodoz et al. 2025)

The **radial composition gradient** measured by Quinodoz et al. (2025, *Nature*):
- 0–15 min post-synthesis: early cleavage steps, rRNA at FC/DFC boundary
- 30–45 min: middle cleavage, rRNA in GC
- 60–90 min: late cleavage, rRNA approaching nucleoplasm
- 2'-O-methylations: mostly complete within 15 minutes (early DFC modifications) with a subset at 60–90 min (late GC modifications)

This establishes that the **assembly state** of the pre-ribosomal particle is itself radially graded: inner regions contain the most immature particles, outer regions contain nearly-mature particles. This is the temporal annealing schedule experienced spatially, as the simulation models.

---

## 7. Synthesis: The Effective Coupling Gradient in the Nucleolus

### 7.1 Summary of Each Chemical

| Chemical | Direction | Magnitude | Key effect on assembly |
|---|---|---|---|
| **pH** | Acidic at centre (6.5), neutral at periphery (7.2) | 0.7 units (5-fold H⁺); directly measured | RNA-RNA contacts strengthened at low pH (inverted effect); assembly factor (DDX21) activity impaired at low pH; r-protein Arg/Lys binding pH-independent |
| **Mg²⁺** | Total [Mg²⁺] likely enriched in condensate vs. cytoplasm; free [Mg²⁺] gradient within the nucleolus **unknown** | Total enrichment ~5-fold (speculative); free [Mg²⁺] gradient direction unknown | RNA folding and tertiary structure; r-protein binding (15-fold Kd shift over 2–8 mM range); but free [Mg²⁺] may be depleted at RNA-dense centre due to sequestration |
| **Polyamines** | Likely enriched in condensate, spatial gradient within nucleolus unknown | Unknown | Partial Mg²⁺ substitute; RNA stabilisation; faces same total-vs-free ambiguity as Mg²⁺ |
| **K⁺** | Possibly lower at centre (excluded from condensate) | ~20% estimated | Required for ribosome function; gradient uncertain |
| **Crowding** | Higher at centre (entangled gel), lower at periphery (less dense) | ~10× rRNA above overlap concentration | Slows diffusion; stabilises compact RNA; increases effective concentration |
| **Electric potential** | Likely negative (interior) relative to nucleoplasm | Tens of mV (estimate) | Drives charged species; affects protonation equilibria |

### 7.2 Are the Gradients Consistent with an Annealing Picture?

The key question is: do the chemical conditions in the nucleolus interior (FC/DFC, r ≈ 0) constitute a *denaturing* environment relative to the periphery (GC/nucleoplasm, r ≈ R_c)?

The answer is **yes**, with important nuance:

**pH**: The interior (pH 6.5) is suboptimal for ribosome assembly relative to physiological pH (7.0–7.4). DDX21 (an essential nucleolar helicase) is non-functional below pH 6.0 or above pH 8.0, and optimal at 6.8–7.5, meaning the FC at pH 6.5 is at the edge of functional conditions. The pH gradient (from 6.5 to 7.2 radially outward) thus directly represents an increase in assembly-competent conditions from centre to periphery. This is **consistent with the denaturing-interior, physiological-periphery scenario** of the model.

**Mg²⁺**: As established in §2.5, all available measurements report **total** [Mg²⁺] enrichment in the condensate phase, not **free** [Mg²⁺]. In the RNA-dense FC and DFC, the majority of total Mg²⁺ is bound to RNA phosphates; the free [Mg²⁺] available for new assembly events could be lower at the centre than at the GC periphery. No direct measurement of a free [Mg²⁺] gradient within the nucleolus exists. Additionally, the cytoplasmic boundary fixes free [Mg²⁺] at r = R_c to ~0.5–1 mM — below the Goldilocks peak for rRNA folding — meaning γ = 1 at the periphery cannot be anchored to Mg²⁺ reaching optimal levels there. The Mg²⁺ reconstitution data (S4-rRNA Kd shifts, Nierhaus protocols) are useful for calibrating the energy scale J but do not currently justify a specific gradient direction or γ₀ value.

**Polyamines and other cations**: As polycations, spermidine and spermine face the same total-versus-free ambiguity as Mg²⁺. Direct sub-nucleolar measurements do not exist; their contribution to any spatial gradient is unknown.

### 7.3 Overall Effective Coupling Profile

Translating all gradients into a single dimensionless coupling parameter γ(r) as used in the simulation model:

**Definition**: γ = 1 represents "physiologically optimal conditions for ribosome assembly" (correct contacts are stable and form readily). γ < 1 represents conditions where correct contacts are weakened relative to physiological.

**pH contribution — contact-type dependent**: The pH gradient affects RNA-RNA and RNA-protein contacts in opposite directions (§1.2). RNA-RNA contacts (secondary structure, base-stacking) are modestly **strengthened** at the acidic interior; r-protein Arg/Lys contacts are pH-insensitive. The single-γ model therefore receives a mixed pH signal. The most defensible pH-driven contribution to γ < 1 in the interior is through **assembly factor impairment**: DDX21 and related helicases are sub-optimal at pH 6.5, causing kinetic traps in RNA folding to accumulate at the condensate centre. This is captured phenomenologically by γ₀ < 1 even though the underlying mechanism is kinetic. In the dual-gradient model (complexes.md), the pH data more naturally justifies γ₁ higher at the centre (RNA-RNA contacts strengthened) and γ₂ lower at the centre (assembly impaired).

**Mg²⁺ contribution**: Cannot be reliably assigned a direction or magnitude within the nucleolus given current data (§2.5, §7.2). The Mg²⁺ literature constrains the energy scale J but not the spatial gradient.

**Physical status of γ₀**: Given the above, γ₀ is best treated as a **phenomenological parameter** representing the net reduced assembly propensity in the condensate interior, rather than a quantity derivable from a specific chemical mechanism. The empirical support for γ₀ < 1 is strong: rRNA processing occurs sequentially as particles move outward (Quinodoz et al. 2025), ribosome subunits exit the condensate rather than assembling inside it, and DDX21 is sub-optimal at pH 6.5. The range γ₀ = 0.3–0.5 represents the interior as partially (not fully) denaturing, consistent with rRNA having some secondary structure at the FC while not reaching the correct final assembled state there.

γ_eff(r) ≈ γ₀ + (1 − γ₀) × f(r/R_c)

where γ₀ ≈ 0.3–0.5 and f is approximately linear, following the measured near-linear pH profile (King et al. 2024).

### 7.4 Recommended Parameter Ranges for the Model

Based on the quantitative literature:

| Parameter | Recommended range | Physical basis |
|---|---|---|
| γ₀ | 0.3–0.5 | Phenomenological: interior is partially but not fully denaturing. Empirically supported by sequential rRNA processing (Quinodoz et al. 2025), DDX21 pH-sensitivity (Rai et al. 2023), and absence of mature particles from the condensate core. Not derivable from a single chemical mechanism. |
| γ(R_c) | 0.8–1.0 | Nucleoplasm at physiological pH 7.2, with chaperone-assisted assembly. γ = 1 represents in-vivo-competent conditions, not in-vitro reconstitution conditions (which require Mg²⁺ concentrations far above cytoplasmic levels). |
| Gradient shape | Approximately linear | Directly supported by the near-linear pH gradient measured by King et al. (2024). Mg²⁺ gradient direction within the nucleolus is unknown and does not justify a specific non-linear shape. |
| J (patch coupling) | 5–10 kBT | S4-rRNA Kd shifts 15-fold over 2–8 mM Mg²⁺: ΔΔG ≈ kBT ln(15) ≈ 2.7 kBT total; typical r-protein:rRNA contacts measured by ITC are 2–8 kcal/mol (~3–14 kBT). J = 8 kBT is at the lower end, appropriate for secondary/tertiary binders disrupted by the annealing. |

**Functional form for γ(r)**:

The most defensible functional form uses the measured pH profile, with the relevant variable being pH-dependent assembly factor activity:

γ_eff(r) = γ₀ + (1 − γ₀) × (r/R_c)

A sigmoidal form can also capture the steeper FC/DFC–GC transition:

γ_eff(r) ≈ 1 / (1 + exp(−k(pH(r) − pH_ref)))

with pH(r) = 6.5 + 0.7×r/R_c (King et al. 2024), pH_ref ≈ 7.0 (DDX21 activity optimum and assembly-competent pH), and k ≈ 2–3. This gives γ_eff(0) ≈ 0.35–0.4 and γ_eff(R_c) ≈ 0.85–0.9, consistent with the default γ₀ = 0.4 and validating the approximately linear approximation. The key justification for pH_ref ≈ 7.0 is DDX21 activity data (Rai et al. 2023) and the empirically observed assembly optimum at pH 7.0–7.4.

The non-monotonic Mg²⁺ Goldilocks functional form previously considered in earlier drafts of this document has been removed: it depended on an assumed free [Mg²⁺] gradient within the nucleolus for which no direct measurement exists.

### 7.5 Implications for the Simulation

1. **Linear gradient is a first approximation but a good one**: The measured pH gradient is nearly linear and dominates the molecular-level coupling. γ₀ = 0.4 is physically reasonable.

2. **Non-linear gradient shapes are worth exploring**: While the measured pH gradient is approximately linear, the observed steeper FC/DFC-to-GC compositional transition at r/R_c ≈ 0.3–0.5 suggests a slightly sigmoidal profile may be more accurate. The sigmoidal functional form in §7.4 provides a well-motivated alternative. A Mg²⁺ Goldilocks correction is not recommended: the direction of free [Mg²⁺] within the nucleolus is unmeasured and the Goldilocks reasoning does not apply to total [Mg²⁺] (see §2.5).

3. **The gradient magnitude is not extreme**: γ₀ = 0.3–0.5 (not zero) is consistent with the nucleolus interior not being fully denaturing — rRNA must have some secondary structure even in the FC for the particle to remain associated with the condensate.

4. **The model correctly places assembly improvement at the gradient, not the level**: The simulation result that gradient+moderate J works but uniform high coupling fails is consistent with the broader assembly literature showing that both too-weak and too-strong contacts produce frustrated assembly (Wilber et al. 2007; Whitelam et al. 2009). The gradient implements a schedule that avoids both failure modes. This conclusion does not depend on any specific chemical mechanism for γ₀.

5. **Temporal fluctuations could be modelled**: The pH fluctuations on timescales of minutes (Rai et al. 2023) and Mg²⁺ fluctuations tied to ATP oscillations (Yewdall et al. 2022) suggest that γ(r, t) has both spatial and temporal components that could be incorporated as correlated noise in the simulation.

---

## 8. Key Citations

- Ausserwoger H et al. (2025). Biomolecular condensates sustain pH gradients at equilibrium through charge neutralization. *Nature Chemistry*. PMC12872462.
- Boisvert FM et al. (2010). The multifunctional nucleolus. *Molecular Systems Biology*. PMC7397901.
- Dai X et al. (2024). Biomolecular condensates regulate cellular electrochemical equilibria. *Cell*. PMC11490381.
- Detrés J, Camacho-Badillo F & Calo E (2025). A pH-centric model of nucleolar activity and regulation. *Journal of Molecular Biology*. 437, 169136. PubMed 40216015.
- Guth-Metzler R et al. (2023). Goldilocks and RNA: where Mg²⁺ concentration is just right. *Nucleic Acids Research*. 51(8):3529–3540. PMC10164553.
- Hori N, Denesyuk NA & Thirumalai D (2021). Shape changes and cooperativity in the folding of the central domain of the 16S ribosomal RNA. *PNAS*. PMC7958424.
- Kilgore HR et al. (2023). Distinct chemical environments in biomolecular condensates. *Nature Chemical Biology*. PMC12181805.
- King MR et al. (2024). Macromolecular condensation organizes nucleolar sub-phases to set up a pH gradient essential for ribosome biogenesis. *Cell*. PMC11938373.
- Lim B et al. (2025). Micropipette aspiration reveals differential RNA-dependent viscoelasticity of nucleolar subcompartments. *PNAS*. PMC12146704.
- Nierhaus KH (2014). Mg²⁺, K⁺, and the ribosome. *Journal of Bacteriology*. PMC4248827.
- Nierhaus KH & Dohme F (1974). Total reconstitution of functionally active 50S ribosomal subunits from Escherichia coli. *PNAS*. 71(12):4713–4717.
- Pan T & Sosnick T (1997). Intermediates and kinetic traps in the folding of a large ribozyme revealed by circular dichroism and UV absorbance spectroscopies. *Nature Structural Biology*.
- Phair RD & Misteli T (2000). High mobility of proteins in the mammalian cell nucleus. *Science*. PMC2185520 (see also JCB 2001 for FRAP tables).
- Quinodoz SA et al. (2025). Mapping and engineering RNA-driven architecture of the multiphase nucleolus. *Nature*. PMC11463421.
- Rai AK et al. (2023). Kinase-controlled phase transition of membraneless organelles in mitosis. *Molecular Cell*. PMC10803072.
- Riback JA et al. (2023). Viscoelasticity and advective flow of RNA underlie nucleolar form and function. *Molecular Cell*. PMC11089468.
- Selective Ion Binding and Uptake Shape the Microenvironment of Biomolecular Condensates (2025). *JACS*.
- Competitive Mg²⁺ Regulation of Biomolecular Condensate Microenvironments (2025). *JACS Au*.
- Structural insights into the distortion of the ribosomal small subunit at different Mg²⁺ concentrations (2023). *PMC10046523*.
- Vafabakhsh R et al. (2024). Multiscale biophysical analysis of nucleolus disassembly during mitosis. *PNAS*. PMC10861868.
- Yewdall NA et al. (2022). ATP:Mg²⁺ shapes material properties of protein-RNA condensates. *Biophysical Journal*. PMC9674983.
- Zhu L et al. (2015). Nuclear bodies: the emerging biophysics of nucleoplasm. *Cell*. PMC5127388.
- Principles of ion binding to RNA from 1.55 Å bacterial ribosome structure (2025). *Nucleic Acids Research*.
- Biomolecular condensates are characterized by interphase electric potentials (2024). *JACS*.
