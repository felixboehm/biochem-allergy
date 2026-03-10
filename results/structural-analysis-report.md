# Structural Analysis of Nutrient-Dependent Allergy Enzymes

## Overview

This report presents a structural and biochemical analysis of seven key
enzymes involved in allergy biochemistry, each dependent on specific
micronutrients (metals, vitamins, trace elements) for their catalytic activity.

The central thesis: **micronutrient deficiencies can systematically disable
the enzymatic machinery that regulates histamine, leukotrienes, tryptophan
metabolism, and oxidative stress — creating a biochemical environment that
amplifies allergic responses.**

---

## 1. Protein Structures Retrieved

| Enzyme | Source | PDB ID | Resolution | Organism |
|--------|--------|--------|------------|----------|
| DAO | RCSB PDB | 3HI7 | 1.80 Å | Homolog/Predicted |
| HNMT | RCSB PDB | 1JQD | 2.28 Å | Human |
| HDC | RCSB PDB | 4E1O | 1.80 Å | Homolog/Predicted |
| IDO1 | RCSB PDB | 2D0T | 2.30 Å | Human |
| ALOX5 | RCSB PDB | 3O8Y | 2.39 Å | Human |
| SOD1 | RCSB PDB | 1PU0 | 1.70 Å | Human |
| GPX1 | RCSB PDB | 2F8A | 1.50 Å | Human |

---

## 2. Cofactor Binding Site Analysis

### DAO — Diamine Oxidase (AOC1/ABP1)

- **Gene:** AOC1
- **Cofactors:** Copper (Cu²⁺), Topaquinone (TPQ, from Tyr + Cu)
- **Required Nutrients:** Copper, Vitamin B6 (PLP)
- **Substrate:** Histamine
- **Reaction:** Histamine → Imidazole acetaldehyde + NH₃ + H₂O₂
- **Function:** Primary extracellular histamine degradation

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| ASN | 460 | B | C | N (TPQ) | 1.33 |
| ASN | 460 | A | C | N (TPQ) | 1.34 |
| HIS | 512 | B | NE2 | CU (CU) | 2.04 |
| HIS | 512 | A | NE2 | CU (CU) | 2.08 |
| ASP | 462 | A | N | C (TPQ) | 2.41 |
| ASP | 462 | B | N | C (TPQ) | 2.42 |
| THR | 685 | A | OG1 | O (TPQ) | 2.69 |
| THR | 685 | B | OG1 | O (TPQ) | 2.69 |
| HIS | 675 | B | CE1 | CU (CU) | 2.97 |
| HIS | 510 | B | CE1 | CU (CU) | 2.99 |

**When Copper, Vitamin B6 (PLP) is/are deficient:**
> Histamine accumulates → pseudo-allergic reactions, food intolerance

**Clinical consequence:** Histamine intolerance, migraine, GI symptoms, skin flushing

---

### HNMT — Histamine N-Methyltransferase

- **Gene:** HNMT
- **Cofactors:** S-adenosylmethionine (SAM)
- **Required Nutrients:** Methionine, Folate (B9), B12, B6
- **Substrate:** Histamine
- **Reaction:** Histamine + SAM → Nτ-methylhistamine + SAH
- **Function:** Intracellular histamine degradation (brain, airways)

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| GLU | 89 | A | OE1 | O (SAH) | 2.58 |
| GLU | 28 | B | OE1 | N (HSM) | 2.58 |
| ASN | 283 | A | OD1 | N (HSM) | 2.64 |
| GLU | 28 | A | OE1 | N (HSM) | 2.69 |
| GLY | 60 | A | O | N (SAH) | 2.81 |
| ILE | 142 | B | O | N (SAH) | 2.87 |
| ASN | 283 | B | OD1 | N (HSM) | 2.87 |
| SER | 120 | A | N | N (SAH) | 2.88 |
| GLN | 143 | B | NE2 | N (HSM) | 2.90 |
| ILE | 142 | A | O | N (SAH) | 2.92 |

**When Methionine, Folate (B9), B12, B6 is/are deficient:**
> Histamine accumulates intracellularly → neurological & respiratory symptoms

**Clinical consequence:** Brain fog, insomnia, anxiety, asthma exacerbation

---

### HDC — Histidine Decarboxylase

- **Gene:** HDC
- **Cofactors:** Pyridoxal 5'-phosphate (PLP, vitamin B6)
- **Required Nutrients:** Vitamin B6 (PLP)
- **Substrate:** L-Histidine
- **Reaction:** L-Histidine → Histamine + CO₂
- **Function:** Histamine biosynthesis

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| SER | 151 | E | OG | O (PLP) | 2.46 |
| SER | 151 | F | OG | O (PLP) | 2.46 |
| SER | 354 | D | OG | O (PLP) | 2.47 |
| SER | 151 | C | OG | O (PLP) | 2.48 |
| SER | 151 | B | OG | O (PLP) | 2.49 |
| SER | 354 | C | OG | O (PLP) | 2.49 |
| SER | 354 | A | OG | O (PLP) | 2.50 |
| SER | 151 | A | OG | O (PLP) | 2.51 |
| SER | 354 | E | OG | O (PLP) | 2.51 |
| SER | 151 | D | OG | O (PLP) | 2.52 |

**When Vitamin B6 (PLP) is/are deficient:**
> With adequate B6: normal histamine production. B6 deficiency paradoxically may increase mast cell histamine release

**Clinical consequence:** Dysregulated histamine production, mast cell instability

---

### IDO1 — Indoleamine 2,3-Dioxygenase 1

- **Gene:** IDO1
- **Cofactors:** Heme/Protoporphyrin IX (contains Fe²⁺)
- **Required Nutrients:** Iron, Heme
- **Substrate:** L-Tryptophan
- **Reaction:** L-Trp + O₂ → N-formylkynurenine
- **Function:** Tryptophan catabolism, immune regulation via kynurenine pathway

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| HIS | 346 | A | NE2 | FE (HEM) | 2.09 |
| HIS | 346 | B | NE2 | FE (HEM) | 2.11 |
| SER | 263 | A | OG | O (HEM) | 2.55 |
| ARG | 343 | B | NH2 | O (HEM) | 2.58 |
| SER | 263 | B | OG | O (HEM) | 2.59 |
| LEU | 384 | B | CD1 | C (HEM) | 3.18 |
| LEU | 384 | A | CD1 | C (HEM) | 3.19 |
| ARG | 343 | A | NE | C (HEM) | 3.30 |
| VAL | 391 | A | CG1 | O (HEM) | 3.33 |
| VAL | 391 | B | CG1 | O (HEM) | 3.34 |

**When Iron, Heme is/are deficient:**
> Iron deficiency → reduced IDO1 → altered Trp metabolism → immune dysregulation

**Clinical consequence:** Immune imbalance, Th1/Th2 skewing, loss of oral tolerance

---

### ALOX5 — Arachidonate 5-Lipoxygenase

- **Gene:** ALOX5
- **Cofactors:** Non-heme Iron (Fe²⁺/Fe³⁺)
- **Required Nutrients:** Iron
- **Substrate:** Arachidonic acid
- **Reaction:** Arachidonic acid → 5-HPETE → LTA₄ → Leukotrienes
- **Function:** Produces leukotrienes (potent inflammatory mediators in asthma/allergy)

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| HIS | 372 | A | NE2 | FE (FE2) | 2.10 |
| HIS | 550 | B | NE2 | FE (FE2) | 2.14 |
| HIS | 372 | B | NE2 | FE (FE2) | 2.16 |
| HIS | 550 | A | NE2 | FE (FE2) | 2.21 |
| ILE | 673 | B | O | FE (FE2) | 2.21 |
| ILE | 673 | A | O | FE (FE2) | 2.24 |
| HIS | 367 | A | NE2 | FE (FE2) | 2.24 |
| HIS | 367 | B | NE2 | FE (FE2) | 2.27 |
| ASN | 554 | B | OD1 | FE (FE2) | 3.18 |
| ASN | 554 | A | OD1 | FE (FE2) | 3.22 |

**When Iron is/are deficient:**
> Iron-dependent: both excess and deficiency alter leukotriene production

**Clinical consequence:** Asthma, bronchoconstriction, allergic inflammation

---

### SOD1 — Superoxide Dismutase 1 (Cu/Zn)

- **Gene:** SOD1
- **Cofactors:** Copper (Cu²⁺, catalytic), Zinc (Zn²⁺, structural)
- **Required Nutrients:** Copper, Zinc
- **Substrate:** Superoxide radical (O₂⁻)
- **Reaction:** 2 O₂⁻ + 2H⁺ → O₂ + H₂O₂
- **Function:** Antioxidant defense — removes superoxide radicals

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| ASP | 83 | G | OD1 | ZN (ZN) | 1.88 |
| ASP | 83 | D | OD1 | ZN (ZN) | 1.90 |
| ASP | 83 | I | OD1 | ZN (ZN) | 1.91 |
| ASP | 83 | H | OD1 | ZN (ZN) | 1.92 |
| ASP | 83 | B | OD1 | ZN (ZN) | 1.93 |
| HIS | 80 | C | ND1 | ZN (ZN) | 1.95 |
| HIS | 46 | G | ND1 | CU (CU1) | 1.95 |
| HIS | 80 | G | ND1 | ZN (ZN) | 1.95 |
| ASP | 83 | J | OD1 | ZN (ZN) | 1.95 |
| HIS | 120 | C | NE2 | CU (CU1) | 1.96 |

**When Copper, Zinc is/are deficient:**
> Cu/Zn deficiency → oxidative stress → mast cell activation, tissue damage

**Clinical consequence:** Oxidative damage amplifies allergic inflammation, epithelial barrier breakdown

---

### GPX1 — Glutathione Peroxidase 1

- **Gene:** GPX1
- **Cofactors:** Selenocysteine (Sec, amino acid with Se)
- **Required Nutrients:** Selenium
- **Substrate:** H₂O₂ / Lipid hydroperoxides
- **Reaction:** 2 GSH + H₂O₂ → GSSG + 2 H₂O
- **Function:** Reduces peroxides, prevents oxidative damage

**Coordinating Residues (closest to cofactor):**

| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |
|---------|----------|-------|------|----------|---------------|
| CYS | 49 | A | SG | Se (SEC) | 0.00 |
| GLN | 80 | A | NE2 | Se (SEC) | 3.20 |
| TRP | 160 | A | NE1 | Se (SEC) | 3.40 |

**When Selenium is/are deficient:**
> Se deficiency → oxidative stress → NF-κB activation → pro-inflammatory cascade

**Clinical consequence:** Amplified allergic inflammation, asthma severity, oxidative tissue damage

---

## 3. Comparative Summary

### Nutrient-Enzyme-Consequence Map

| Nutrient | Enzymes Affected | Pathway Disrupted | Clinical Impact |
|----------|-----------------|-------------------|-----------------|
| B12 | HNMT | Histamine | Reduced SAM production → HNMT dysfunction |
| B6 | HNMT | Histamine |  |
| Copper | DAO, SOD1 | Antioxidant, Histamine | Histamine accumulation (DAO), oxidative stress (SOD1) |
| Folate (B9) | HNMT | Histamine | Reduced SAM production → HNMT dysfunction, epigenetic changes |
| Heme | IDO1 | Tryptophan/Kynurenine | IDO1 dysfunction → tryptophan metabolism disruption |
| Iron | IDO1, ALOX5 | Tryptophan/Kynurenine, Leukotriene | Altered leukotriene production, immune dysregulation, impaired tryptophan metabolism |
| Methionine | HNMT | Histamine | Reduced SAM → impaired histamine methylation (HNMT) |
| Selenium | GPX1 | Antioxidant | Oxidative stress, amplified NF-kB signaling, increased allergy severity |
| Vitamin B6 (PLP) | DAO, HDC | Histamine | Histamine metabolism impairment (DAO, HDC), neurotransmitter imbalance |
| Zinc | SOD1 | Antioxidant | Oxidative stress (SOD1 structural instability), immune dysfunction |

### The Vicious Cycle of Nutrient Depletion in Allergy

1. **Chronic inflammation** increases demand for antioxidant nutrients (Se, Zn, Cu)
2. **Gut inflammation** reduces nutrient absorption (Fe, B6, B12, folate)
3. **Histamine excess** (from impaired DAO/HNMT) further damages gut lining
4. **Oxidative stress** (from impaired SOD1/GPX1) activates mast cells
5. **Mast cell degranulation** releases more histamine → back to step 1

This creates a self-reinforcing cycle where nutrient deficiencies
worsen the very allergic responses that deplete those nutrients.

## 4. Figures

### Structural Figures
- `figures/structure_3d_DAO.png` — 3D Cα trace of DAO
- `figures/structure_3d_HNMT.png` — 3D Cα trace of HNMT
- `figures/structure_3d_HDC.png` — 3D Cα trace of HDC
- `figures/structure_3d_IDO1.png` — 3D Cα trace of IDO1
- `figures/structure_3d_ALOX5.png` — 3D Cα trace of ALOX5
- `figures/structure_3d_SOD1.png` — 3D Cα trace of SOD1
- `figures/structure_3d_GPX1.png` — 3D Cα trace of GPX1

### Schematic Figures
- `figures/schematic_DAO.png` — Domain architecture and deficiency impact for DAO
- `figures/schematic_HNMT.png` — Domain architecture and deficiency impact for HNMT
- `figures/schematic_HDC.png` — Domain architecture and deficiency impact for HDC
- `figures/schematic_IDO1.png` — Domain architecture and deficiency impact for IDO1
- `figures/schematic_ALOX5.png` — Domain architecture and deficiency impact for ALOX5
- `figures/schematic_SOD1.png` — Domain architecture and deficiency impact for SOD1
- `figures/schematic_GPX1.png` — Domain architecture and deficiency impact for GPX1

### Binding Site Figures
- `figures/binding_site_DAO.png` — Cofactor coordination diagram for DAO
- `figures/binding_site_HNMT.png` — Cofactor coordination diagram for HNMT
- `figures/binding_site_HDC.png` — Cofactor coordination diagram for HDC
- `figures/binding_site_IDO1.png` — Cofactor coordination diagram for IDO1
- `figures/binding_site_ALOX5.png` — Cofactor coordination diagram for ALOX5
- `figures/binding_site_SOD1.png` — Cofactor coordination diagram for SOD1
- `figures/binding_site_GPX1.png` — Cofactor coordination diagram for GPX1

### Summary Figures
- `figures/comparative_summary.png` — All 7 enzymes compared side by side
- `figures/nutrient_enzyme_network.png` — Nutrient-enzyme-consequence network

## 5. Methods

- **Structure retrieval:** RCSB PDB REST API and AlphaFold EBI API (v4)
- **Structure parsing:** BioPython PDB/mmCIF parsers
- **Binding site analysis:** NeighborSearch with 3.5 Å distance cutoff
- **Visualization:** matplotlib (3D Cα traces, 2D schematics)
- **All figures:** 300 DPI, white background
