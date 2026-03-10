# Allergy Biochemistry: Inflammation, Nutrition, and Medication Effects

A comprehensive biochemical reference covering IgE-mediated hypersensitivity, chronic allergic inflammation, nutrient status impacts, and the pharmacological and metabolic consequences of common allergy medications.

---

## 1. Allergy Biochemistry Fundamentals

### 1.1 IgE-Mediated Type I Hypersensitivity: The Sensitization Cascade

The allergic response begins with **sensitization**, a multi-step immunological priming:

1. **Allergen uptake by APCs**: Dendritic cells in mucosal/epithelial tissue capture allergen proteins, process them, and present peptide fragments on MHC class II molecules.
2. **Th2 polarization**: Naive CD4+ T cells encountering allergen-loaded DCs in the presence of IL-4 (from basophils, NKT cells, or other innate sources) differentiate into Th2 effector cells. GATA-3 is the master transcription factor driving Th2 commitment.
3. **Cytokine-driven B cell class switching**: Th2 cells secrete **IL-4 and IL-13**, which act on B cells to induce class-switch recombination (CSR) at the immunoglobulin heavy chain locus, switching from IgM/IgG to **IgE** production. This requires:
   - IL-4/IL-13 signaling through IL-4Ralpha/STAT6
   - CD40L-CD40 co-stimulation (T cell-B cell contact)
   - Activation-induced cytidine deaminase (AID) for CSR at the Smu-to-Sepsilon switch regions
4. **IgE binding to mast cells**: Secreted IgE binds with extraordinarily high affinity (Ka ~10^10 M^-1) to **FcepsilonRI** (the high-affinity IgE receptor) on mast cell and basophil surfaces. FcepsilonRI is a tetrameric receptor (alphabetagamma2) where the alpha chain binds IgE, the beta chain amplifies signaling, and the gamma chain dimer contains ITAMs essential for signal transduction.
5. **Primed state**: Mast cells now sit loaded with allergen-specific IgE, awaiting re-exposure. This primed state can persist for weeks to months given the long half-life of IgE on mast cells (vs. ~2 days free in serum).

### 1.2 Mast Cell Degranulation: The Signaling Cascade

Upon re-exposure, multivalent allergen crosslinks two or more IgE-FcepsilonRI complexes, triggering a precisely orchestrated intracellular signaling cascade:

1. **Receptor aggregation and Lyn activation**: Crosslinking brings FcepsilonRI complexes into proximity. The Src-family kinase **Lyn**, constitutively associated with the beta chain, is activated (assisted by CD45 phosphatase activity). Lyn phosphorylates **ITAMs** (immunoreceptor tyrosine-based activation motifs) on both the beta and gamma chains.
2. **Syk recruitment and activation**: Phosphorylated gamma-chain ITAMs recruit **Syk** kinase via its tandem SH2 domains. Lyn then phosphorylates and fully activates Syk. Syk is the critical amplifying kinase of the cascade.
3. **Scaffold assembly (LAT/SLP-76)**: Activated Syk phosphorylates the transmembrane adaptor **LAT** (linker for activation of T cells) and cytosolic adaptor **SLP-76**. These serve as scaffolds recruiting:
   - PLCgamma1 and PLCgamma2
   - Grb2, Gads, SOS
   - Vav1 (guanine nucleotide exchange factor)
   - PI3K (phosphoinositide 3-kinase)
4. **PLCgamma activation and second messengers**: **PLCgamma** cleaves membrane PIP2 (phosphatidylinositol 4,5-bisphosphate) into two critical second messengers:
   - **IP3** (inositol 1,4,5-trisphosphate) --> binds IP3 receptors on the ER --> releases intracellular Ca2+ stores --> triggers store-operated calcium entry (SOCE) via STIM1/Orai1 channels --> sustained Ca2+ elevation
   - **DAG** (diacylglycerol) --> activates **PKC** (protein kinase C) isoforms --> contributes to granule fusion machinery
5. **Calcium flux and degranulation**: The sustained Ca2+ rise drives:
   - SNARE complex assembly (VAMP/synaptobrevin, SNAP-23, syntaxin-4) for granule-plasma membrane fusion
   - Cytoskeletal reorganization (actin disassembly at cortex to allow granule access)
   - Compound exocytosis: granule-to-granule fusion enabling rapid bulk release
6. **Parallel signaling arms**:
   - **PI3K pathway**: generates PIP3 --> recruits/activates Btk, PLCgamma, and Akt --> amplifies calcium signaling and promotes survival
   - **MAPK cascades**: Ras-Raf-MEK-ERK and JNK/p38 pathways --> activate transcription factors (Elk-1, c-Fos, c-Jun) for de novo mediator synthesis
   - **NF-kappaB activation**: essential for cytokine gene transcription (TNF-alpha, IL-6, IL-13)

### 1.3 Mediator Categories

#### Preformed mediators (released within seconds to minutes):

| Mediator | Storage | Key Actions |
|----------|---------|-------------|
| **Histamine** | Granule matrix (bound to heparin/serglycin) | Vasodilation (H1R on endothelium --> NO/PGI2), vascular permeability increase (endothelial cell contraction), smooth muscle contraction, pruritus (sensory nerve stimulation), gastric acid secretion (H2R) |
| **Tryptase** | Major granule serine protease (MCT and MCTC subtypes) | Cleaves fibrinogen, activates PAR-2 on epithelium (inflammation amplification), activates MMPs, degrades neuropeptides like VIP and CGRP |
| **Chymase** | In MCTC type mast cells; electron-dense granule regions | Converts angiotensin I --> angiotensin II, degrades extracellular matrix, activates IL-1beta and TGF-beta, processes MMP-9 |
| **TNF-alpha** | Both preformed (granule-stored) and newly synthesized | Endothelial activation (E-selectin, VCAM-1 upregulation), neutrophil recruitment, NF-kappaB activation in target cells |
| **Heparin** | Serglycin proteoglycan in granule matrix | Anticoagulant, packages/stabilizes other granule mediators (histamine, proteases), binds growth factors |

#### Newly synthesized lipid mediators (minutes to hours):

Generated from **arachidonic acid** released from membrane phospholipids by **cPLA2** (cytosolic phospholipase A2):

| Mediator | Synthetic Pathway | Key Actions |
|----------|-------------------|-------------|
| **PGD2** | COX-1/COX-2 --> PGH2 --> hematopoietic PGD synthase | Bronchoconstriction (DP1, TP receptors), vasodilation, eosinophil/Th2 chemotaxis (CRTH2/DP2), nasal congestion |
| **LTC4** | 5-LOX --> LTA4 --> LTC4 synthase | Parent cysteinyl leukotriene; exported via MRP1 transporter |
| **LTD4** | Extracellular conversion: LTC4 --> gamma-GT removes glutamate | Potent bronchoconstrictor (~1000x histamine), mucus secretion, vascular permeability. Acts via CysLT1 receptor |
| **LTE4** | LTD4 --> dipeptidase removes glycine | Longest-lived cysteinyl leukotriene; sustained bronchoconstriction, eosinophil recruitment |
| **PAF** (platelet-activating factor) | Acetylation of lyso-PAF by LPCAT2 | Platelet aggregation, neutrophil activation, bronchoconstriction, hypotension; extremely potent at picomolar concentrations |

#### Newly transcribed cytokines and chemokines (hours):

- **IL-4, IL-13**: Th2 amplification, B cell IgE switching, mucus production
- **IL-5**: Eosinophil differentiation, recruitment, survival
- **IL-6**: Acute phase response, B cell differentiation
- **IL-8 (CXCL8)**: Neutrophil chemotaxis
- **TNF-alpha** (de novo pool supplements preformed stores)
- **CCL2 (MCP-1), CCL3 (MIP-1alpha)**: Monocyte/macrophage recruitment
- **GM-CSF**: Eosinophil and neutrophil survival/activation

### 1.4 Late-Phase Response (4-12 hours post-exposure)

The late-phase reaction represents a transition from acute mast cell-driven events to a broader inflammatory infiltrate:

1. **Eosinophil recruitment**: Driven by:
   - **IL-5** (from mast cells and Th2 cells): eosinophil maturation, bone marrow release, tissue survival
   - **Eotaxin (CCL11)**: selective eosinophil chemotaxis via CCR3
   - **PGD2** via CRTH2 receptor on eosinophils
   - Endothelial adhesion molecules (VCAM-1, P-selectin) upregulated by TNF-alpha and IL-4
2. **Basophil influx**: Additional histamine and IL-4/IL-13 source
3. **Th2 lymphocyte recruitment**: Amplifies cytokine milieu, sustains IgE production
4. **Neutrophil infiltration**: Driven by LTB4, IL-8, PAF
5. **Tissue consequences**: Edema, mucus hypersecretion, sustained smooth muscle constriction, tissue remodeling initiation

---

## 2. How Allergy Drives Chronic Inflammation

### 2.1 Persistent Th2 Skewing and Systemic Effects

In atopic individuals, the immune system is chronically biased toward a **Type 2 inflammatory profile**. This is not merely localized but has systemic consequences:

- **GATA-3+ Th2 cells** dominate the CD4+ T cell pool in affected tissues, continuously producing IL-4, IL-5, and IL-13
- **IL-4/IL-13 create positive feedback loops**: they promote further Th2 differentiation while suppressing Th1 (IFN-gamma) and Treg (FOXP3+) responses
- **Systemic IgE elevation**: chronic Th2 cytokine exposure drives ongoing B cell class switching, elevating total and specific IgE
- **Bone marrow eosinophil production**: IL-5 drives eosinopoiesis; circulating eosinophilia is a hallmark of chronic allergic disease
- **ILC2 cells** (innate lymphoid cells type 2) respond to epithelial alarmins (TSLP, IL-25, IL-33) and amplify the Th2 cytokine milieu independently of adaptive immunity

### 2.2 Eosinophilic Tissue Infiltration and Cytotoxic Damage

Eosinophils are the primary effector cells in chronic allergic tissue damage. Their granules contain four major cationic proteins:

| Granule Protein | Mechanism of Damage |
|-----------------|---------------------|
| **MBP** (major basic protein) | Directly cytotoxic to epithelium; disrupts lipid bilayers, increases membrane permeability; activates mast cells and basophils (creating amplification loops); stimulates acetylcholine release (vagal hyperreactivity) |
| **ECP** (eosinophil cationic protein) | Pore-forming toxin (RNase superfamily); cytotoxic to epithelial cells, neurons, and myocardium; activates mast cell degranulation |
| **EPO** (eosinophil peroxidase) | Generates reactive oxygen species (HOBr from H2O2 + Br-); oxidative damage to lipids, proteins, DNA in surrounding tissue |
| **EDN** (eosinophil-derived neurotoxin) | RNase activity; neurotoxic; antiviral; activates dendritic cells via TLR2 |

In chronic allergic disease, persistent eosinophil infiltration causes:
- Epithelial desquamation and denudation
- Subepithelial fibrosis (via TGF-beta from eosinophils)
- Smooth muscle hypertrophy/hyperplasia
- Goblet cell metaplasia and mucus gland hypertrophy
- Basement membrane thickening
- Angiogenesis (VEGF release)

### 2.3 Epithelial Barrier Dysfunction

IL-4 and IL-13 directly compromise epithelial barrier integrity through multiple mechanisms:

- **Tight junction disassembly**: IL-13 downregulates the expression of **occludin**, **ZO-1** (zonula occludens-1), and **beta-catenin** in epithelial cells. Simultaneously, it upregulates the pore-forming tight junction protein **claudin-2**, which creates cation-selective paracellular channels that increase permeability.
- **Epithelial apoptosis**: IL-4/IL-13 induce epithelial cell apoptosis, creating physical gaps in the barrier.
- **Goblet cell hyperplasia**: IL-13-driven mucin overproduction (MUC5AC) alters the mucosal surface composition but does not compensate for barrier loss.
- **Reduced antimicrobial defense**: Th2 cytokines suppress production of antimicrobial peptides (beta-defensins, LL-37), increasing susceptibility to secondary infections.

The consequence is a **"leaky" epithelium** that allows greater allergen penetration into submucosal tissue, creating a positive feedback loop: more allergen exposure --> more mast cell/Th2 activation --> more barrier damage --> more allergen penetration.

### 2.4 Mucosal Remodeling and Impaired Nutrient Absorption

In the gastrointestinal tract specifically, chronic allergic inflammation drives:

- **Villous atrophy**: chronic eosinophilic infiltration damages intestinal villi, reducing absorptive surface area
- **Crypt hyperplasia**: compensatory proliferation that does not restore normal absorptive capacity
- **Submucosal fibrosis**: TGF-beta-driven collagen deposition thickens the lamina propria, increasing the diffusion distance for nutrients
- **Altered transporter expression**: inflammatory cytokines downregulate nutrient transporters (e.g., SGLT1 for glucose, DMT1 for iron, NPC1L1 for cholesterol)
- **Protein-losing enteropathy**: barrier disruption allows plasma protein leak into the gut lumen
- **Disrupted claudin expression**: loss of claudins 2 and 15 impairs paracellular Na+ flow critical for nutrient co-transport, which in severe cases can cause fatal malnutrition

### 2.5 Cross-Talk Between Allergic Inflammation and Metabolic Pathways

- **NF-kappaB as a shared node**: NF-kappaB drives both allergic inflammatory gene transcription and metabolic inflammation (insulin resistance, hepatic acute phase response)
- **IL-4/IL-13 effects on macrophage polarization**: promote M2 macrophage phenotype, which in adipose tissue links to altered insulin sensitivity
- **Histamine effects on metabolism**: H1R activation stimulates glycogenolysis; H3R in the hypothalamus regulates appetite and energy expenditure
- **Oxidative stress overlap**: ROS from eosinophil degranulation and mast cell activation damage mitochondrial function, impair beta-oxidation, and exhaust antioxidant reserves (glutathione, superoxide dismutase, catalase)
- **HPA axis stimulation**: chronic inflammation activates the hypothalamic-pituitary-adrenal axis, elevating endogenous cortisol, which has its own metabolic consequences (see Section 4)

---

## 3. Allergy's Effect on Nutrient Status

### 3.1 Histamine Metabolism: DAO Cofactor Consumption

Histamine is metabolized by two primary enzymes:

1. **DAO (diamine oxidase)** -- predominant in intestinal mucosa, kidney, placenta
   - Catalyzes oxidative deamination: histamine --> imidazole acetaldehyde + NH3 + H2O2
   - DAO is a **copper-containing amine oxidase** requiring three critical cofactors:

| Cofactor | Role in DAO Function | Depletion Mechanism in Allergy |
|----------|----------------------|-------------------------------|
| **Vitamin B6 (pyridoxal 5'-phosphate, PLP)** | Essential coenzyme for DAO catalytic activity | Consumed stoichiometrically during histamine metabolism; demand increases dramatically with high histamine load. PLP is also consumed by histidine decarboxylase (making more histamine) and by the >140 other PLP-dependent enzymes competing for the same pool |
| **Vitamin C (ascorbic acid)** | Cofactor for DAO; also directly degrades histamine via oxidation | Chronic mast cell activation generates massive ROS --> ascorbate consumed as antioxidant. Blood histamine levels inversely correlate with vitamin C levels. Ascorbate also required for collagen synthesis (impaired in chronically inflamed tissue) |
| **Copper (Cu2+)** | Prosthetic group in DAO active site (topaquinone cofactor requires Cu for biogenesis) | Copper deficiency directly reduces DAO production. Zinc supplementation (common in allergy) can antagonize copper absorption. Inflammatory states alter ceruloplasmin/copper metabolism |

2. **HNMT (histamine N-methyltransferase)** -- intracellular, especially in liver, kidney, CNS
   - Catalyzes: histamine + SAM --> N-methylhistamine + SAH
   - Consumes **S-adenosylmethionine (SAM)**, linking histamine metabolism to methylation capacity
   - SAM depends on folate, B12, and methionine availability
   - Chronic histamine load depletes SAM, potentially impairing other methylation-dependent processes (DNA methylation, creatine synthesis, phosphatidylcholine synthesis, neurotransmitter metabolism)

**The vicious cycle**: High histamine --> DAO/HNMT work overtime --> cofactors depleted --> DAO/HNMT activity falls --> histamine accumulates further --> symptoms worsen --> more mast cell activation from histamine-mediated inflammation.

### 3.2 Oxidative Stress and Antioxidant Depletion

Mast cell activation and eosinophil degranulation are potent generators of reactive oxygen species:

- **EPO** (eosinophil peroxidase): H2O2 + Br- --> HOBr (hypobromous acid) -- highly reactive oxidant
- **Mast cell-derived ROS**: NADPH oxidase activation during degranulation produces superoxide (O2-)
- **Leukotriene synthesis**: 5-LOX generates radical intermediates during LTA4 synthesis
- **DAO reaction product**: H2O2 from histamine degradation itself contributes to oxidative burden

Antioxidants consumed under chronic allergic oxidative stress:

| Antioxidant | Mechanism of Depletion | Consequences of Depletion |
|-------------|------------------------|---------------------------|
| **Glutathione (GSH)** | Consumed by glutathione peroxidase neutralizing H2O2 and lipid peroxides | Impaired xenobiotic detoxification, reduced leukotriene C4 synthesis (GSH is a substrate for LTC4 synthase, creating a paradoxical anti-inflammatory effect), mitochondrial vulnerability |
| **Vitamin C** | Direct scavenging of ROS; regeneration of vitamin E | Impaired collagen synthesis, reduced DAO function, scurvy-like symptoms with severe depletion |
| **Vitamin E (alpha-tocopherol)** | Lipid peroxyl radical chain-breaking in membranes | Membrane vulnerability to lipid peroxidation, neurological dysfunction |
| **Selenium (via glutathione peroxidase)** | Increased demand for selenoenzyme activity | Reduced GPx activity, thyroid dysfunction (selenoprotein deiodinases) |
| **Zinc** | SOD (Cu/Zn-SOD) cofactor; consumed under oxidative stress; redistributed to acute phase proteins | Impaired immune regulation (zinc is critical for Treg function), wound healing impairment, taste/smell dysfunction |

### 3.3 Dietary Restrictions and Nutritional Gaps

Allergen avoidance diets create predictable nutrient deficits:

| Avoided Allergen | Nutrients at Risk |
|------------------|-------------------|
| **Cow's milk** | Calcium, vitamin D, vitamin B2 (riboflavin), iodine, phosphorus, high-quality protein |
| **Eggs** | Vitamin D, B12, selenium, choline, complete protein with all essential amino acids |
| **Fish/shellfish** | Omega-3 fatty acids (EPA/DHA), vitamin D, iodine, selenium, zinc |
| **Tree nuts** | Vitamin E, magnesium, manganese, copper, monounsaturated fats |
| **Wheat** | B vitamins (thiamine, folate, niacin), iron, fiber, zinc |
| **Soy** | Complete protein, iron, calcium, B vitamins, isoflavones |
| **Multiple allergens** (common in children with EoE) | Compounded risk of protein-calorie malnutrition, micronutrient deficiencies, growth failure |

### 3.4 Gut Inflammation and Absorption Impairment

Food allergies and eosinophilic GI disorders create a direct barrier to nutrient absorption:

**Eosinophilic esophagitis (EoE)**:
- Dysphagia and feeding avoidance --> reduced caloric intake
- Stricture formation --> mechanical obstruction
- Elimination diets (often removing 6+ food groups) --> severe nutritional restriction
- In children: growth failure, weight loss, micronutrient deficiencies

**Eosinophilic gastroenteritis (EGE)**:
- **Mucosal involvement** (most common): malabsorption, protein-losing enteropathy, iron deficiency from chronic intestinal blood loss, folate deficiency
- **Muscular involvement**: obstruction, impaired motility --> maldigestion
- **Serosal involvement**: ascites with protein loss
- Documented deficiencies: iron, folate, vitamin D, albumin (hypoalbuminemia), fat-soluble vitamins

**Mechanism of malabsorption in allergic gut disease**:
1. Eosinophilic mucosal infiltration --> direct enterocyte damage
2. MBP/ECP release --> epithelial cytotoxicity
3. Tight junction disruption (IL-4/IL-13) --> loss of selective permeability
4. Villous blunting --> reduced surface area
5. Inflammation-driven increase in intestinal transit --> reduced contact time
6. Altered bile acid metabolism --> impaired fat/fat-soluble vitamin absorption

---

## 4. Glucocorticoid Treatment: Mechanism and Metabolic Consequences

### 4.1 Anti-Inflammatory Mechanism

Glucocorticoids (prednisone, prednisolone, dexamethasone, budesonide, etc.) exert their anti-inflammatory effects through both genomic and non-genomic mechanisms:

#### Genomic mechanisms (onset: hours)

**Transrepression** (primary anti-inflammatory mechanism):
- Glucocorticoid diffuses into the cell and binds the **glucocorticoid receptor (GR)** in the cytoplasm, causing dissociation of HSP90/HSP70 chaperone complex
- Ligand-bound GR translocates to the nucleus
- **GR monomers** physically interact with (tether to) activated **NF-kappaB** (p65/RelA subunit) and **AP-1** (c-Jun/c-Fos), preventing them from binding their target gene promoters
- This suppresses transcription of pro-inflammatory genes: COX-2, iNOS, TNF-alpha, IL-1beta, IL-6, IL-8, ICAM-1, E-selectin, MMP-9, and hundreds of others
- GR also recruits **HDAC2** (histone deacetylase 2) to inflammatory gene promoters, causing chromatin compaction and gene silencing

**Transactivation** (anti-inflammatory protein induction):
- GR homodimers bind **GRE** (glucocorticoid response elements) in promoters of anti-inflammatory genes:
  - **Annexin A1 (lipocortin-1)**: inhibits **cPLA2**, blocking release of arachidonic acid from membrane phospholipids. This is upstream of both the COX and LOX pathways, thereby suppressing ALL prostaglandins AND ALL leukotrienes simultaneously
  - **MAPK phosphatase-1 (MKP-1/DUSP1)**: dephosphorylates and inactivates p38 MAPK and JNK, reducing inflammatory signaling
  - **IkappaB-alpha**: sequesters NF-kappaB in the cytoplasm
  - **SLPI** (secretory leukocyte protease inhibitor): anti-inflammatory protease inhibitor
  - **IL-10**: anti-inflammatory cytokine (indirect upregulation)

#### Non-genomic mechanisms (onset: seconds to minutes):
- Direct membrane effects on mast cells: stabilization, reduced degranulation
- Interference with arachidonic acid release via membrane-associated GR
- Rapid effects on endothelial permeability (vasoconstriction)
- Modulation of T cell receptor signaling via membrane-proximal events

### 4.2 Metabolic Side Effects: The Price of Immunosuppression

Glucocorticoids produce profound metabolic disruption, particularly with systemic (oral/IV) and chronic use:

#### 4.2.1 Glucose and Insulin Metabolism

- **PEPCK (phosphoenolpyruvate carboxykinase) upregulation**: GR directly transactivates the PEPCK gene promoter, increasing the rate-limiting enzyme of gluconeogenesis
- **G6Pase (glucose-6-phosphatase) upregulation**: increases hepatic glucose output by dephosphorylating G6P for release into blood
- **Insulin resistance**: GCs impair insulin signaling at multiple levels:
  - Reduced IRS-1 phosphorylation
  - Decreased GLUT4 translocation in muscle and adipose tissue
  - Inhibition of Akt/PKB pathway
  - In the liver: GR activation inhibits the IR pathway and Akt activity, inducing FOXO1, which further stimulates PEPCK and G6Pase expression
- **Substrate mobilization**: protein catabolism (see below) and lipolysis provide amino acids and glycerol as gluconeogenic substrates
- **Net effect**: hyperglycemia, hyperinsulinemia, steroid-induced diabetes in susceptible individuals

#### 4.2.2 Protein Catabolism and Muscle Wasting

- GCs upregulate **ubiquitin-proteasome pathway** components (MuRF1, MAFbx/atrogin-1) in skeletal muscle
- Suppress **mTOR/S6K1** anabolic signaling
- Inhibit amino acid transport into muscle cells
- Upregulate **glutamine synthetase** --> amino acid mobilization to liver for gluconeogenesis
- **Clinical result**: proximal myopathy, sarcopenia, reduced functional capacity
- Amino acids consumed in gluconeogenesis are unavailable for:
  - Glutathione synthesis (glycine, cysteine, glutamate)
  - Neurotransmitter synthesis
  - Immune cell proliferation
  - Tissue repair

#### 4.2.3 Calcium and Bone Metabolism

Multiple mechanisms converge to cause glucocorticoid-induced osteoporosis (GIO):

| Mechanism | Biochemical Detail |
|-----------|-------------------|
| **Reduced intestinal Ca2+ absorption** | GCs suppress transcellular calcium transport by downregulating TRPV6 (calcium channel), calbindin-D9k (intracellular transport), and PMCA1b (basolateral extrusion). They also reduce VDR expression in enterocytes |
| **Increased renal Ca2+ excretion** | GCs decrease renal tubular calcium reabsorption (reduced TRPV5 and calbindin-D28k expression) |
| **Osteoblast suppression** | GCs induce osteoblast apoptosis and inhibit osteoblast differentiation (suppress Wnt/beta-catenin signaling by upregulating DKK1 and sclerostin). Reduce BMP-2, osteocalcin, and type I collagen expression |
| **Osteoclast activation** (early) | GCs increase RANKL and decrease OPG expression in osteoblasts, shifting the RANKL/OPG ratio to favor osteoclastogenesis |
| **Osteocyte apoptosis** | Direct cytotoxic effect on osteocytes, disrupting the mechanosensory network |
| **Net effect** | Rapid bone loss (especially trabecular), vertebral fractures, avascular necrosis of femoral head |

#### 4.2.4 Magnesium Wasting

- GCs increase **renal magnesium excretion** through reduced reabsorption in the loop of Henle and distal convoluted tubule (downregulation of TRPM6 channel)
- Bone resorption releases Mg2+ but this is lost renally rather than retained
- Magnesium depletion consequences:
  - Impaired ATP production (Mg-ATP is the true substrate for kinases)
  - Worsened insulin resistance (Mg2+ is required for insulin receptor tyrosine kinase activity)
  - Increased neuromuscular excitability (cramps, arrhythmias)
  - Impaired PTH secretion and action (perpetuating calcium disorders)
  - Exacerbation of airway hyperreactivity in asthma (smooth muscle calcium handling)

#### 4.2.5 Potassium Loss

- GCs have variable **mineralocorticoid activity** (cortisol binds MR; prednisone has moderate MR activity; dexamethasone has minimal MR activity)
- Mineralocorticoid effect: activates ENaC (epithelial sodium channel) and ROMK (renal outer medullary K+ channel) in the collecting duct --> Na+ reabsorption, K+ and H+ secretion
- Additionally, GC-induced insulin secretion drives cellular K+ uptake, masking true total-body K+ depletion until it becomes severe
- Hypokalemia consequences: muscle weakness (compounding myopathy), cardiac arrhythmias, metabolic alkalosis, impaired insulin secretion (K+ needed for beta-cell function)

#### 4.2.6 Vitamin D Accelerated Catabolism

- GCs induce **CYP24A1** (25-hydroxyvitamin D 24-hydroxylase) -- the enzyme that catabolizes both 25(OH)D3 (calcidiol) and 1,25(OH)2D3 (calcitriol) into inactive metabolites via 24-hydroxylation
- Note: Some studies in osteoblast cell lines show GCs may suppress CYP24A1 in certain contexts, but the net systemic effect in patients on chronic GCs is reduced vitamin D levels, as demonstrated by NHANES data showing significant association between GC use and low 25(OH)D
- GCs also reduce VDR expression in target tissues, creating functional vitamin D resistance even if levels are adequate
- Combined with reduced intestinal calcium absorption, this creates a "double hit" on calcium homeostasis

#### 4.2.7 Zinc Depletion

- GC-induced HPA axis disruption alters zinc homeostasis (the HPA axis participates in zinc level regulation)
- Increased urinary zinc excretion
- Redistribution of zinc to hepatic metallothionein during the acute-phase response
- Zinc depletion consequences particularly relevant to allergy:
  - Impaired Treg cell function (zinc required for FOXP3 stability)
  - Reduced Cu/Zn-SOD antioxidant activity
  - Impaired wound healing (zinc-dependent MMPs and collagen synthesis)
  - Compromised taste and smell (zinc-dependent)
  - Impaired DAO function (zinc participates in DAO regulation)

#### 4.2.8 B Vitamin Increased Demand

- Gluconeogenesis is a B-vitamin-intensive process:
  - **B1 (thiamine)**: pyruvate dehydrogenase complex, alpha-ketoglutarate dehydrogenase
  - **B2 (riboflavin)**: FAD-dependent enzymes in gluconeogenesis and electron transport
  - **B3 (niacin)**: NAD+/NADH cycling in gluconeogenic and TCA reactions
  - **B5 (pantothenic acid)**: CoA synthesis for acetyl-CoA, succinyl-CoA
  - **B6 (PLP)**: aminotransferases converting amino acids to gluconeogenic intermediates (already depleted by histamine metabolism -- double hit)
  - **Biotin**: pyruvate carboxylase cofactor (first committed step of gluconeogenesis)
  - **Folate**: one-carbon metabolism supporting the increased nucleotide turnover from cell proliferation/repair
- Increased protein catabolism also demands B6 for transamination of mobilized amino acids

#### 4.2.9 HPA Axis Suppression

- Exogenous GCs suppress CRH (hypothalamus) and ACTH (anterior pituitary) via negative feedback
- Chronic use causes **adrenal cortex atrophy** (zona fasciculata and reticularis)
- Abrupt withdrawal risks **adrenal crisis**: hypotension, hypoglycemia, hyponatremia, hyperkalemia, shock
- Recovery of the HPA axis can take 6-12 months after prolonged GC therapy
- During suppression, the body cannot mount an appropriate cortisol stress response, increasing vulnerability to infection, surgery, and physiological stress

#### 4.2.10 Skin and Connective Tissue

- GCs inhibit fibroblast proliferation and collagen synthesis (types I and III)
- Suppress hyaluronic acid production
- Reduce glycosaminoglycan content in dermis
- Inhibit angiogenesis
- **Clinical result**: skin thinning (dermal atrophy), purpura, striae, delayed wound healing, increased infection risk

### 4.3 The Glucocorticoid Paradox

This creates a fundamental therapeutic tension:

> **Glucocorticoids effectively suppress allergic inflammation** (via NF-kappaB/AP-1 transrepression, annexin A1-mediated PLA2 inhibition, eosinophil apoptosis) **but simultaneously create nutrient depletions** (Ca, Mg, K, Zn, vitamin D, B vitamins, vitamin C) **that themselves promote inflammation and impair the body's ability to resolve it.**

Specific examples of the paradox:
- GCs deplete magnesium --> magnesium deficiency increases NF-kappaB activity and IL-6 production --> more inflammation
- GCs deplete zinc --> zinc deficiency impairs Treg function --> Th2 skewing worsens --> more allergic inflammation
- GCs deplete vitamin D --> vitamin D deficiency impairs antimicrobial peptide production and Treg induction --> increased infection and reduced immune tolerance
- GCs cause insulin resistance --> hyperglycemia promotes glycation of proteins, oxidative stress, and NF-kappaB activation --> metabolic inflammation
- GCs deplete B6 --> impaired DAO function --> histamine accumulates --> worsened allergy symptoms requiring more GCs

---

## 5. Antihistamine and Related Drug Effects

### 5.1 H1 Receptor Blockers

#### Mechanism: Inverse Agonism

Modern understanding recognizes that H1 "antihistamines" are not simple competitive antagonists but **inverse agonists**. The H1 receptor has constitutive (basal) activity -- it signals at a low level even without histamine bound. H1 antihistamines:
- Bind to a different site on the receptor than histamine
- **Stabilize the receptor in its inactive conformation**
- Reduce constitutive signaling below baseline, in addition to preventing histamine activation
- This is pharmacologically distinct from neutral antagonism (which would only block histamine without affecting basal activity)

#### What H1 Blockers Block:

| Effect | Mechanism | Clinical Manifestation Blocked |
|--------|-----------|-------------------------------|
| Vasodilation | H1R on endothelium --> eNOS --> NO production | Flushing, nasal congestion (partially) |
| Vascular permeability | H1R on endothelial cells --> contraction --> gap formation | Edema, urticaria, rhinorrhea |
| Pruritus | H1R on sensory nerve C-fibers --> depolarization | Itching (skin, eyes, nose) |
| Bronchoconstriction | H1R on airway smooth muscle --> Gq --> PLC --> Ca2+ --> contraction | Histamine-mediated airway narrowing |
| Smooth muscle contraction | H1R on GI/genitourinary smooth muscle | Cramping (partially) |

#### What H1 Blockers Do NOT Block:

- **Mast cell degranulation itself**: antihistamines do not prevent the release of histamine or other mediators; they only block downstream signaling at H1R
- **Late-phase inflammatory response**: eosinophil recruitment, Th2 cytokine production, and tissue remodeling are not histamine-dependent processes
- **Leukotriene effects**: bronchoconstriction, mucus production, and eosinophil chemotaxis driven by CysLTs are unaffected by H1 blockers
- **Prostaglandin effects**: PGD2-mediated vasodilation, pain, and eosinophil recruitment are unaffected
- **PAF effects**: systemic hypotension, platelet aggregation, and bronchoconstriction via PAF are unaffected
- **Tryptase/chymase tissue destruction**: protease-mediated effects continue
- **H2R, H3R, H4R-mediated effects**: each receptor subtype mediates distinct effects (gastric acid, neurotransmission, eosinophil chemotaxis respectively) not blocked by H1-selective drugs

#### First-Generation H1 Antihistamines (diphenhydramine, chlorpheniramine, hydroxyzine, promethazine):

**Poor selectivity** -- these drugs bind not only H1R but also:
- **Muscarinic acetylcholine receptors (mAChR)**: anticholinergic effects
  - Reduced GI motility and secretions --> impaired digestion and nutrient mixing
  - Dry mouth --> impaired starch digestion (reduced salivary amylase)
  - Urinary retention
  - Tachycardia
  - Mydriasis and blurred vision (cycloplegia)
- **Serotonin receptors (5-HT)**: appetite changes, sedation
- **Alpha-adrenergic receptors**: orthostatic hypotension
- **CNS H1 receptors**: freely cross BBB --> sedation, cognitive impairment, impaired psychomotor function
- **hERG potassium channels**: cardiac QT prolongation risk (terfenadine was withdrawn for this reason)

**Nutritional implications of anticholinergic effects**:
- Reduced gastric acid secretion (muscarinic M1/M3) --> impaired protein digestion and mineral absorption
- Reduced intestinal motility --> constipation, altered gut microbiome, reduced nutrient exposure to absorptive surface
- Reduced pancreatic secretion --> impaired fat and protein digestion

#### Second-Generation H1 Antihistamines (cetirizine, loratadine, fexofenadine, desloratadine, levocetirizine, bilastine):

- **High H1R selectivity**: negligible binding to muscarinic, serotonergic, or adrenergic receptors
- **Minimal BBB penetration**: substrates for P-glycoprotein efflux pump at the blood-brain barrier
- **No significant anticholinergic effects**: cetirizine and fexofenadine showed no anticholinergic activity even at supratherapeutic concentrations
- **Minimal nutritional impact**: no significant direct effects on GI motility, gastric acid, or nutrient absorption
- Some second-generation antihistamines (cetirizine, desloratadine) demonstrate modest **anti-inflammatory properties** beyond pure H1 blockade: reduced ICAM-1 expression, reduced eosinophil chemotaxis, reduced NF-kappaB activation. This is thought to be related to their inverse agonist activity at constitutively active H1R

### 5.2 H2 Receptor Blockers (famotidine, ranitidine, cimetidine)

H2 blockers are used in allergy primarily for:
- Refractory urticaria (additive effect with H1 blockers)
- Mast cell activation syndrome (MCAS) -- reducing gastric symptoms
- Systemic mastocytosis

#### Mechanism:
- Competitive (inverse) agonism at **H2 receptors** on gastric parietal cells
- H2R --> Gs --> adenylyl cyclase --> cAMP --> PKA --> H+/K+-ATPase (proton pump) activation
- H2 blockers reduce this cascade --> reduced gastric acid secretion (but less potently than PPIs)

#### Nutritional Consequences of Reduced Gastric Acid:

| Nutrient | Mechanism of Impaired Absorption | Clinical Significance |
|----------|----------------------------------|----------------------|
| **Vitamin B12** | Gastric acid + pepsin required to release B12 from food proteins. Reduced acid --> B12 remains protein-bound --> cannot bind intrinsic factor. Studies show ~50% of elderly patients on prolonged acid suppression have deficient/insufficient B12 levels. Note: supplemental B12 (free-form) is not affected | Megaloblastic anemia, peripheral neuropathy, cognitive impairment, elevated homocysteine |
| **Iron (non-heme)** | Gastric acid converts Fe3+ (ferric) to Fe2+ (ferrous), the absorbable form. Also, acid pH keeps iron soluble. Reduced acid --> iron remains insoluble Fe3+ | Iron deficiency anemia (particularly relevant: already depleted by eosinophilic GI bleeding in allergic patients) |
| **Calcium** | Acid pH required for dissolution of calcium salts (especially calcium carbonate). Reduced acid --> calcium remains insoluble in the gut lumen | Osteoporosis risk (compounding GC effects if co-administered) |
| **Magnesium** | Similar pH-dependent solubility as calcium | Hypomagnesemia (compounding GC renal wasting if co-administered) |
| **Folate** | Acid pH assists in deconjugation of food polyglutamate folate to absorbable monoglutamate form | Megaloblastic anemia, elevated homocysteine, impaired DNA synthesis |
| **Zinc** | Acid environment facilitates zinc solubilization from food matrices | Impaired immune function, wound healing (compounding GC and allergy-driven depletion) |

### 5.3 Leukotriene Receptor Antagonists (Montelukast, Zafirlukast)

#### Mechanism:
- Selective competitive antagonist at the **CysLT1 receptor**
- Blocks the effects of LTC4, LTD4, and LTE4 (but not their synthesis)
- CysLT1R is Gq-coupled --> PLC --> IP3/DAG --> Ca2+ --> smooth muscle contraction, mucus secretion, edema

#### What Montelukast Blocks:
- Leukotriene-driven bronchoconstriction
- Mucus hypersecretion
- Airway edema
- Eosinophil recruitment to airways (CysLTs are eosinophil chemoattractants)
- Late-phase inflammatory response component (where H1 antihistamines fail)

#### Neuropsychiatric Effects:
- **2020 FDA boxed warning**: agitation, aggression, depression, suicidal ideation, hallucinations, insomnia, anxiousness
- **Mechanism remains incompletely understood** but hypotheses include:
  - **CysLT1 and CysLT2 receptors are expressed in the CNS**: neurons, astrocytes, microglia, cerebrovascular endothelium
  - CysLTs normally modulate neuroinflammation, blood-brain barrier integrity, and neuronal survival
  - Blocking CysLT1 in the CNS may disrupt leukotriene-mediated neuroprotective signaling
  - Altered serotonergic and dopaminergic signaling (though direct receptor binding at these receptors is not established)
  - Montelukast crosses the BBB; CNS concentrations may be sufficient for receptor occupancy
  - Interestingly, preclinical studies show montelukast has neuroprotective effects in some models (Alzheimer's, stroke), suggesting the neuropsychiatric effects may reflect disruption of homeostatic CysLT signaling in susceptible individuals

#### Nutritional implications:
- Montelukast itself has no known direct nutrient-depleting effects
- By reducing leukotriene-driven inflammation, it may indirectly improve gut barrier function and nutrient absorption in patients with allergic GI inflammation

### 5.4 Mast Cell Stabilizers (Cromolyn Sodium, Nedocromil)

#### Mechanism:
- **Prevent mast cell degranulation upstream** rather than blocking downstream mediator effects
- Primary mechanism: inhibition of calcium influx into mast cells
  - Cromolyn binds calcium-binding proteins on the mast cell membrane, forming a ternary complex with Ca2+ that blocks calcium entry
  - Also inhibits intermediate-conductance chloride channels on the cytoplasmic side of mast cells, which indirectly affects calcium flux through membrane potential changes
  - The net effect is prevention of the Ca2+-dependent granule fusion machinery from activating
- By preventing degranulation, cromolyn blocks the release of ALL preformed mediators (histamine, tryptase, chymase, TNF-alpha, heparin) and reduces de novo synthesis of lipid mediators
- Also inhibits activation of other inflammatory cells: eosinophils, monocytes, neutrophils (though less potently)

#### Advantages:
- No systemic absorption (used topically: inhaled, nasal, ophthalmic, oral for GI)
- Oral cromolyn (Gastrocrom) acts locally in the GI tract -- useful for mast cell-mediated food allergy, MCAS
- No significant nutrient-depleting effects
- No CNS effects

#### Limitations:
- Must be used **prophylactically** (before allergen exposure); ineffective once degranulation has occurred
- Requires frequent dosing (short half-life)
- Less potent than glucocorticoids for severe inflammation
- Exact mechanism still debated

---

## 6. The Combined Picture: Compounded Nutrient Challenges in the Treated Allergic Patient

### 6.1 Scenario: The Allergic Patient on Glucocorticoids + Antihistamines

Consider a patient with moderate-severe allergic disease (e.g., allergic asthma with eosinophilic esophagitis) being treated with:
- **Oral prednisone** (or equivalent systemic GC)
- **Cetirizine** (H1 blocker)
- **Famotidine** (H2 blocker for GI symptoms)
- **Montelukast** (leukotriene antagonist)
- **Elimination diet** (dairy-free, egg-free, wheat-free)

### 6.2 Convergent Nutrient Depletion Map

The following table maps how each factor contributes to depletion of the same critical nutrients -- creating compounded deficits:

| Nutrient | Allergy Itself | Glucocorticoid | H2 Blocker | Diet Restriction | Net Risk |
|----------|---------------|----------------|------------|------------------|----------|
| **Calcium** | Gut malabsorption from inflammation | Reduced intestinal absorption + increased renal excretion + bone resorption | Reduced acid-dependent dissolution | Dairy-free diet removes primary Ca source | **CRITICAL** -- 4-way hit |
| **Vitamin D** | Reduced absorption (fat-soluble vitamin, impaired in gut inflammation) | Accelerated catabolism (CYP24A1), reduced VDR expression | Minimal direct effect | Dairy/egg/fish avoidance removes dietary sources | **CRITICAL** -- 3-way hit |
| **Magnesium** | Gut inflammation impairs absorption; oxidative stress depletes | Increased renal excretion (TRPM6 downregulation) | Reduced acid-dependent solubilization | Nut/whole grain avoidance if applicable | **SEVERE** -- 3-way hit |
| **Zinc** | Consumed by oxidative stress (SOD), redistributed to acute phase; DAO function demands zinc | Increased urinary excretion; HPA disruption | Reduced acid-dependent solubilization | Meat/shellfish/wheat avoidance | **SEVERE** -- 3-way hit |
| **Iron** | GI blood loss (eosinophilic infiltration); reduced absorption in inflamed gut | Minimal direct effect (some impaired utilization) | Reduced ferric-to-ferrous conversion | Meat/wheat avoidance | **SEVERE** -- 3-way hit |
| **Vitamin B12** | Malabsorption from inflamed ileum/stomach | Increased demand from catabolism pathways | Reduced acid-dependent release from food proteins | Dairy/egg avoidance removes major sources | **SEVERE** -- 4-way hit |
| **Vitamin B6 (PLP)** | Consumed by DAO for histamine metabolism; consumed by HNMT indirectly (via SAM cycle) | Increased demand from gluconeogenesis-related transamination | Minimal direct effect | Potentially from dietary restriction | **SEVERE** -- 2-way hit with high baseline demand |
| **Vitamin C** | Consumed by ROS scavenging; consumed as DAO cofactor; inversely correlated with histamine levels | Increased urinary loss; increased oxidative stress from metabolic disruption | Minimal direct effect | May be reduced if fruits restricted | **MODERATE-SEVERE** |
| **Potassium** | Diarrhea from GI inflammation | Mineralocorticoid-effect renal wasting | Minimal direct effect | Fruit/vegetable restriction if applicable | **MODERATE-SEVERE** |
| **Folate** | Malabsorption from inflamed gut | Increased demand from cell turnover/catabolism | Reduced acid-dependent deconjugation | Wheat/legume avoidance | **MODERATE-SEVERE** |
| **Copper** | Required for DAO (depleted by high demand); zinc supplementation antagonizes copper absorption | Minimal direct effect | Minimal direct effect | Shellfish/nut avoidance | **MODERATE** |
| **Selenium** | Consumed by glutathione peroxidase under oxidative stress | Increased urinary excretion | Minimal direct effect | Fish/Brazil nut avoidance | **MODERATE** |
| **Protein** | Protein-losing enteropathy; malabsorption | Muscle catabolism increases protein turnover but depletes muscle stores | Reduced protein digestion (less pepsin activity) | Multiple allergen avoidance reduces protein sources | **MODERATE-SEVERE** |

### 6.3 Biochemical Cascade of Compounded Depletion

The interconnected consequences create self-reinforcing pathological cycles:

```
CYCLE 1: The Calcium-Bone-Vitamin D Spiral
    Glucocorticoid
        --> reduced intestinal Ca absorption
        --> increased renal Ca excretion
        --> accelerated vitamin D catabolism (CYP24A1)
        --> reduced VDR expression
    + H2 blocker --> reduced Ca dissolution
    + Dairy-free diet --> no dietary Ca source
    + Gut inflammation --> malabsorption of Ca and vitamin D
        --> severe negative calcium balance
        --> PTH elevation (secondary hyperparathyroidism)
        --> accelerated bone resorption to maintain serum Ca
        --> osteoporosis, fracture risk
        --> Mg depletion impairs PTH secretion --> worsens everything

CYCLE 2: The Histamine-Nutrient Depletion Feedback Loop
    Chronic mast cell activation
        --> massive histamine release
        --> DAO works overtime --> consumes B6, vitamin C, copper
        --> HNMT works overtime --> consumes SAM (deplets folate, B12, methionine)
        --> DAO cofactor depletion --> DAO activity falls
        --> histamine accumulates
        --> MORE mast cell activation (histamine amplifies via H4R on mast cells)
        --> MORE histamine to degrade --> MORE cofactor consumption
    + Glucocorticoid B6 demand (transamination) --> LESS B6 for DAO
    + Oxidative stress consuming vitamin C --> LESS C for DAO
    + Zinc depletion from GCs --> impaired DAO regulation
        --> histamine intolerance on top of true allergy

CYCLE 3: The Inflammation-Nutrient-Inflammation Spiral
    Allergic inflammation
        --> nutrient depletion (Zn, Mg, vitamin D, antioxidants)
    Zinc deficiency
        --> impaired Treg (FOXP3 instability) --> loss of immune tolerance
        --> Th2 skewing worsens
    Magnesium deficiency
        --> NF-kappaB disinhibition --> increased pro-inflammatory cytokines
        --> impaired insulin signaling --> metabolic inflammation
    Vitamin D deficiency
        --> impaired Treg induction (vitamin D drives FOXP3 expression)
        --> reduced antimicrobial peptides --> increased infections
        --> infections trigger more inflammation
    Antioxidant depletion
        --> uncontrolled oxidative stress --> NF-kappaB activation
        --> tissue damage --> more inflammatory cell recruitment
    ALL OF THESE --> more inflammation --> more nutrient depletion --> cycle continues

CYCLE 4: The Glucocorticoid Metabolic Trap
    Glucocorticoid therapy
        --> insulin resistance + hyperglycemia
        --> hyperglycemia --> protein glycation --> AGE formation
        --> AGEs activate RAGE receptor --> NF-kappaB activation --> inflammation
        --> protein catabolism --> amino acid depletion
        --> glutathione synthesis impaired (less cysteine, glycine, glutamate)
        --> reduced antioxidant capacity --> more oxidative stress
        --> muscle wasting --> reduced physical activity --> further bone loss
        --> HPA suppression --> need continued GC therapy to avoid adrenal crisis
        --> the very drug treating inflammation becomes a source of it
```

### 6.4 Clinical Implications

**For the clinician managing the allergic patient on multiple medications:**

1. **Monitor**: serum calcium (ionized), 25(OH)D, magnesium, zinc, ferritin/iron studies, B12, folate, potassium, albumin, fasting glucose/HbA1c, bone density (DEXA)
2. **Supplement proactively** (not reactively):
   - Calcium: 1000-1200 mg/day (divided doses, citrate form preferred if on H2 blockers -- citrate does not require acid for dissolution, unlike carbonate)
   - Vitamin D: 1000-4000 IU/day (monitor 25(OH)D, target 40-60 ng/mL)
   - Magnesium: 300-400 mg/day (glycinate or citrate for bioavailability)
   - Zinc: 15-30 mg/day (balanced with 1-2 mg copper to prevent copper deficiency)
   - B-complex: particularly B6, B12, folate
   - Vitamin C: 500-1000 mg/day (divided doses)
3. **Minimize systemic GC exposure**: use topical/inhaled routes (budesonide, fluticasone) with lower systemic bioavailability when possible
4. **Consider mast cell stabilizers** (cromolyn) as steroid-sparing agents
5. **Dietary counseling**: ensure allergen-free diet is nutritionally complete; consider registered dietitian referral, especially for children on multiple food elimination diets
6. **Address the root**: biologic therapies (omalizumab/anti-IgE, dupilumab/anti-IL-4Ralpha, mepolizumab/anti-IL-5) can reduce disease burden without the metabolic costs of systemic GCs

---

## Summary Table: Medication Mechanisms and Gaps

| Drug Class | What It Targets | What It Misses | Nutrient Cost |
|------------|-----------------|----------------|---------------|
| **Glucocorticoids** | NF-kappaB, AP-1, PLA2 (via annexin A1) --> suppresses nearly all inflammatory pathways | Does not address Th2 skewing root cause; creates metabolic disease | Ca, Mg, K, Zn, vitamin D, B vitamins, vitamin C, protein |
| **H1 antihistamines (2nd gen)** | H1R-mediated vasodilation, permeability, pruritus | Late-phase response, leukotrienes, prostaglandins, mast cell degranulation, tissue remodeling | Minimal |
| **H1 antihistamines (1st gen)** | Same as 2nd gen + CNS sedation | Same as 2nd gen | Impaired GI motility/absorption (anticholinergic) |
| **H2 blockers** | Gastric acid secretion (and some mast cell-mediated GI symptoms) | Does not address systemic allergic inflammation | B12, iron, calcium, magnesium, folate, zinc |
| **Montelukast** | CysLT1-mediated bronchoconstriction, mucus, eosinophil recruitment | Does not block histamine, prostaglandins, or mast cell degranulation; neuropsychiatric risk | None known |
| **Cromolyn** | Mast cell degranulation (upstream prevention) | Only prophylactic; no effect on already-released mediators; less potent than GCs | None |

---

## Sources

- [Mechanisms of mast cell signaling in anaphylaxis - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2788154/)
- [IgE and mast cells in allergic disease - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3597223/)
- [Mast Cell Mediators: Their Differential Release and the Secretory Pathways Involved - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4231949/)
- [Biological implications of preformed mast cell mediators - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11114649/)
- [Type I Hypersensitivity Reaction - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK560561/)
- [Effector mechanisms in allergic reactions - Immunobiology - NCBI](https://www.ncbi.nlm.nih.gov/books/NBK27112/)
- [Involvement and repair of epithelial barrier dysfunction in allergic diseases - Frontiers](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1348272/full)
- [IL-4 and IL-13 cause barrier dysfunction in human airway epithelial cells - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3875607/)
- [Nasal Epithelial Barrier Integrity and Tight Junctions Disruption in Allergic Rhinitis - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8176953/)
- [Type 2 immunity in allergic diseases - Nature](https://www.nature.com/articles/s41423-025-01261-2)
- [Eosinophilic Inflammation in Allergic Asthma - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3627984/)
- [Malnutrition in Eosinophilic Gastrointestinal Disorders - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7824578/)
- [Eosinophilic Gastroenteritis - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK547729/)
- [Nutritional recommendations for patients undergoing prolonged glucocorticoid therapy - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9080102/)
- [Prednisone Nutrient Depletion - Dr. Megan](https://prednisonepharmacist.com/the-evidence/prednisone-nutrient-depletion/)
- [Molecular Mechanisms of Glucocorticoid-Induced Insulin Resistance - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7827500/)
- [Selective modulation of the glucocorticoid receptor: transrepression of NF-kappaB and AP-1 - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3889831/)
- [Annexin A1 - Wikipedia](https://en.wikipedia.org/wiki/Annexin_A1)
- [Effects of glucocorticoids on CYP24A1 in Saos-2 cells and primary human osteoblasts](https://www.sciencedirect.com/science/article/abs/pii/S0303720719302278)
- [Association of Glucocorticoid Use and Low 25-Hydroxyvitamin D Levels: NHANES 2001-2006 - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3232615/)
- [Calcium and vitamin D for corticosteroid-induced osteoporosis - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7046131/)
- [H1 Antihistamines: Current Status and Future Directions - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3650962/)
- [Pharmacology of Antihistamines - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3666185/)
- [Antihistamines - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK538188/)
- [Proton pump inhibitors and risk of vitamin and mineral deficiency - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4110863/)
- [Famotidine and Vitamin B12 Interactions - Drugs.com](https://www.drugs.com/drug-interactions/famotidine-with-vitamin-b12-1066-0-754-3756.html?professional=1)
- [Leukotriene Receptor Antagonists - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK554445/)
- [Neuromodulatory effects of leukotriene receptor antagonists - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0014299924004436)
- [Cromolyn Sodium - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK557473/)
- [Mast cell stabilizers: from pathogenic roles to targeting therapies - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324444/)
- [DAO Deficiency and Histamine - MTHFR Support](https://www.mthfrsupport.com.au/2016/09/dao-deficiency-and-histamine-the-unlikely-connection/)
- [Histamine intolerance - BIOGENA](https://biogena.com/en/knowledge/guide/histamine-intolerance_bba_4070385)
