# Gene Expression Mining for Allergy-Nutrient Connections

## Overview

This analysis mines gene expression data from NCBI GEO datasets and published literature to identify how allergic rhinitis (AR) dysregulates nutrient-dependent enzymes and pathways. The goal: connect transcriptomic changes in allergic nasal tissue to increased or altered micronutrient demands.

## Data Sources

### GEO Datasets Analyzed

| Dataset | Title | Tissue | Comparison | Platform | Samples |
|---------|-------|--------|------------|----------|---------|
| **GSE101720** | GENEBRO: Gene expression profiling of lung epithelium | Nasal epithelium (brushings) | AR + Rhinitis vs Controls | RNA-seq | 26 nasal (7 AR, 10 R, 9 C) |
| **GSE261239** | Nasal mucosa gene expression upon birch pollen provocation | Nasal mucosa (scrapings) | Allergic vs Non-allergic | RNA-seq (QuantSeq) | 46 (22 allergic, 24 non-allergic) |

### Literature Curation
43 gene expression entries curated from published differential expression studies in allergic rhinitis nasal tissue (meta-analyses, individual transcriptomic studies).

**Important note on dataset interpretation**: GSE101720 measures baseline nasal epithelial expression (not allergen-challenged), while GSE261239 includes timepoints after allergen provocation. Literature values often reflect allergen-challenged or seasonal exposure states. This means some alarmin/cytokine genes (TSLP, IL33, etc.) may appear lower in baseline GEO data but are strongly upregulated upon allergen exposure.

## Target Gene Sets (71 genes across 10 functional categories)

| Category | Genes | Relevance |
|----------|-------|-----------|
| Histamine metabolism | HDC, AOC1, HNMT, HRH1-4 | Histamine production/degradation/signaling |
| Nutrient transporters | SLC41A1, SLC30A1, SLC39A14, SLC11A2, SLC40A1, TRPM6/7 | Mineral absorption/distribution |
| Inflammation/NF-kB | NFKB1, RELA, PTGS2, ALOX5, TNF, IL1B, IL6 | Core inflammatory machinery |
| Antioxidant defense | SOD1/2, GPX1/4, CAT, TXN, TXNRD1, GSR, GCLC/M | ROS neutralization |
| Mast cell/allergy | FCER1A, KIT, TPSAB1, CMA1, MS4A2 | IgE signaling, mast cell markers |
| Eosinophil | CCL11, CCR3, EPX, RNASE2/3, IL5RA, PRG2 | Eosinophil recruitment/function |
| Alarmin/ILC2 | TSLP, IL33, IL25, GATA3, IL5, IL13 | Epithelial alarmins, Th2 cytokines |
| Tryptophan metabolism | TPH1, IDO1, TDO2, KMO, KYNU, HAAO | Serotonin/kynurenine pathways |
| Methylation | MTHFR, MTR, MTRR, MAT1A, AHCY, CBS | One-carbon metabolism |
| Vitamin D | VDR, CYP27B1, CYP24A1, CYP2R1 | Local vitamin D activation/signaling |

## Key Findings

### 1. Tryptophan Steal Confirmed by IDO1 Upregulation

**IDO1** (indoleamine 2,3-dioxygenase) is among the most consistently upregulated genes across datasets (mean log2FC = +1.13). This enzyme:
- Requires **iron (heme)** as cofactor
- Shunts tryptophan into the kynurenine pathway, away from serotonin synthesis
- Is induced by IFN-gamma in allergic inflammation
- The downstream kynurenine pathway enzymes **KMO** (FAD/B2-dependent) and **KYNU** (B6-dependent) are also upregulated

**Metabolic consequence**: Tryptophan depletion, increased demand for iron, B2, B6. Potential serotonin deficit. Neurotoxic kynurenine metabolites accumulate.

### 2. Vitamin D Axis is Disrupted in Three Ways

| Gene | Direction | Cofactors | Effect |
|------|-----------|-----------|--------|
| VDR | Down | Vitamin D3, Zinc | Reduced receptor = impaired signaling even with adequate vitamin D |
| CYP27B1 | Down | Iron (heme) | Reduced local activation of 25(OH)D to active 1,25(OH)2D |
| CYP24A1 | Up | Iron (heme) | Increased inactivation of active vitamin D |
| CYP2R1 | Down | Iron (heme) | Reduced initial hydroxylation step |

**Metabolic consequence**: Triple hit on local vitamin D signaling -- less activation, more degradation, less receptor. All CYP enzymes require iron-heme, creating competition with other upregulated heme enzymes (IDO1, PTGS2, EPX).

### 3. Histamine: Double Hit from Production Up + Degradation Down

| Gene | Direction | Cofactors | Effect |
|------|-----------|-----------|--------|
| HDC | Up | Vitamin B6 (PLP) | More histamine production |
| AOC1 (DAO) | Down | Copper, B6 | Less histamine degradation (oxidative pathway) |
| HNMT | Down | SAM | Less histamine degradation (methylation pathway) |
| HRH1 | Up | -- | More histamine receptor expression |

**Metabolic consequence**: Histamine accumulates from both increased production and decreased clearance. B6 is consumed by upregulated HDC, while copper deficiency impairs DAO. SAM depletion impairs HNMT.

### 4. Antioxidant System: Upregulated but Substrate-Starved

Antioxidant enzyme genes are upregulated (a stress response), but this increases cofactor demand:

| Gene | FC | Cofactors | Demand |
|------|-----|-----------|--------|
| GPX1 | +0.57 | Selenium | Increased Se demand |
| SOD2 | +0.40 | Manganese | Increased Mn demand |
| SOD1 | +0.38 | Copper, Zinc | Increased Cu/Zn demand |
| GPX4 | +0.19 | Selenium | Increased Se demand |
| TXNRD1 | +0.27 (lit) | Selenium, B2 | Increased Se+B2 demand |
| GCLC | +0.22 (lit) | Cysteine | Increased cysteine demand for glutathione |

Meanwhile, **CAT** (catalase, iron-heme dependent) is downregulated, and **TXN** (thioredoxin, zinc-dependent) is strongly downregulated (FC = -0.97).

**Metabolic consequence**: The antioxidant system is trying to compensate but likely running out of cofactors (Se, Mn, Cu, Zn). Catalase loss means more H2O2 accumulation.

### 5. Magnesium Transport is Impaired

| Gene | Direction | Function |
|------|-----------|----------|
| TRPM6 | Down | Epithelial Mg channel (master regulator of Mg absorption) |
| SLC41A1 | Down | Mg transporter |
| TRPM7 | Down | Ubiquitous Mg/Ca channel-kinase |

**Metabolic consequence**: Reduced magnesium absorption at the mucosal level. Mg is needed for >300 enzymes including NF-kB regulation, and its deficit promotes inflammation.

### 6. Mast Cell Proteases Consume Zinc

**TPSAB1** (tryptase) is the most upregulated gene with nutrient cofactors (mean FC = +1.40). Tryptase and **CMA1** (chymase) are zinc-dependent serine proteases. Massive mast cell degranulation in AR tissue sequesters zinc into enzymatic pools, potentially depleting it from other zinc-dependent processes (SOD1, VDR zinc fingers, GATA3 zinc fingers, NF-kB zinc fingers).

### 7. Iron-Heme Competition

Multiple upregulated enzymes compete for the limited heme iron pool:

| Gene | FC | Function |
|------|-----|----------|
| IDO1 | +1.13 | Tryptophan catabolism |
| PTGS2 (COX-2) | +0.96 | Prostaglandin synthesis |
| EPX | +0.77 (lit) | Eosinophil peroxidase |
| ALOX5 | +0.50 (lit) | Leukotriene synthesis |

Meanwhile, iron-heme enzymes in protective pathways are downregulated:
- CYP27B1 (vitamin D activation)
- CAT (catalase, H2O2 clearance)
- CBS (transsulfuration)

**Metabolic consequence**: Inflammatory heme enzymes win the competition, starving protective/metabolic heme enzymes.

## Nutrient Demand Summary

### Nutrients Under INCREASED Demand

| Nutrient | Upregulated Dependent Genes | Impact |
|----------|---------------------------|--------|
| Iron (heme) | IDO1, PTGS2, EPX, ALOX5 | Inflammatory enzymes consume heme |
| Selenium | GPX1, GPX4, TXNRD1 | Antioxidant defense needs more Se |
| Zinc | SOD1, TPSAB1, CMA1, NFKB1, RELA | Proteases + TFs consume zinc |
| Manganese | SOD2 | MnSOD upregulated |
| Copper | SOD1 | Cu/ZnSOD upregulated |
| Vitamin B6 | HDC, KYNU, CBS | Histamine + kynurenine + transsulfuration |
| FAD (B2) | KMO, TXNRD1 | Kynurenine pathway + antioxidant |

### Nutrients with IMPAIRED Absorption/Signaling

| Nutrient | Mechanism | Downregulated Gene |
|----------|-----------|-------------------|
| Magnesium | Transporter downregulation | TRPM6, SLC41A1, TRPM7 |
| Vitamin D | Triple axis disruption | VDR, CYP27B1 down; CYP24A1 up |
| SAM/Methyl donors | HNMT down; impaired methylation | HNMT, MTRR |

## Novel Connections Discovered

1. **Heme iron competition model**: IDO1 and PTGS2 upregulation creates a "heme sink" that starves protective CYP enzymes (CYP27B1 for vitamin D, CYP2R1) and catalase. This links allergic inflammation directly to functional iron-heme deficiency in specific pathways.

2. **Zinc partitioning crisis**: Upregulation of zinc-dependent mast cell proteases (TPSAB1 +1.4, CMA1), NF-kB subunits, and SOD1 simultaneously increases zinc demand in inflammatory compartments while VDR and TXN (both zinc-dependent) are downregulated -- suggesting zinc is being competitively redistributed.

3. **B6 tug-of-war**: HDC (histamine synthesis, up) and KYNU (kynurenine metabolism, up) both require B6-PLP as cofactor. Meanwhile, AOC1/DAO (histamine degradation, also B6-dependent) is downregulated. B6 is being consumed by inflammatory pathways at the expense of protective ones.

4. **Magnesium absorption deficit**: All three major epithelial Mg transport/channel genes are downregulated. Since nasal and intestinal epithelia share transport mechanisms, this may reflect a systemic pattern that compounds Mg deficiency.

5. **Antioxidant futility cycle**: Antioxidant enzyme expression is upregulated (compensatory response) but their cofactors (Se, Mn, Cu, Zn) are likely being depleted by the very inflammation they're trying to counteract. This creates a futile cycle where more enzyme is made but cannot function.

## Methods

- GEO data downloaded via NCBI HTTPS API
- GSE101720: Pre-normalized expression matrix; nasal samples selected; Welch's t-test comparing all rhinitis (AR+R) vs controls
- GSE261239: VST-transformed counts; allergic vs non-allergic comparison across all timepoints
- Literature curation from published AR transcriptomic studies
- FDR correction via Benjamini-Hochberg per dataset
- Consensus direction determined by majority vote across datasets (threshold: |log2FC| > 0.3)

## Output Files

### Tables
- `all_gene_expression_results.csv` -- Full results from all sources (155 entries)
- `consensus_gene_expression.csv` -- Consensus across datasets per gene
- `nutrient_demand_implications.csv` -- Genes with cofactor demand annotations
- `GSE101720_AR_vs_Control.csv` -- Dataset-specific results
- `GSE101720_AllRhinitis_vs_Control.csv` -- Rhinitis (all) vs controls
- `GSE101720_Rhinitis_vs_Control.csv` -- Rhinitis-only vs controls

### Figures
- `gene_expression_heatmap.png` -- Heatmap across functional categories and datasets
- `gene_dotplot_nutrients.png` -- Dot plot of fold changes with cofactor annotations
- `pathway_nutrient_network.png` -- Network showing nutrient cofactors linked to dysregulated genes
- `nutrient_demand_summary.png` -- Bar chart of nutrient demand burden

## Script
`scripts/geo_allergy_gene_mining.py`
