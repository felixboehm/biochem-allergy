# Allergy Biochemistry Research

Research repository for investigating the biochemical mechanisms underlying allergic reactions, including IgE-mediated hypersensitivity, mast cell degranulation, histamine pathways, and immune signaling cascades.

## Research Areas

### Immunoglobulin E (IgE) Biology
- IgE structure and receptor binding (FcεRI, FcεRII/CD23)
- Class-switch recombination to IgE
- Allergen-IgE crosslinking mechanisms

### Mast Cell and Basophil Biochemistry
- Degranulation signaling pathways (Syk, LAT, PLCγ)
- Preformed mediators (histamine, tryptase, chymase)
- De novo synthesized mediators (prostaglandins, leukotrienes, cytokines)

### Histamine Pathways
- Histamine synthesis (histidine decarboxylase)
- Histamine receptors (H1R–H4R) and downstream signaling
- Histamine metabolism (DAO, HNMT)

### Cytokine Networks
- Th2 cytokines (IL-4, IL-5, IL-13) in allergic inflammation
- Alarmin signaling (TSLP, IL-25, IL-33)
- Regulatory T cell modulation

### Allergen Biochemistry
- Allergen protein families and structural motifs
- Protease allergens and epithelial barrier disruption
- Cross-reactivity and molecular mimicry

## Project Structure

```
allergy-biochemistry/
├── data/                # Raw and processed datasets
├── notebooks/           # Jupyter analysis notebooks
├── scripts/             # Analysis and pipeline scripts
├── results/             # Figures, tables, and reports
├── literature/          # Literature reviews and notes
└── docs/                # Documentation and protocols
```

## Tools & Databases

This project leverages Claude Scientific Skills for:
- **PubMed / bioRxiv** — Literature search and review
- **UniProt / PDB / AlphaFold** — Protein structure analysis
- **KEGG / Reactome** — Pathway mapping
- **ChEMBL / PubChem / DrugBank** — Drug and compound data
- **STRING** — Protein-protein interaction networks
- **ClinVar / ClinicalTrials.gov** — Clinical variant and trial data
- **GEO / GTEx** — Gene expression datasets
- **RDKit / BioPython / Scanpy** — Computational analysis

## Getting Started

```bash
# Create virtual environment
uv venv .venv
source .venv/bin/activate

# Install core dependencies
uv pip install biopython rdkit scanpy matplotlib seaborn pandas numpy
```

## License

MIT
