# Allergy Biochemistry Research

Understanding the full biochemical chain from food intake and nutrient absorption through inflammation, allergy, and medication effects — with a focus on what goes wrong, why recovery is slow, and what the body actually needs to heal.

## Start Here

If you're new to this repo, read in this order:

1. **[docs/README.md](docs/README.md)** — Documentation index with reading paths for different questions
2. **[Allergy-to-NARE Transition](docs/allergy-to-nare-transition.md)** — If your main question is "why do I still feel terrible after allergy season?"
3. **[Inflammation & Nutrient Cycles](docs/inflammation-biochemistry-nutrient-cycles.md)** — The core document explaining how inflammation depletes nutrients and creates self-sustaining vicious cycles, with a symptom → biochemistry mapping table

## Research Results

All findings are in `docs/` as detailed markdown documents (~250 KB total):

| Document | What it answers |
|----------|----------------|
| [Hormonal Cascade of Food Intake](docs/hormonal-biochemical-cascade-food-intake.md) | What happens biochemically when you eat? Insulin, glucagon, incretins, satiety hormones, cortisol, thyroid |
| [Nutrient Requirements & Inflammation](docs/nutrient-requirements-and-inflammation.md) | How much of each nutrient do you need? How does inflammation change that? Quantitative tables, depletion mechanisms |
| [Inflammation & Nutrient Depletion Cycles](docs/inflammation-biochemistry-nutrient-cycles.md) | How does inflammation become chronic? What vicious cycles form? Symptom-to-biochemistry mapping |
| [Allergy, Medications & Nutrition](docs/allergy-inflammation-nutrition-medications.md) | How do allergy, cortisol, and antihistamines compound nutrient depletion? |
| [Allergy-to-NARE Transition](docs/allergy-to-nare-transition.md) | Why do symptoms persist after pollen is gone? How to distinguish lingering allergy from NARE |
| [NARE Full Reference](docs/NARE-non-allergic-rhinitis-eosinophilia-syndrome.md) | Complete NARE pathophysiology, ILC2/alarmin axis, diagnosis, treatment |
| [Post-Allergy Recovery](docs/post-allergy-recovery-biochemistry.md) | What the body needs to return to baseline, 12-nutrient supplementation protocol, 12-week recovery timeline |

## Project Structure

```
allergy-biochemistry/
├── docs/                 # Research documents (main output)
│   ├── README.md         # Index with reading paths and cross-references
│   └── *.md              # Detailed biochemistry documents
├── data/                 # Raw and processed datasets (future)
│   ├── raw/
│   └── processed/
├── notebooks/            # Jupyter analysis notebooks (future)
├── scripts/              # Analysis and pipeline scripts (future)
├── results/              # Figures, tables, and reports (future)
│   ├── figures/
│   └── tables/
├── literature/           # Literature reviews and notes (future)
└── .claude/skills/       # 161 scientific skills (K-Dense) scoped to this project
```

**`docs/`** is where all current research lives. Other directories are scaffolded for future computational analysis (pathway modeling, gene expression datasets, protein structure work).

## How Research Is Done

This project uses **Claude Code** with **[Claude Scientific Skills](https://github.com/K-Dense-AI/claude-scientific-skills)** (161 skills installed at `.claude/skills/`). The skills provide optimized access to scientific databases and packages:

### Databases used
- **PubMed / bioRxiv** — Literature search, systematic reviews
- **UniProt / PDB / AlphaFold** — Protein structure and function
- **KEGG / Reactome / STRING** — Pathway mapping, protein interaction networks
- **ChEMBL / PubChem / DrugBank** — Drug and compound data
- **ClinVar / ClinicalTrials.gov** — Clinical variant and trial data
- **GEO / GTEx** — Gene expression datasets
- **gnomAD / Ensembl** — Genomic variation and annotation

### Computational tools (for future analysis)
- **BioPython / RDKit** — Sequence and molecular analysis
- **Scanpy / AnnData** — Single-cell transcriptomics
- **matplotlib / seaborn / plotly** — Visualization
- **statsmodels / scikit-learn** — Statistical modeling

### Research workflow
1. Define research question
2. AI-assisted literature synthesis across multiple databases
3. Biochemical pathway mapping with enzyme/receptor specificity
4. Cross-referencing across documents to build connected understanding
5. Document findings with full pathway detail (enzyme names, receptor subtypes, signaling cascades)
6. Future: validate key claims computationally using expression data, protein structures, pathway models

## Research Areas

### Core Chain: Food → Hormones → Nutrients → Inflammation → Allergy
- Hormonal response to food intake (insulin/glucagon, incretins, mTOR/AMPK)
- Nutrient requirements and biochemical roles (minerals, vitamins, cofactors)
- Inflammation biochemistry (NF-κB, COX/LOX, cytokines, resolution)
- Chronic inflammation and nutrient depletion vicious cycles

### Allergy-Specific
- IgE-mediated hypersensitivity (FcεRI signaling, mast cell degranulation)
- Mast cell and basophil mediators (histamine, tryptase, leukotrienes)
- Th2/ILC2 immune circuits and eosinophilic inflammation
- Epithelial barrier dysfunction and mucosal remodeling
- NARE (Non-Allergic Rhinitis with Eosinophilia Syndrome)
- Allergy-to-NARE transition mechanisms

### Medication Effects on Biochemistry
- Glucocorticoids: anti-inflammatory mechanism vs metabolic side effects
- Antihistamines: H1/H2 blockade and nutrient absorption consequences
- Leukotriene receptor antagonists
- The medication paradox: treating inflammation while depleting recovery nutrients

### Recovery
- Post-allergy nutrient debt and replenishment
- Supplementation protocols (forms, doses, timing, interactions)
- Recovery timeline mapping

## Getting Started (Computational Work — Optional)

Only needed if you want to run Python analysis scripts, Jupyter notebooks, or query databases programmatically. Not required for reading the research docs.

```bash
# Requires uv (Python package manager): https://docs.astral.sh/uv/
uv venv .venv
source .venv/bin/activate

# Install core dependencies
uv pip install biopython rdkit scanpy matplotlib seaborn pandas numpy
```

## License

MIT
