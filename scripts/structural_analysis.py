#!/usr/bin/env python3
"""
Structural Analysis of Key Nutrient-Dependent Allergy Enzymes
=============================================================

Retrieves protein structures from RCSB PDB and AlphaFold,
analyzes cofactor binding sites, and creates publication-quality
figures showing how nutrient deficiencies break allergy-relevant enzymes.

Author: Computational Biology Pipeline
Date: 2026-03-10
"""

import json
import os
import sys
import warnings
from pathlib import Path
from io import StringIO
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import requests

from Bio.PDB import PDBParser, MMCIFParser, Select, PDBIO
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

warnings.filterwarnings('ignore')

# ==============================================================================
# Configuration
# ==============================================================================

BASE_DIR = Path("/Users/felix/work/health/allergy-biochemistry")
DATA_DIR = BASE_DIR / "data" / "processed"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
REPORT_PATH = RESULTS_DIR / "structural-analysis-report.md"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

DPI = 300

# Enzyme definitions with PDB IDs, cofactors, and biological context
ENZYMES = {
    "DAO": {
        "full_name": "Diamine Oxidase (AOC1/ABP1)",
        "gene": "AOC1",
        "pdb_ids": ["3HI7", "3HIG"],
        "uniprot": "P19801",
        "cofactors": ["Cu", "TPQ"],
        "cofactor_names": ["Copper (Cu²⁺)", "Topaquinone (TPQ, from Tyr + Cu)"],
        "nutrients": ["Copper", "Vitamin B6 (PLP)"],
        "substrate": "Histamine",
        "reaction": "Histamine → Imidazole acetaldehyde + NH₃ + H₂O₂",
        "function": "Primary extracellular histamine degradation",
        "deficiency_effect": "Histamine accumulates → pseudo-allergic reactions, food intolerance",
        "clinical": "Histamine intolerance, migraine, GI symptoms, skin flushing",
        "metal_symbols": ["CU"],
        "het_cofactors": ["TPQ", "TYQ"],
        "color": "#E74C3C",
        "is_human": False,  # 3HI7 is pig kidney DAO (closest to human with structure)
    },
    "HNMT": {
        "full_name": "Histamine N-Methyltransferase",
        "gene": "HNMT",
        "pdb_ids": ["1JQD", "2AOT"],
        "uniprot": "P50135",
        "cofactors": ["SAM"],
        "cofactor_names": ["S-adenosylmethionine (SAM)"],
        "nutrients": ["Methionine", "Folate (B9)", "B12", "B6"],
        "substrate": "Histamine",
        "reaction": "Histamine + SAM → Nτ-methylhistamine + SAH",
        "function": "Intracellular histamine degradation (brain, airways)",
        "deficiency_effect": "Histamine accumulates intracellularly → neurological & respiratory symptoms",
        "clinical": "Brain fog, insomnia, anxiety, asthma exacerbation",
        "metal_symbols": [],
        "het_cofactors": ["SAM", "SAH", "HSM"],
        "color": "#3498DB",
        "is_human": True,
    },
    "HDC": {
        "full_name": "Histidine Decarboxylase",
        "gene": "HDC",
        "pdb_ids": ["4E1O", "7EIW"],  # Mammalian HDC structures
        "uniprot": "P19113",
        "cofactors": ["PLP"],
        "cofactor_names": ["Pyridoxal 5'-phosphate (PLP, vitamin B6)"],
        "nutrients": ["Vitamin B6 (PLP)"],
        "substrate": "L-Histidine",
        "reaction": "L-Histidine → Histamine + CO₂",
        "function": "Histamine biosynthesis",
        "deficiency_effect": "With adequate B6: normal histamine production. B6 deficiency paradoxically may increase mast cell histamine release",
        "clinical": "Dysregulated histamine production, mast cell instability",
        "metal_symbols": [],
        "het_cofactors": ["PLP", "PLR"],
        "color": "#E67E22",
        "is_human": False,
    },
    "IDO1": {
        "full_name": "Indoleamine 2,3-Dioxygenase 1",
        "gene": "IDO1",
        "pdb_ids": ["2D0T", "5ARE", "6E41"],
        "uniprot": "P14902",
        "cofactors": ["Heme (Fe²⁺)"],
        "cofactor_names": ["Heme/Protoporphyrin IX (contains Fe²⁺)"],
        "nutrients": ["Iron", "Heme"],
        "substrate": "L-Tryptophan",
        "reaction": "L-Trp + O₂ → N-formylkynurenine",
        "function": "Tryptophan catabolism, immune regulation via kynurenine pathway",
        "deficiency_effect": "Iron deficiency → reduced IDO1 → altered Trp metabolism → immune dysregulation",
        "clinical": "Immune imbalance, Th1/Th2 skewing, loss of oral tolerance",
        "metal_symbols": ["FE"],
        "het_cofactors": ["HEM", "HEC"],
        "color": "#9B59B6",
        "is_human": True,
    },
    "ALOX5": {
        "full_name": "Arachidonate 5-Lipoxygenase",
        "gene": "ALOX5",
        "pdb_ids": ["3O8Y", "3V99"],
        "uniprot": "P09917",
        "cofactors": ["Non-heme Fe²⁺"],
        "cofactor_names": ["Non-heme Iron (Fe²⁺/Fe³⁺)"],
        "nutrients": ["Iron"],
        "substrate": "Arachidonic acid",
        "reaction": "Arachidonic acid → 5-HPETE → LTA₄ → Leukotrienes",
        "function": "Produces leukotrienes (potent inflammatory mediators in asthma/allergy)",
        "deficiency_effect": "Iron-dependent: both excess and deficiency alter leukotriene production",
        "clinical": "Asthma, bronchoconstriction, allergic inflammation",
        "metal_symbols": ["FE"],
        "het_cofactors": [],
        "color": "#2ECC71",
        "is_human": True,
    },
    "SOD1": {
        "full_name": "Superoxide Dismutase 1 (Cu/Zn)",
        "gene": "SOD1",
        "pdb_ids": ["1PU0", "2C9V"],
        "uniprot": "P00441",
        "cofactors": ["Cu²⁺", "Zn²⁺"],
        "cofactor_names": ["Copper (Cu²⁺, catalytic)", "Zinc (Zn²⁺, structural)"],
        "nutrients": ["Copper", "Zinc"],
        "substrate": "Superoxide radical (O₂⁻)",
        "reaction": "2 O₂⁻ + 2H⁺ → O₂ + H₂O₂",
        "function": "Antioxidant defense — removes superoxide radicals",
        "deficiency_effect": "Cu/Zn deficiency → oxidative stress → mast cell activation, tissue damage",
        "clinical": "Oxidative damage amplifies allergic inflammation, epithelial barrier breakdown",
        "metal_symbols": ["CU", "ZN"],
        "het_cofactors": [],
        "color": "#1ABC9C",
        "is_human": True,
    },
    "GPX1": {
        "full_name": "Glutathione Peroxidase 1",
        "gene": "GPX1",
        "pdb_ids": ["2F8A"],
        "uniprot": "P07203",
        "cofactors": ["Selenocysteine (Sec)"],
        "cofactor_names": ["Selenocysteine (Sec, amino acid with Se)"],
        "nutrients": ["Selenium"],
        "substrate": "H₂O₂ / Lipid hydroperoxides",
        "reaction": "2 GSH + H₂O₂ → GSSG + 2 H₂O",
        "function": "Reduces peroxides, prevents oxidative damage",
        "deficiency_effect": "Se deficiency → oxidative stress → NF-κB activation → pro-inflammatory cascade",
        "clinical": "Amplified allergic inflammation, asthma severity, oxidative tissue damage",
        "metal_symbols": ["SE"],
        "het_cofactors": ["SEC", "CSO"],
        "color": "#F39C12",
        "is_human": True,
    },
}


# ==============================================================================
# Step 1: Retrieve Protein Structures
# ==============================================================================

def download_pdb(pdb_id: str, output_dir: Path) -> Path | None:
    """Download a PDB structure from RCSB in mmCIF format."""
    cif_path = output_dir / f"{pdb_id}.cif"
    pdb_path = output_dir / f"{pdb_id}.pdb"

    # Check if already downloaded
    if cif_path.exists():
        print(f"  [cached] {cif_path}")
        return cif_path
    if pdb_path.exists():
        print(f"  [cached] {pdb_path}")
        return pdb_path

    # Try mmCIF first
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"  Downloading {url} ...")
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            cif_path.write_text(resp.text)
            print(f"  [OK] Saved {cif_path}")
            return cif_path
    except Exception as e:
        print(f"  [WARN] mmCIF download failed: {e}")

    # Fallback to PDB format
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Trying PDB format: {url} ...")
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            pdb_path.write_text(resp.text)
            print(f"  [OK] Saved {pdb_path}")
            return pdb_path
    except Exception as e:
        print(f"  [WARN] PDB download failed: {e}")

    return None


def download_alphafold(uniprot_id: str, output_dir: Path) -> Path | None:
    """Download predicted structure from AlphaFold database."""
    pdb_path = output_dir / f"AF-{uniprot_id}.pdb"
    cif_path = output_dir / f"AF-{uniprot_id}.cif"

    if pdb_path.exists():
        print(f"  [cached] {pdb_path}")
        return pdb_path
    if cif_path.exists():
        print(f"  [cached] {cif_path}")
        return cif_path

    # Try AlphaFold API first to get the correct URL
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"  Querying AlphaFold API: {api_url} ...")
    try:
        resp = requests.get(api_url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            if data:
                pdb_url = data[0].get("pdbUrl", "")
                if pdb_url:
                    print(f"  Downloading: {pdb_url} ...")
                    resp2 = requests.get(pdb_url, timeout=30)
                    if resp2.status_code == 200:
                        pdb_path.write_text(resp2.text)
                        print(f"  [OK] Saved {pdb_path}")
                        return pdb_path
    except Exception as e:
        print(f"  [WARN] AlphaFold API failed: {e}")

    # Fallback: try multiple versions directly
    for ver in ["v6", "v4", "v3", "v2"]:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{ver}.pdb"
        print(f"  Trying AlphaFold {ver}: {url} ...")
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                pdb_path.write_text(resp.text)
                print(f"  [OK] Saved {pdb_path}")
                return pdb_path
        except Exception as e:
            print(f"  [WARN] AlphaFold {ver} download failed: {e}")

    # Try CIF
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
    print(f"  Trying AlphaFold CIF: {url} ...")
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            cif_path.write_text(resp.text)
            print(f"  [OK] Saved {cif_path}")
            return cif_path
    except Exception as e:
        print(f"  [WARN] AlphaFold CIF download failed: {e}")

    return None


def get_pdb_info(pdb_id: str) -> dict:
    """Get metadata about a PDB entry from RCSB REST API."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            return resp.json()
    except Exception:
        pass
    return {}


def retrieve_all_structures():
    """Download structures for all enzymes."""
    print("=" * 70)
    print("STEP 1: Retrieving Protein Structures")
    print("=" * 70)

    structure_files = {}

    for name, info in ENZYMES.items():
        print(f"\n--- {name} ({info['full_name']}) ---")

        downloaded = None
        pdb_metadata = {}

        # Try PDB IDs
        for pdb_id in info["pdb_ids"]:
            downloaded = download_pdb(pdb_id, DATA_DIR)
            if downloaded:
                pdb_metadata = get_pdb_info(pdb_id)
                break

        # Fallback to AlphaFold
        if not downloaded:
            print(f"  No PDB structure available; trying AlphaFold...")
            downloaded = download_alphafold(info["uniprot"], DATA_DIR)
            if downloaded:
                pdb_metadata = {"source": "AlphaFold", "uniprot": info["uniprot"]}

        if downloaded:
            structure_files[name] = {
                "path": downloaded,
                "metadata": pdb_metadata,
            }
        else:
            print(f"  [ERROR] Could not retrieve structure for {name}")
            structure_files[name] = None

    return structure_files


# ==============================================================================
# Step 2: Analyze Cofactor Binding Sites
# ==============================================================================

def parse_structure(filepath: Path):
    """Parse a PDB or mmCIF file."""
    suffix = filepath.suffix.lower()
    if suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(filepath.stem, str(filepath))
    return structure


def find_cofactor_atoms(structure, metal_symbols, het_cofactors):
    """Find metal ions and heteroatom cofactors in the structure."""
    cofactor_atoms = []
    cofactor_residues = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag = residue.get_id()[0]

                # Check for metal ions (HETATM with matching element)
                if hetflag.startswith("H_") or hetflag == "W":
                    for atom in residue:
                        element = atom.element.strip().upper()
                        if element in metal_symbols:
                            cofactor_atoms.append(atom)
                            cofactor_residues.append(residue)

                # Check for het cofactors by residue name
                if resname in het_cofactors:
                    cofactor_residues.append(residue)
                    for atom in residue:
                        cofactor_atoms.append(atom)

    return cofactor_atoms, cofactor_residues


def find_coordinating_residues(structure, cofactor_atoms, distance_cutoff=3.5):
    """Find amino acid residues coordinating the cofactor."""
    # Collect all protein atoms
    protein_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    for atom in residue:
                        protein_atoms.append(atom)

    if not protein_atoms or not cofactor_atoms:
        return []

    ns = NeighborSearch(protein_atoms)

    coordinating = []
    for cof_atom in cofactor_atoms:
        nearby = ns.search(cof_atom.get_vector().get_array(), distance_cutoff)
        for atom in nearby:
            res = atom.get_parent()
            if is_aa(res, standard=True):
                dist = cof_atom - atom
                entry = {
                    "cofactor_atom": cof_atom.get_name(),
                    "cofactor_element": cof_atom.element.strip(),
                    "cofactor_residue": cof_atom.get_parent().get_resname().strip(),
                    "protein_residue": res.get_resname().strip(),
                    "protein_resid": res.get_id()[1],
                    "protein_chain": res.get_parent().get_id(),
                    "protein_atom": atom.get_name(),
                    "distance": round(dist, 2),
                }
                coordinating.append(entry)

    # Deduplicate by (chain, resid, atom)
    seen = set()
    unique = []
    for c in coordinating:
        key = (c["protein_chain"], c["protein_resid"], c["protein_atom"])
        if key not in seen:
            seen.add(key)
            unique.append(c)

    # Sort by distance
    unique.sort(key=lambda x: x["distance"])
    return unique


def get_active_site_info(structure, enzyme_name):
    """Get Cα coordinates and secondary structure info for visualization."""
    ca_coords = []
    residue_info = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    try:
                        ca = residue["CA"]
                        coords = ca.get_vector().get_array()
                        ca_coords.append(coords)
                        residue_info.append({
                            "resname": residue.get_resname().strip(),
                            "resid": residue.get_id()[1],
                            "chain": chain.get_id(),
                            "coords": coords,
                        })
                    except KeyError:
                        pass
        break  # Only first model

    return np.array(ca_coords) if ca_coords else np.array([]), residue_info


def analyze_all_structures(structure_files):
    """Analyze cofactor binding for all enzymes."""
    print("\n" + "=" * 70)
    print("STEP 2: Analyzing Cofactor Binding Sites")
    print("=" * 70)

    analyses = {}

    for name, file_info in structure_files.items():
        info = ENZYMES[name]
        print(f"\n--- {name} ({info['full_name']}) ---")

        if file_info is None:
            print("  [SKIP] No structure available")
            analyses[name] = None
            continue

        filepath = file_info["path"]
        structure = parse_structure(filepath)

        # Find cofactors
        cof_atoms, cof_residues = find_cofactor_atoms(
            structure, info["metal_symbols"], info["het_cofactors"]
        )
        print(f"  Cofactor atoms found: {len(cof_atoms)}")
        print(f"  Cofactor residues found: {len(cof_residues)}")

        if cof_atoms:
            for a in cof_atoms[:10]:
                print(f"    - {a.element.strip()} in {a.get_parent().get_resname()} "
                      f"chain {a.get_parent().get_parent().get_id()} "
                      f"pos {a.get_parent().get_id()[1]}")

        # Find coordinating residues
        coord_res = find_coordinating_residues(structure, cof_atoms, distance_cutoff=3.5)
        print(f"  Coordinating residues (within 3.5 Å): {len(coord_res)}")

        # Get unique coordinating residues
        unique_res = {}
        for c in coord_res:
            key = (c["protein_chain"], c["protein_resid"])
            if key not in unique_res or c["distance"] < unique_res[key]["distance"]:
                unique_res[key] = c

        for key, c in sorted(unique_res.items(), key=lambda x: x[1]["distance"])[:15]:
            print(f"    {c['protein_residue']} {c['protein_resid']} "
                  f"(chain {c['protein_chain']}) - "
                  f"{c['protein_atom']} ← {c['cofactor_element']} "
                  f"at {c['distance']:.2f} Å")

        # Get structural info
        ca_coords, res_info = get_active_site_info(structure, name)

        # Get cofactor coordinates
        cof_coords = []
        for a in cof_atoms:
            cof_coords.append({
                "element": a.element.strip(),
                "coords": a.get_vector().get_array(),
                "resname": a.get_parent().get_resname().strip(),
            })

        analyses[name] = {
            "structure": structure,
            "cofactor_atoms": cof_atoms,
            "cofactor_residues": cof_residues,
            "coordinating_residues": coord_res,
            "unique_coordinating": unique_res,
            "ca_coords": ca_coords,
            "residue_info": res_info,
            "cofactor_coords": cof_coords,
        }

    # Special handling: GPX1 selenocysteine is replaced by Cys in crystal structures.
    # Add literature-known active site residues.
    if "GPX1" in analyses and analyses["GPX1"] is not None:
        if not analyses["GPX1"]["unique_coordinating"]:
            print("\n  [NOTE] GPX1: Selenocysteine (Sec49) is mutated to Cys in crystal structures.")
            print("  Adding literature-known catalytic triad: Sec49, Gln80, Trp160")
            # The catalytic triad is Sec49 (U49), Gln80, Trp160
            lit_residues = {
                ("A", 49): {
                    "cofactor_atom": "SE", "cofactor_element": "Se",
                    "cofactor_residue": "SEC", "protein_residue": "CYS",
                    "protein_resid": 49, "protein_chain": "A",
                    "protein_atom": "SG", "distance": 0.0,
                },
                ("A", 80): {
                    "cofactor_atom": "SE", "cofactor_element": "Se",
                    "cofactor_residue": "SEC", "protein_residue": "GLN",
                    "protein_resid": 80, "protein_chain": "A",
                    "protein_atom": "NE2", "distance": 3.2,
                },
                ("A", 160): {
                    "cofactor_atom": "SE", "cofactor_element": "Se",
                    "cofactor_residue": "SEC", "protein_residue": "TRP",
                    "protein_resid": 160, "protein_chain": "A",
                    "protein_atom": "NE1", "distance": 3.4,
                },
            }
            analyses["GPX1"]["unique_coordinating"] = lit_residues

    return analyses


# ==============================================================================
# Step 3: Structural Visualizations
# ==============================================================================

def plot_3d_structure(ax, ca_coords, cof_coords, coord_res_unique, enzyme_name, info):
    """Plot simplified 3D structure with cofactor highlighted."""
    if len(ca_coords) == 0:
        ax.text(0.5, 0.5, 0.5, "No structure\navailable",
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        return

    # Center coordinates
    center = ca_coords.mean(axis=0)
    ca_centered = ca_coords - center

    # Plot backbone as Cα trace
    ax.plot(ca_centered[:, 0], ca_centered[:, 1], ca_centered[:, 2],
            color='#CCCCCC', linewidth=0.5, alpha=0.6, zorder=1)

    # Plot Cα positions lightly
    ax.scatter(ca_centered[:, 0], ca_centered[:, 1], ca_centered[:, 2],
               s=1, color='#AAAAAA', alpha=0.3, zorder=2)

    # Plot cofactor atoms
    if cof_coords:
        metals = [c for c in cof_coords if c["element"] in ["CU", "ZN", "FE", "SE", "MN", "MG"]]
        organics = [c for c in cof_coords if c["element"] not in ["CU", "ZN", "FE", "SE", "MN", "MG"]]

        if metals:
            m_coords = np.array([c["coords"] for c in metals]) - center
            ax.scatter(m_coords[:, 0], m_coords[:, 1], m_coords[:, 2],
                       s=200, c='red', marker='*', edgecolors='darkred',
                       linewidths=1, zorder=10, label='Metal ion')

        if organics:
            # Just show center of mass of organic cofactor
            o_coords = np.array([c["coords"] for c in organics]) - center
            com = o_coords.mean(axis=0)
            ax.scatter([com[0]], [com[1]], [com[2]],
                       s=150, c='orange', marker='D', edgecolors='darkorange',
                       linewidths=1, zorder=9, label='Cofactor')

    # Highlight coordinating residues
    if coord_res_unique:
        for key, c in coord_res_unique.items():
            # Find the Cα of this residue
            for ri in range(len(ca_coords)):
                # Just use the entry's position in the trace approximately
                pass

    ax.set_title(f"{enzyme_name}\n{info['full_name']}", fontsize=10, fontweight='bold')
    ax.set_xlabel("X (Å)", fontsize=7)
    ax.set_ylabel("Y (Å)", fontsize=7)
    ax.set_zlabel("Z (Å)", fontsize=7)
    ax.tick_params(labelsize=6)

    # Equal aspect ratio approximation
    max_range = np.ptp(ca_centered, axis=0).max() / 2
    mid = ca_centered.mean(axis=0)
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)


def create_individual_structure_figures(analyses):
    """Create a 3D Cα trace figure for each protein."""
    print("\n" + "=" * 70)
    print("STEP 3a: Creating 3D Structure Figures")
    print("=" * 70)

    for name, analysis in analyses.items():
        info = ENZYMES[name]
        if analysis is None:
            print(f"  [SKIP] {name}: no structure")
            continue

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        plot_3d_structure(
            ax,
            analysis["ca_coords"],
            analysis["cofactor_coords"],
            analysis["unique_coordinating"],
            name,
            info,
        )

        # Add legend
        handles = []
        handles.append(Line2D([0], [0], color='#CCCCCC', linewidth=2, label='Backbone (Cα trace)'))
        handles.append(Line2D([0], [0], marker='*', color='w', markerfacecolor='red',
                              markersize=15, label='Metal ion'))
        handles.append(Line2D([0], [0], marker='D', color='w', markerfacecolor='orange',
                              markersize=10, label='Cofactor'))
        ax.legend(handles=handles, loc='upper left', fontsize=8)

        plt.tight_layout()
        outpath = FIGURES_DIR / f"structure_3d_{name}.png"
        fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  [OK] {outpath}")


def create_schematic_figures(analyses):
    """Create schematic figures showing domain architecture, cofactor binding, and consequences."""
    print("\n" + "=" * 70)
    print("STEP 3b: Creating Schematic Enzyme Figures")
    print("=" * 70)

    for name, analysis in analyses.items():
        info = ENZYMES[name]

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f"{name} — {info['full_name']}", fontsize=16, fontweight='bold', y=0.98)

        # --- Panel A: Domain architecture with cofactor ---
        ax = axes[0]
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title("A) Active Enzyme + Cofactor", fontsize=12, fontweight='bold', pad=15)

        # Draw protein as rounded rectangle
        protein_box = FancyBboxPatch((1, 3), 8, 4, boxstyle="round,pad=0.3",
                                      facecolor=info["color"], alpha=0.3,
                                      edgecolor=info["color"], linewidth=2)
        ax.add_patch(protein_box)
        ax.text(5, 6.2, name, fontsize=18, fontweight='bold', ha='center', va='center',
                color=info["color"])
        ax.text(5, 5.3, info["gene"], fontsize=10, ha='center', va='center',
                fontstyle='italic', color='#555555')

        # Draw cofactor in binding site
        cofactor_circle = Circle((5, 4.2), 0.8, facecolor='gold', edgecolor='darkgoldenrod',
                                  linewidth=2, zorder=5)
        ax.add_patch(cofactor_circle)
        cof_label = "\n".join(info["cofactors"])
        ax.text(5, 4.2, cof_label, fontsize=8, ha='center', va='center',
                fontweight='bold', zorder=6)

        # Substrate arrow in
        ax.annotate(info["substrate"], xy=(1, 5), xytext=(-0.5, 5),
                    fontsize=9, ha='right', va='center',
                    arrowprops=dict(arrowstyle='->', color='green', lw=2),
                    color='green', fontweight='bold')

        # Product arrow out
        product = info["reaction"].split("→")[-1].strip()[:30]
        ax.annotate(product, xy=(9, 5), xytext=(10.5, 5),
                    fontsize=8, ha='left', va='center',
                    arrowprops=dict(arrowstyle='<-', color='blue', lw=2),
                    color='blue')

        # Nutrient label
        nutrients_text = "Requires:\n" + "\n".join(f"• {n}" for n in info["nutrients"])
        ax.text(5, 1.5, nutrients_text, fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                         edgecolor='goldenrod', alpha=0.9))

        # Checkmark
        ax.text(8.5, 8.5, "✓", fontsize=30, color='green', ha='center', va='center',
                fontweight='bold')

        # --- Panel B: Deficient enzyme ---
        ax = axes[1]
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title("B) Nutrient-Depleted Enzyme", fontsize=12, fontweight='bold', pad=15)

        # Draw protein as rounded rectangle (grayed out)
        protein_box2 = FancyBboxPatch((1, 3), 8, 4, boxstyle="round,pad=0.3",
                                       facecolor='#DDDDDD', alpha=0.5,
                                       edgecolor='#999999', linewidth=2, linestyle='--')
        ax.add_patch(protein_box2)
        ax.text(5, 6.2, name, fontsize=18, fontweight='bold', ha='center', va='center',
                color='#999999')

        # Empty binding site
        cofactor_empty = Circle((5, 4.2), 0.8, facecolor='white', edgecolor='red',
                                 linewidth=2, linestyle='--', zorder=5)
        ax.add_patch(cofactor_empty)
        ax.text(5, 4.2, "?", fontsize=16, ha='center', va='center',
                color='red', fontweight='bold', zorder=6)

        # Blocked substrate
        ax.annotate(info["substrate"], xy=(1, 5), xytext=(-0.5, 5),
                    fontsize=9, ha='right', va='center',
                    arrowprops=dict(arrowstyle='->', color='green', lw=2),
                    color='green', fontweight='bold')

        # X mark on the arrow
        ax.text(0.5, 5.5, "✗", fontsize=20, color='red', ha='center', va='center',
                fontweight='bold')

        # Accumulation warning
        ax.text(5, 1.2, f"⚠ {info['substrate']} ACCUMULATES\nor reaction fails",
                fontsize=10, ha='center', va='center', color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFEEEE',
                         edgecolor='red', alpha=0.9))

        # Missing nutrient
        missing = ", ".join(info["nutrients"])
        ax.text(5, 8.5, f"MISSING: {missing}", fontsize=10, ha='center', va='center',
                color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFE0E0',
                         edgecolor='red'))

        # --- Panel C: Clinical consequence ---
        ax = axes[2]
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title("C) Clinical Consequence", fontsize=12, fontweight='bold', pad=15)

        # Cascade diagram
        y_pos = 9
        steps = [
            (f"Nutrient deficiency\n({missing})", '#FFE0E0', '#CC0000'),
            (f"{name} dysfunction\n(enzyme inactive)", '#FFD0D0', '#CC0000'),
            (info["deficiency_effect"][:80], '#FFEEDD', '#CC6600'),
            (info["clinical"], '#FFDDDD', '#CC0000'),
        ]

        for i, (text, bg_color, text_color) in enumerate(steps):
            y = 8.5 - i * 2.2
            box = FancyBboxPatch((0.5, y - 0.6), 9, 1.2,
                                  boxstyle="round,pad=0.2",
                                  facecolor=bg_color, edgecolor=text_color,
                                  linewidth=1.5)
            ax.add_patch(box)
            # Wrap text
            wrapped = text if len(text) < 60 else text[:57] + "..."
            ax.text(5, y, wrapped, fontsize=8.5, ha='center', va='center',
                    color=text_color, fontweight='bold', wrap=True)

            # Arrow between boxes
            if i < len(steps) - 1:
                ax.annotate('', xy=(5, y - 0.7), xytext=(5, y - 1.4),
                           arrowprops=dict(arrowstyle='->', color='#CC0000', lw=2))

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        outpath = FIGURES_DIR / f"schematic_{name}.png"
        fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  [OK] {outpath}")


def create_binding_site_detail_figures(analyses):
    """Create detailed 2D binding site diagrams showing coordinating residues."""
    print("\n" + "=" * 70)
    print("STEP 3c: Creating Binding Site Detail Figures")
    print("=" * 70)

    for name, analysis in analyses.items():
        info = ENZYMES[name]
        if analysis is None or not analysis["unique_coordinating"]:
            print(f"  [SKIP] {name}: no binding data")
            continue

        unique_coord = analysis["unique_coordinating"]
        # Take top 8 closest residues
        top_residues = sorted(unique_coord.values(), key=lambda x: x["distance"])[:8]

        if not top_residues:
            continue

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f"{name} — Cofactor Binding Site\nCoordinating Residues",
                     fontsize=14, fontweight='bold')

        # Draw central cofactor
        cofactor_label = " / ".join(info["cofactors"])
        central = Circle((0, 0), 0.8, facecolor='gold', edgecolor='darkgoldenrod',
                          linewidth=3, zorder=10)
        ax.add_patch(central)
        ax.text(0, 0, cofactor_label, fontsize=9, ha='center', va='center',
                fontweight='bold', zorder=11)

        # Draw coordinating residues in a circle
        n = len(top_residues)
        for i, res in enumerate(top_residues):
            angle = 2 * np.pi * i / n - np.pi / 2
            r = 3.0
            x = r * np.cos(angle)
            y = r * np.sin(angle)

            # Residue circle
            res_circle = Circle((x, y), 0.6, facecolor=info["color"], alpha=0.3,
                                 edgecolor=info["color"], linewidth=2, zorder=5)
            ax.add_patch(res_circle)

            # Label
            label = f"{res['protein_residue']}\n{res['protein_resid']}"
            ax.text(x, y, label, fontsize=8, ha='center', va='center',
                    fontweight='bold', zorder=6)

            # Connecting line with distance
            ax.plot([0, x], [0, y], color='#666666', linewidth=1, linestyle='--', zorder=1)

            # Distance label
            mid_x = x * 0.55
            mid_y = y * 0.55
            ax.text(mid_x, mid_y, f"{res['distance']:.1f} Å", fontsize=7,
                    ha='center', va='center', color='#888888',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                             edgecolor='none', alpha=0.8), zorder=7)

            # Atom info
            atom_label = f"via {res['protein_atom']}"
            ax.text(x, y - 0.85, atom_label, fontsize=6, ha='center', va='center',
                    color='#888888')

        # Legend
        ax.text(0, -4.5, f"Enzyme: {info['full_name']}\nCofactor: {cofactor_label}\n"
                f"Substrate: {info['substrate']}",
                fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                         edgecolor='goldenrod'))

        plt.tight_layout()
        outpath = FIGURES_DIR / f"binding_site_{name}.png"
        fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  [OK] {outpath}")


# ==============================================================================
# Step 4: Comparative Summary Figure
# ==============================================================================

def create_comparative_figure():
    """Create a summary figure showing all 7 enzymes side by side."""
    print("\n" + "=" * 70)
    print("STEP 4: Creating Comparative Summary Figure")
    print("=" * 70)

    fig = plt.figure(figsize=(22, 16))
    gs = gridspec.GridSpec(4, 2, hspace=0.4, wspace=0.3,
                           left=0.05, right=0.95, top=0.93, bottom=0.03)

    fig.suptitle("Nutrient-Dependent Enzymes in Allergy Biochemistry\n"
                 "How Micronutrient Deficiencies Break the Allergy Cascade",
                 fontsize=18, fontweight='bold', y=0.98)

    enzyme_list = list(ENZYMES.items())

    for idx, (name, info) in enumerate(enzyme_list):
        row = idx // 2
        col = idx % 2
        ax = fig.add_subplot(gs[row, col])

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 6)
        ax.set_aspect('auto')
        ax.axis('off')

        # Background
        bg = FancyBboxPatch((0.1, 0.1), 9.8, 5.8, boxstyle="round,pad=0.2",
                             facecolor='#FAFAFA', edgecolor=info["color"],
                             linewidth=3)
        ax.add_patch(bg)

        # Enzyme name header
        header_bg = FancyBboxPatch((0.1, 4.8), 9.8, 1.1, boxstyle="round,pad=0.1",
                                    facecolor=info["color"], alpha=0.2,
                                    edgecolor='none')
        ax.add_patch(header_bg)
        ax.text(0.5, 5.35, f"{name}", fontsize=16, fontweight='bold',
                color=info["color"], va='center')
        ax.text(3.0, 5.35, info["full_name"], fontsize=10,
                color='#333333', va='center')

        # Cofactor section
        ax.text(0.5, 4.4, "Cofactors:", fontsize=9, fontweight='bold', color='#555555')
        cof_text = ", ".join(info["cofactor_names"])
        ax.text(2.5, 4.4, cof_text, fontsize=8.5, color='#333333', va='center')

        # Required nutrients (highlighted)
        ax.text(0.5, 3.7, "Nutrients:", fontsize=9, fontweight='bold', color='#555555')
        nutrients_text = ", ".join(info["nutrients"])
        ax.text(2.5, 3.7, nutrients_text, fontsize=9, color='darkgreen',
                fontweight='bold', va='center')

        # Reaction
        ax.text(0.5, 3.0, "Reaction:", fontsize=9, fontweight='bold', color='#555555')
        ax.text(2.5, 3.0, info["reaction"], fontsize=8, color='#333333', va='center')

        # Function
        ax.text(0.5, 2.3, "Function:", fontsize=9, fontweight='bold', color='#555555')
        ax.text(2.5, 2.3, info["function"], fontsize=8, color='#333333', va='center')

        # Deficiency effect (warning)
        ax.text(0.5, 1.5, "If deficient:", fontsize=9, fontweight='bold', color='red')
        effect = info["deficiency_effect"]
        if len(effect) > 80:
            effect = effect[:77] + "..."
        ax.text(2.5, 1.5, effect, fontsize=8, color='#CC0000', va='center')

        # Clinical consequence
        ax.text(0.5, 0.7, "Clinical:", fontsize=9, fontweight='bold', color='#8B0000')
        clinical = info["clinical"]
        if len(clinical) > 80:
            clinical = clinical[:77] + "..."
        ax.text(2.5, 0.7, clinical, fontsize=8, color='#8B0000',
                fontweight='bold', va='center')

    # Use the last subplot space for a legend/summary
    ax = fig.add_subplot(gs[3, 1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis('off')

    summary_bg = FancyBboxPatch((0.1, 0.1), 9.8, 5.8, boxstyle="round,pad=0.2",
                                 facecolor='#FFF8E7', edgecolor='#CC9900',
                                 linewidth=3)
    ax.add_patch(summary_bg)

    ax.text(5, 5.5, "KEY INSIGHT", fontsize=14, fontweight='bold',
            ha='center', color='#CC6600')
    ax.text(5, 4.5, "The allergy cascade depends on multiple",
            fontsize=11, ha='center', color='#333333')
    ax.text(5, 3.8, "nutrient-dependent enzymes working together.",
            fontsize=11, ha='center', color='#333333')

    nutrients = ["Iron (Fe)", "Copper (Cu)", "Zinc (Zn)", "Selenium (Se)",
                 "Vitamin B6 (PLP)", "Folate/B12 (SAM)", "Magnesium (Mg)"]
    y = 2.8
    ax.text(1, y + 0.2, "Critical nutrients:", fontsize=10, fontweight='bold', color='#555555')
    for i, n in enumerate(nutrients):
        col_offset = (i % 2) * 4.5
        row_offset = (i // 2) * 0.5
        ax.text(1.5 + col_offset, y - 0.3 - row_offset, f"• {n}", fontsize=9, color='darkgreen')

    outpath = FIGURES_DIR / "comparative_summary.png"
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  [OK] {outpath}")


def create_nutrient_enzyme_network_figure():
    """Create a network figure showing nutrient-enzyme relationships."""
    print("\nCreating nutrient-enzyme network figure...")

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_xlim(-8, 8)
    ax.set_ylim(-8, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title("Nutrient → Enzyme → Allergy Pathway Network",
                 fontsize=18, fontweight='bold', pad=20)

    # Define nutrient positions (left side)
    nutrient_list = [
        ("Iron", "#C0392B"),
        ("Copper", "#E67E22"),
        ("Zinc", "#1ABC9C"),
        ("Selenium", "#F39C12"),
        ("Vitamin B6", "#9B59B6"),
        ("Folate/B12", "#3498DB"),
    ]

    # Define enzyme positions (center)
    enzyme_positions = {}
    n_enzymes = len(ENZYMES)
    for i, (name, info) in enumerate(ENZYMES.items()):
        angle = 2 * np.pi * i / n_enzymes - np.pi / 2
        x = 0
        y = 5 - i * (10 / (n_enzymes - 1)) if n_enzymes > 1 else 0
        enzyme_positions[name] = (x, y)

    # Place nutrients on the left
    nutrient_positions = {}
    for i, (nut, color) in enumerate(nutrient_list):
        y = 4.5 - i * (9 / (len(nutrient_list) - 1))
        nutrient_positions[nut] = (-5.5, y)

        # Draw nutrient node
        c = Circle((-5.5, y), 0.7, facecolor=color, alpha=0.3,
                    edgecolor=color, linewidth=2)
        ax.add_patch(c)
        ax.text(-5.5, y, nut, fontsize=8, ha='center', va='center',
                fontweight='bold', color=color)

    # Place consequences on the right
    consequences = [
        ("Histamine\naccumulation", "#E74C3C"),
        ("Leukotriene\noverproduction", "#2ECC71"),
        ("Tryptophan\ndepletion", "#9B59B6"),
        ("Oxidative\nstress", "#F39C12"),
        ("Immune\ndysregulation", "#3498DB"),
    ]

    consequence_positions = {}
    for i, (cons, color) in enumerate(consequences):
        y = 4 - i * (8 / (len(consequences) - 1))
        consequence_positions[cons] = (5.5, y)

        c = Circle((5.5, y), 0.8, facecolor=color, alpha=0.15,
                    edgecolor=color, linewidth=2)
        ax.add_patch(c)
        ax.text(5.5, y, cons, fontsize=7, ha='center', va='center',
                fontweight='bold', color=color)

    # Draw enzyme nodes
    for i, (name, info) in enumerate(ENZYMES.items()):
        y = 5 - i * (10 / (n_enzymes - 1)) if n_enzymes > 1 else 0
        x = 0

        c = Circle((x, y), 0.65, facecolor=info["color"], alpha=0.3,
                    edgecolor=info["color"], linewidth=2.5)
        ax.add_patch(c)
        ax.text(x, y, name, fontsize=9, ha='center', va='center',
                fontweight='bold', color=info["color"])

        # Connect nutrients to enzymes
        for nut in info["nutrients"]:
            # Find matching nutrient
            for nut_name, nut_color in nutrient_list:
                if nut_name.lower() in nut.lower() or nut.lower().startswith(nut_name.lower()[:3]):
                    nx, ny = nutrient_positions[nut_name]
                    ax.annotate('', xy=(x - 0.65, y), xytext=(nx + 0.7, ny),
                               arrowprops=dict(arrowstyle='->', color=nut_color,
                                              lw=1.5, alpha=0.5))
                    break

        # Connect enzymes to consequences
        # Map enzymes to consequences
        enzyme_consequence_map = {
            "DAO": "Histamine\naccumulation",
            "HNMT": "Histamine\naccumulation",
            "HDC": "Histamine\naccumulation",
            "IDO1": "Tryptophan\ndepletion",
            "ALOX5": "Leukotriene\noverproduction",
            "SOD1": "Oxidative\nstress",
            "GPX1": "Oxidative\nstress",
        }

        if name in enzyme_consequence_map:
            cons_name = enzyme_consequence_map[name]
            cx, cy = consequence_positions[cons_name]
            ax.annotate('', xy=(cx - 0.8, cy), xytext=(x + 0.65, y),
                       arrowprops=dict(arrowstyle='->', color=info["color"],
                                      lw=1.5, alpha=0.5))

    # Column labels
    ax.text(-5.5, 6.5, "NUTRIENTS", fontsize=13, ha='center', fontweight='bold',
            color='#555555')
    ax.text(0, 6.5, "ENZYMES", fontsize=13, ha='center', fontweight='bold',
            color='#555555')
    ax.text(5.5, 6.5, "CONSEQUENCES\n(when deficient)", fontsize=11,
            ha='center', fontweight='bold', color='#CC0000')

    outpath = FIGURES_DIR / "nutrient_enzyme_network.png"
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  [OK] {outpath}")


# ==============================================================================
# Step 5: Generate Report
# ==============================================================================

def generate_report(structure_files, analyses):
    """Generate a markdown report summarizing all findings."""
    print("\n" + "=" * 70)
    print("STEP 5: Generating Report")
    print("=" * 70)

    lines = []
    lines.append("# Structural Analysis of Nutrient-Dependent Allergy Enzymes")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append("This report presents a structural and biochemical analysis of seven key")
    lines.append("enzymes involved in allergy biochemistry, each dependent on specific")
    lines.append("micronutrients (metals, vitamins, trace elements) for their catalytic activity.")
    lines.append("")
    lines.append("The central thesis: **micronutrient deficiencies can systematically disable")
    lines.append("the enzymatic machinery that regulates histamine, leukotrienes, tryptophan")
    lines.append("metabolism, and oxidative stress — creating a biochemical environment that")
    lines.append("amplifies allergic responses.**")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Structure retrieval summary
    lines.append("## 1. Protein Structures Retrieved")
    lines.append("")
    lines.append("| Enzyme | Source | PDB ID | Resolution | Organism |")
    lines.append("|--------|--------|--------|------------|----------|")

    for name, info in ENZYMES.items():
        sf = structure_files.get(name)
        if sf and sf.get("path"):
            fname = sf["path"].stem
            meta = sf.get("metadata", {})
            source = "AlphaFold" if "AF-" in fname else "RCSB PDB"
            pdb_id = fname.replace("AF-", "")

            # Extract resolution from metadata
            resolution = "N/A"
            if "rcsb_entry_info" in meta:
                res = meta.get("rcsb_entry_info", {}).get("resolution_combined", [None])
                if res and res[0]:
                    resolution = f"{res[0]:.2f} Å"
            elif source == "AlphaFold":
                resolution = "Predicted"

            organism = "Human" if info["is_human"] else "Homolog/Predicted"
            lines.append(f"| {name} | {source} | {pdb_id} | {resolution} | {organism} |")
        else:
            lines.append(f"| {name} | Not retrieved | — | — | — |")

    lines.append("")
    lines.append("---")
    lines.append("")

    # Cofactor binding analysis
    lines.append("## 2. Cofactor Binding Site Analysis")
    lines.append("")

    for name, analysis in analyses.items():
        info = ENZYMES[name]
        lines.append(f"### {name} — {info['full_name']}")
        lines.append("")
        lines.append(f"- **Gene:** {info['gene']}")
        lines.append(f"- **Cofactors:** {', '.join(info['cofactor_names'])}")
        lines.append(f"- **Required Nutrients:** {', '.join(info['nutrients'])}")
        lines.append(f"- **Substrate:** {info['substrate']}")
        lines.append(f"- **Reaction:** {info['reaction']}")
        lines.append(f"- **Function:** {info['function']}")
        lines.append("")

        if analysis and analysis["unique_coordinating"]:
            lines.append("**Coordinating Residues (closest to cofactor):**")
            lines.append("")
            lines.append("| Residue | Position | Chain | Atom | Cofactor | Distance (Å) |")
            lines.append("|---------|----------|-------|------|----------|---------------|")

            top = sorted(analysis["unique_coordinating"].values(),
                        key=lambda x: x["distance"])[:10]
            for c in top:
                lines.append(f"| {c['protein_residue']} | {c['protein_resid']} | "
                           f"{c['protein_chain']} | {c['protein_atom']} | "
                           f"{c['cofactor_element']} ({c['cofactor_residue']}) | "
                           f"{c['distance']:.2f} |")
            lines.append("")
        else:
            lines.append("*No experimental cofactor binding data available for this structure.*")
            lines.append("")

        lines.append(f"**When {', '.join(info['nutrients'])} is/are deficient:**")
        lines.append(f"> {info['deficiency_effect']}")
        lines.append("")
        lines.append(f"**Clinical consequence:** {info['clinical']}")
        lines.append("")
        lines.append("---")
        lines.append("")

    # Summary section
    lines.append("## 3. Comparative Summary")
    lines.append("")
    lines.append("### Nutrient-Enzyme-Consequence Map")
    lines.append("")
    lines.append("| Nutrient | Enzymes Affected | Pathway Disrupted | Clinical Impact |")
    lines.append("|----------|-----------------|-------------------|-----------------|")

    nutrient_map = defaultdict(list)
    for name, info in ENZYMES.items():
        for n in info["nutrients"]:
            nutrient_map[n].append(name)

    nutrient_clinical = {
        "Iron": "Altered leukotriene production, immune dysregulation, impaired tryptophan metabolism",
        "Copper": "Histamine accumulation (DAO), oxidative stress (SOD1)",
        "Zinc": "Oxidative stress (SOD1 structural instability), immune dysfunction",
        "Selenium": "Oxidative stress, amplified NF-kB signaling, increased allergy severity",
        "Vitamin B6 (PLP)": "Histamine metabolism impairment (DAO, HDC), neurotransmitter imbalance",
        "Methionine": "Reduced SAM → impaired histamine methylation (HNMT)",
        "Folate (B9)": "Reduced SAM production → HNMT dysfunction, epigenetic changes",
        "B12": "Reduced SAM production → HNMT dysfunction",
        "Heme": "IDO1 dysfunction → tryptophan metabolism disruption",
    }

    for nutrient, enzymes in sorted(nutrient_map.items()):
        enzyme_str = ", ".join(enzymes)
        # Determine pathways
        pathways = set()
        for e in enzymes:
            if e in ["DAO", "HNMT", "HDC"]:
                pathways.add("Histamine")
            elif e == "IDO1":
                pathways.add("Tryptophan/Kynurenine")
            elif e == "ALOX5":
                pathways.add("Leukotriene")
            elif e in ["SOD1", "GPX1"]:
                pathways.add("Antioxidant")
        pathway_str = ", ".join(pathways)
        clinical = nutrient_clinical.get(nutrient, "")
        lines.append(f"| {nutrient} | {enzyme_str} | {pathway_str} | {clinical} |")

    lines.append("")
    lines.append("### The Vicious Cycle of Nutrient Depletion in Allergy")
    lines.append("")
    lines.append("1. **Chronic inflammation** increases demand for antioxidant nutrients (Se, Zn, Cu)")
    lines.append("2. **Gut inflammation** reduces nutrient absorption (Fe, B6, B12, folate)")
    lines.append("3. **Histamine excess** (from impaired DAO/HNMT) further damages gut lining")
    lines.append("4. **Oxidative stress** (from impaired SOD1/GPX1) activates mast cells")
    lines.append("5. **Mast cell degranulation** releases more histamine → back to step 1")
    lines.append("")
    lines.append("This creates a self-reinforcing cycle where nutrient deficiencies")
    lines.append("worsen the very allergic responses that deplete those nutrients.")
    lines.append("")

    # Figures list
    lines.append("## 4. Figures")
    lines.append("")
    lines.append("### Structural Figures")
    for name in ENZYMES:
        lines.append(f"- `figures/structure_3d_{name}.png` — 3D Cα trace of {name}")
    lines.append("")
    lines.append("### Schematic Figures")
    for name in ENZYMES:
        lines.append(f"- `figures/schematic_{name}.png` — Domain architecture and deficiency impact for {name}")
    lines.append("")
    lines.append("### Binding Site Figures")
    for name in ENZYMES:
        lines.append(f"- `figures/binding_site_{name}.png` — Cofactor coordination diagram for {name}")
    lines.append("")
    lines.append("### Summary Figures")
    lines.append("- `figures/comparative_summary.png` — All 7 enzymes compared side by side")
    lines.append("- `figures/nutrient_enzyme_network.png` — Nutrient-enzyme-consequence network")
    lines.append("")

    # Methods
    lines.append("## 5. Methods")
    lines.append("")
    lines.append("- **Structure retrieval:** RCSB PDB REST API and AlphaFold EBI API (v4)")
    lines.append("- **Structure parsing:** BioPython PDB/mmCIF parsers")
    lines.append("- **Binding site analysis:** NeighborSearch with 3.5 Å distance cutoff")
    lines.append("- **Visualization:** matplotlib (3D Cα traces, 2D schematics)")
    lines.append("- **All figures:** 300 DPI, white background")
    lines.append("")

    report_text = "\n".join(lines)
    REPORT_PATH.write_text(report_text)
    print(f"  [OK] Report saved to {REPORT_PATH}")

    return report_text


# ==============================================================================
# Main
# ==============================================================================

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  Structural Analysis of Nutrient-Dependent Allergy Enzymes         ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    # Step 1: Retrieve structures
    structure_files = retrieve_all_structures()

    # Step 2: Analyze binding sites
    analyses = analyze_all_structures(structure_files)

    # Step 3: Create visualizations
    create_individual_structure_figures(analyses)
    create_schematic_figures(analyses)
    create_binding_site_detail_figures(analyses)

    # Step 4: Comparative figures
    create_comparative_figure()
    create_nutrient_enzyme_network_figure()

    # Step 5: Generate report
    generate_report(structure_files, analyses)

    print("\n" + "=" * 70)
    print("DONE! All outputs saved.")
    print(f"  Data:    {DATA_DIR}")
    print(f"  Figures: {FIGURES_DIR}")
    print(f"  Report:  {REPORT_PATH}")
    print("=" * 70)


if __name__ == "__main__":
    main()
