# PDBSiteScan-PLACER Pipeline

A Python tool that integrates [PDBSiteScan](https://wwwmgs.bionet.nsc.ru/cgi-bin/mgs/fastprot/pdbsitescan/) binding site predictions with [PLACER](https://github.com/baker-laboratory/PLACER) for protein-ligand docking analysis.

## Overview

This pipeline:
1. Parses PDBSiteScan output to identify potential binding sites
2. Downloads reference PDB structures from RCSB
3. Extracts and superimposes ligands onto target protein
4. Runs PLACER for ligand pose optimization
5. Outputs confidence metrics for each predicted binding site

## Installation

### Prerequisites

- Python 3.8+
- [PLACER](https://github.com/baker-laboratory/PLACER) installed and configured
- NumPy

### Setup

```bash
git clone https://github.com/your-username/pdbsitescan-placer.git
cd pdbsitescan-placer

## Usage

### Basic Usage

```bash
python pdbsitescan_placer.py -s <sitescan_output.txt> -p <protein.pdb> -o <output_dir>
```

### Arguments

| Argument | Short | Required | Default | Description |
|----------|-------|----------|---------|-------------|
| `--sitescan` | `-s` | Yes | - | PDBSiteScan output file (txt format) |
| `--pdb` | `-p` | Yes | - | Input PDB file to analyze |
| `--output` | `-o` | No | `./placer_output` | Output directory for results |
| `--nsamples` | `-n` | No | `1` | Number of PLACER samples (recommended: 50-100) |
| `--placer` | - | No | `/Users/venzel/models/PLACER-main` | Path to PLACER installation |

### Example

```bash
# Run with default settings
python pdbsitescan_placer.py -s 1rm8_pdbsitescan_test.txt -p 1rm8.pdb -o ./results

# Run with more samples for better accuracy
python pdbsitescan_placer.py -s sitescan.txt -p protein.pdb -o ./results -n 50

# Specify custom PLACER path
python pdbsitescan_placer.py -s sitescan.txt -p protein.pdb --placer /path/to/PLACER
```

## Output

```
output_dir/
├── METRICS_EXPLANATION.txt      # Guide to interpreting metrics
├── summary_results.csv          # Combined results table
├── cleaned_protein.pdb          # Input PDB with non-protein atoms removed
├── site_1_<PDBID>_<LIGAND>/
│   ├── <pdbid>.pdb              # Downloaded reference structure
│   ├── ligand_<name>.pdb        # Extracted and transformed ligand
│   ├── protein_ligand_combined.pdb
│   ├── *_model.pdb              # PLACER output structure
│   └── *.csv                    # PLACER metrics
├── site_2_.../
└── ...
```
## Workflow

```
┌─────────────────────┐
│  PDBSiteScan Output │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐     ┌─────────────────────┐
│   Parse binding     │     │   Clean input PDB   │
│   site information  │     │   (remove HETATM)   │
└──────────┬──────────┘     └──────────┬──────────┘
           │                           │
           ▼                           │
┌─────────────────────┐                │
│  Download reference │                │
│  PDB from RCSB      │                │
└──────────┬──────────┘                │
           │                           │
           ▼                           │
┌─────────────────────┐                │
│  Extract ligand &   │                │
│  superimpose sites  │◄───────────────┘
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│   Combine protein   │
│   and ligand        │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│    Run PLACER       │
│    optimization     │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  Output metrics &   │
│  optimized models   │
└─────────────────────┘
```

## PDBSiteScan Input Format

The input file should be a text export from PDBSiteScan containing entries like:

```
ID 1OFYCC3 Include in structure alignment
PDBID 1OFY
SITE_TYPE COMPOUND_BINDING
SITE_DESCR GOL BINDING SITE FOR CHAIN B
LIGAND GLYCEROL;
LIGAND_FORMUL C3 H8 O3;
Max Dist    1.53
RMSD        1.12
Site chain in the protein:     A    A    A
Site positions in the protein: 153  219  221
Site residues in the protein:  R    D    H
Site chain in the PDBSite:     A    B    B
Site positions in the PDBSite: 256  77   81
Site residues in the PDBSite:  R    D    H
END
```

## References

- **PLACER**: [PLACER: Protein-Ligand Atomistic Conformational Ensemble Resolver](https://www.biorxiv.org/content/10.1101/2024.09.25.614868v1)
- **PDBSiteScan web-server**: https://wwwmgs.bionet.nsc.ru/mgs/systems/fastprot/pdbsitescan.html
- **PDBSiteScan paper**: [PDBSiteScan: a program for searching for active, binding and posttranslational modification sites in the 3D structures of proteins](https://doi.org/10.1093/nar/gkh439)

## License

MIT License

---

*This tool was co-developed with [Claude](https://claude.ai), Anthropic's AI assistant.*
