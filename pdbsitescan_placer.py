#!/usr/bin/env python3
"""
PDBSiteScan to PLACER Pipeline

This script:
1. Parses PDBSiteScan output to get binding site information
2. Cleans input PDB (removes non-protein atoms)
3. Downloads reference PDBs from RCSB and extracts ligands
4. Superimposes binding sites and transforms ligand coordinates
5. Runs PLACER for each ligand with the binding site residues
6. Outputs confidence metrics table
"""

import argparse
import os
import re
import sys
import subprocess
import urllib.request
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import tempfile
import shutil
import numpy as np

# Default PLACER path
DEFAULT_PLACER_PATH = "/Users/venzel/models/PLACER-main"


@dataclass
class SiteEntry:
    """Represents a single site entry from PDBSiteScan output."""
    id: str
    pdb_id: str
    site_type: str
    site_descr: str
    ligand_name: str
    ligand_formula: str
    max_dist: float
    rmsd: float
    protein_chains: List[str]
    protein_positions: List[int]
    protein_residues: List[str]
    pdbsite_chains: List[str]
    pdbsite_positions: List[int]
    pdbsite_residues: List[str]


def parse_pdbsitescan(filepath: str) -> List[SiteEntry]:
    """Parse PDBSiteScan output file and return list of site entries."""
    entries = []

    with open(filepath, 'r') as f:
        content = f.read()

    # Split by END markers
    blocks = re.split(r'\bEND\b', content)

    for block in blocks:
        block = block.strip()
        if not block:
            continue

        lines = {}
        for line in block.split('\n'):
            # Remove line numbers and arrows if present
            line = re.sub(r'^\s*\d+→', '', line).strip()
            if not line:
                continue

            # Parse key-value pairs
            if line.startswith('ID '):
                lines['ID'] = line[3:].split()[0]
            elif line.startswith('PDBID '):
                lines['PDBID'] = line[6:].strip()
            elif line.startswith('SITE_TYPE '):
                lines['SITE_TYPE'] = line[10:].strip()
            elif line.startswith('SITE_DESCR '):
                lines['SITE_DESCR'] = line[11:].strip()
            elif line.startswith('LIGAND ') and 'LIGAND_FORMUL' not in line:
                lines['LIGAND'] = line[7:].strip().rstrip(';')
            elif line.startswith('LIGAND_FORMUL '):
                lines['LIGAND_FORMUL'] = line[14:].strip().rstrip(';')
            elif line.startswith('Max Dist'):
                try:
                    lines['MAX_DIST'] = float(line.split()[-1])
                except:
                    lines['MAX_DIST'] = 0.0
            elif line.startswith('RMSD'):
                try:
                    lines['RMSD'] = float(line.split()[-1])
                except:
                    lines['RMSD'] = 0.0
            elif line.startswith('Site chain in the protein:'):
                lines['PROTEIN_CHAINS'] = line.replace('Site chain in the protein:', '').split()
            elif line.startswith('Site positions in the protein:'):
                positions = line.replace('Site positions in the protein:', '').split()
                lines['PROTEIN_POSITIONS'] = [int(p) for p in positions]
            elif line.startswith('Site residues in the protein:'):
                lines['PROTEIN_RESIDUES'] = line.replace('Site residues in the protein:', '').split()
            elif line.startswith('Site chain in the PDBSite:'):
                lines['PDBSITE_CHAINS'] = line.replace('Site chain in the PDBSite:', '').split()
            elif line.startswith('Site positions in the PDBSite:'):
                positions = line.replace('Site positions in the PDBSite:', '').split()
                lines['PDBSITE_POSITIONS'] = [int(p) for p in positions]
            elif line.startswith('Site residues in the PDBSite:'):
                lines['PDBSITE_RESIDUES'] = line.replace('Site residues in the PDBSite:', '').split()

        # Create entry if we have minimum required fields
        if 'PDBID' in lines and 'PROTEIN_POSITIONS' in lines:
            entry = SiteEntry(
                id=lines.get('ID', ''),
                pdb_id=lines.get('PDBID', ''),
                site_type=lines.get('SITE_TYPE', ''),
                site_descr=lines.get('SITE_DESCR', ''),
                ligand_name=lines.get('LIGAND', ''),
                ligand_formula=lines.get('LIGAND_FORMUL', ''),
                max_dist=lines.get('MAX_DIST', 0.0),
                rmsd=lines.get('RMSD', 0.0),
                protein_chains=lines.get('PROTEIN_CHAINS', []),
                protein_positions=lines.get('PROTEIN_POSITIONS', []),
                protein_residues=lines.get('PROTEIN_RESIDUES', []),
                pdbsite_chains=lines.get('PDBSITE_CHAINS', []),
                pdbsite_positions=lines.get('PDBSITE_POSITIONS', []),
                pdbsite_residues=lines.get('PDBSITE_RESIDUES', [])
            )
            entries.append(entry)

    return entries


def clean_pdb(input_pdb: str, output_pdb: str) -> None:
    """Clean PDB file by removing non-protein atoms (HETATM, waters, etc.)."""
    protein_lines = []

    with open(input_pdb, 'r') as f:
        for line in f:
            # Keep only ATOM records (protein atoms) and essential headers
            if line.startswith('ATOM'):
                protein_lines.append(line)
            elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE',
                                   'EXPDTA', 'AUTHOR', 'CRYST1', 'ORIGX',
                                   'SCALE', 'MODEL', 'ENDMDL', 'TER')):
                protein_lines.append(line)

    # Add END if not present
    if protein_lines and not protein_lines[-1].startswith('END'):
        protein_lines.append('END\n')

    with open(output_pdb, 'w') as f:
        f.writelines(protein_lines)


def get_ca_coords(pdb_path: str, residues: List[Tuple[str, int]]) -> Tuple[np.ndarray, List[Tuple[str, int]]]:
    """
    Extract CA coordinates for specified residues.
    Returns coordinates array and list of found residues.
    """
    coords = []
    found_residues = []

    # Build set of target residues
    target_set = set((chain, resnum) for chain, resnum in residues)

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name != 'CA':
                    continue

                chain = line[21].strip()
                try:
                    resnum = int(line[22:26].strip())
                except:
                    continue

                if (chain, resnum) in target_set:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    found_residues.append((chain, resnum))

    return np.array(coords), found_residues


def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation matrix using Kabsch algorithm.
    P and Q are Nx3 matrices of corresponding points.
    Returns rotation matrix R and translation vector t such that Q ≈ P @ R + t
    """
    # Center the points
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation
    R = Vt.T @ U.T

    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Compute translation
    t = centroid_Q - centroid_P @ R

    return R, t


def transform_ligand_pdb(ligand_pdb: str, output_pdb: str, R: np.ndarray, t: np.ndarray,
                         new_chain: str = 'L') -> None:
    """
    Transform ligand coordinates using rotation R and translation t.
    Also changes chain ID to avoid conflicts with protein chains.
    """
    transformed_lines = []
    resnum_counter = 1
    seen_residues = set()

    with open(ligand_pdb, 'r') as f:
        for line in f:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Transform
                coord = np.array([x, y, z])
                new_coord = coord @ R + t

                # Get current residue info
                old_chain = line[21]
                old_resnum = line[22:26].strip()
                resname = line[17:20].strip()
                res_key = f"{old_chain}_{old_resnum}_{resname}"

                # Assign new residue number (keep same for same residue)
                if res_key not in seen_residues:
                    seen_residues.add(res_key)

                # Rebuild line with new coordinates and new chain
                new_line = (line[:21] + new_chain + f"{resnum_counter:4d}" +
                           line[26:30] + f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}" + line[54:])
                transformed_lines.append(new_line)
            else:
                transformed_lines.append(line)

    with open(output_pdb, 'w') as f:
        f.writelines(transformed_lines)


def download_pdb(pdb_id: str, output_path: str) -> bool:
    """Download PDB file from RCSB."""
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        # Try mmCIF format
        url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
        try:
            cif_path = output_path.replace('.pdb', '.cif')
            urllib.request.urlretrieve(url_cif, cif_path)
            return True
        except Exception as e2:
            print(f"Error downloading CIF {pdb_id}: {e2}")
            return False


def extract_ligand_from_pdb(pdb_path: str, site_descr: str, output_ligand_path: str) -> Optional[str]:
    """
    Extract ligand from PDB file based on site description.
    Returns the 3-letter ligand code if successful.
    """
    # Parse ligand code from site description (e.g., "GOL BINDING SITE FOR CHAIN B")
    ligand_code = site_descr.split()[0] if site_descr else None

    if not ligand_code:
        print(f"Could not determine ligand code from: {site_descr}")
        return None

    ligand_lines = []

    # Determine file type
    is_cif = pdb_path.endswith('.cif')

    if is_cif:
        # For mmCIF files, extract HETATM records for the ligand
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') and ligand_code in line:
                    ligand_lines.append(line)
    else:
        # For PDB files
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    # Check if this is our ligand (residue name at positions 17-20)
                    resname = line[17:20].strip()
                    if resname == ligand_code:
                        ligand_lines.append(line)

    if ligand_lines:
        with open(output_ligand_path, 'w') as f:
            f.writelines(ligand_lines)
            f.write('END\n')
        return ligand_code
    else:
        print(f"Ligand {ligand_code} not found in {pdb_path}")
        return None


def combine_protein_and_ligand(protein_pdb: str, ligand_pdb: str, output_pdb: str) -> None:
    """Combine cleaned protein PDB with ligand PDB."""
    with open(output_pdb, 'w') as out:
        # Write protein atoms
        with open(protein_pdb, 'r') as prot:
            for line in prot:
                if line.startswith('END'):
                    continue
                out.write(line)

        # Write ligand atoms
        with open(ligand_pdb, 'r') as lig:
            for line in lig:
                out.write(line)


def run_placer(input_pdb: str, ligand_spec: str, nsamples: int,
               output_dir: str, suffix: str, placer_path: str) -> Optional[str]:
    """
    Run PLACER on the input PDB with specified ligand to predict.
    ligand_spec should be in format: chain-name3-resno (e.g., L-GOL-1)
    Returns path to output CSV if successful.
    """
    if not ligand_spec:
        print("No ligand specification provided")
        return None

    # Convert to absolute paths
    input_pdb_abs = os.path.abspath(input_pdb)
    output_dir_abs = os.path.abspath(output_dir)

    cmd = [
        sys.executable,
        os.path.join(placer_path, "run_PLACER.py"),
        "--ifile", input_pdb_abs,
        "--odir", output_dir_abs,
        "--nsamples", str(nsamples),
        "--predict_ligand", ligand_spec,
        "--suffix", suffix,
        "--rerank", "prmsd"
    ]

    print(f"Running PLACER: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=placer_path)
        print(result.stdout)
        if result.returncode != 0:
            print(f"PLACER error (stderr): {result.stderr}")
            return None
        if result.stderr:
            print(f"PLACER stderr: {result.stderr}")

        # Find output CSV
        base_name = os.path.splitext(os.path.basename(input_pdb))[0]
        csv_pattern = f"{base_name}_{suffix}.csv"
        csv_path = os.path.join(output_dir_abs, csv_pattern)

        if os.path.exists(csv_path):
            return csv_path

        # Try to find any CSV in output dir
        for f in os.listdir(output_dir_abs):
            if f.endswith('.csv') and suffix in f:
                return os.path.join(output_dir_abs, f)

        return None
    except Exception as e:
        print(f"Error running PLACER: {e}")
        return None


def parse_placer_csv(csv_path: str) -> Dict:
    """Parse PLACER output CSV and return confidence metrics."""
    metrics = {}

    try:
        with open(csv_path, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            return metrics

        headers = lines[0].strip().split(',')

        # Get metrics from first (best) model
        values = lines[1].strip().split(',')

        for h, v in zip(headers, values):
            try:
                metrics[h] = float(v)
            except:
                metrics[h] = v
    except Exception as e:
        print(f"Error parsing CSV {csv_path}: {e}")

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description='Run PLACER analysis based on PDBSiteScan results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sitescan', '-s', required=True,
                        help='PDBSiteScan output file (txt format)')
    parser.add_argument('--pdb', '-p', required=True,
                        help='Input PDB file to analyze')
    parser.add_argument('--output', '-o', default='./placer_output',
                        help='Output directory for results')
    parser.add_argument('--nsamples', '-n', type=int, default=1,
                        help='Number of samples for PLACER')
    parser.add_argument('--placer', type=str, default=DEFAULT_PLACER_PATH,
                        help='Path to PLACER installation directory')

    args = parser.parse_args()

    # Validate PLACER path
    if not os.path.exists(os.path.join(args.placer, "run_PLACER.py")):
        print(f"Error: PLACER not found at {args.placer}")
        print("Please specify correct path with --placer argument")
        return

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Parse PDBSiteScan output
    print(f"Parsing PDBSiteScan output: {args.sitescan}")
    site_entries = parse_pdbsitescan(args.sitescan)
    print(f"Found {len(site_entries)} site entries")

    if not site_entries:
        print("No site entries found. Exiting.")
        return

    # Clean input PDB
    cleaned_pdb = os.path.join(args.output, 'cleaned_protein.pdb')
    print(f"Cleaning PDB: {args.pdb}")
    clean_pdb(args.pdb, cleaned_pdb)

    # Process each site entry
    results = []

    for i, entry in enumerate(site_entries):
        print(f"\n{'='*60}")
        print(f"Processing site {i+1}/{len(site_entries)}: {entry.id}")
        print(f"  PDB: {entry.pdb_id}, Ligand: {entry.ligand_name}")
        print(f"  Site residues in protein: {list(zip(entry.protein_chains, entry.protein_positions, entry.protein_residues))}")
        print(f"  Site residues in PDBSite: {list(zip(entry.pdbsite_chains, entry.pdbsite_positions, entry.pdbsite_residues))}")

        # Create site-specific output directory
        site_dir = os.path.join(args.output, f"site_{i+1}_{entry.pdb_id}_{entry.ligand_name.replace(' ', '_')[:10]}")
        os.makedirs(site_dir, exist_ok=True)

        # Download reference PDB
        ref_pdb_path = os.path.join(site_dir, f"{entry.pdb_id.lower()}.pdb")
        print(f"  Downloading reference PDB: {entry.pdb_id}")

        if not download_pdb(entry.pdb_id, ref_pdb_path):
            print(f"  Failed to download PDB {entry.pdb_id}, skipping...")
            continue

        # Check if CIF was downloaded instead
        ref_cif_path = ref_pdb_path.replace('.pdb', '.cif')
        if os.path.exists(ref_cif_path) and not os.path.exists(ref_pdb_path):
            ref_pdb_path = ref_cif_path

        # Extract ligand
        ligand_pdb_raw = os.path.join(site_dir, f"ligand_{entry.ligand_name.replace(' ', '_')[:10]}_raw.pdb")
        print(f"  Extracting ligand from reference PDB")

        ligand_code = extract_ligand_from_pdb(ref_pdb_path, entry.site_descr, ligand_pdb_raw)

        if not ligand_code:
            print(f"  Failed to extract ligand, skipping...")
            continue

        # Superimpose binding sites and transform ligand
        print(f"  Superimposing binding sites...")

        # Get CA coordinates from reference PDB (PDBSite residues)
        ref_residues = list(zip(entry.pdbsite_chains, entry.pdbsite_positions))
        ref_coords, ref_found = get_ca_coords(ref_pdb_path, ref_residues)

        # Get CA coordinates from input protein (protein residues)
        target_residues = list(zip(entry.protein_chains, entry.protein_positions))
        target_coords, target_found = get_ca_coords(cleaned_pdb, target_residues)

        if len(ref_coords) < 3 or len(target_coords) < 3:
            print(f"  Not enough matching residues for superposition (ref: {len(ref_coords)}, target: {len(target_coords)})")
            print(f"  Skipping...")
            continue

        # Match residues by order (they should correspond)
        min_len = min(len(ref_coords), len(target_coords))
        ref_coords = ref_coords[:min_len]
        target_coords = target_coords[:min_len]

        # Compute transformation (ref -> target)
        R, t = kabsch_rotation(ref_coords, target_coords)

        # Transform ligand
        ligand_pdb = os.path.join(site_dir, f"ligand_{entry.ligand_name.replace(' ', '_')[:10]}.pdb")
        transform_ligand_pdb(ligand_pdb_raw, ligand_pdb, R, t)
        print(f"  Ligand transformed to target coordinate system")

        # Combine protein and ligand
        combined_pdb = os.path.join(site_dir, "protein_ligand_combined.pdb")
        print(f"  Combining protein and ligand")
        combine_protein_and_ligand(cleaned_pdb, ligand_pdb, combined_pdb)

        # Format target residues for display
        target_res_list = []
        for chain, pos in zip(entry.protein_chains, entry.protein_positions):
            target_res_list.append(f"{chain}-{pos}")

        # Format ligand specification for PLACER (chain-name3-resno)
        # Ligand is on chain 'L' with resnum 1
        ligand_spec = f"L-{ligand_code}-1"

        # Run PLACER
        suffix = f"{entry.pdb_id}_{ligand_code}"
        print(f"  Running PLACER with predict_ligand: {ligand_spec}")

        csv_path = run_placer(combined_pdb, ligand_spec, args.nsamples, site_dir, suffix, args.placer)

        # Collect results
        result = {
            'site_id': entry.id,
            'ref_pdb': entry.pdb_id,
            'ligand': entry.ligand_name,
            'ligand_code': ligand_code,
            'site_residues': ', '.join(target_res_list),
            'ligand_pdb': ligand_pdb,
            'placer_csv': csv_path
        }

        if csv_path:
            metrics = parse_placer_csv(csv_path)
            result.update(metrics)

        results.append(result)

    # Output summary table
    print(f"\n{'='*60}")
    print("SUMMARY OF RESULTS")
    print('='*60)

    # Write summary CSV
    summary_csv = os.path.join(args.output, 'summary_results.csv')

    if results:
        # Get all keys
        all_keys = set()
        for r in results:
            all_keys.update(r.keys())
        all_keys = sorted(all_keys)

        with open(summary_csv, 'w') as f:
            f.write(','.join(all_keys) + '\n')
            for r in results:
                row = [str(r.get(k, '')) for k in all_keys]
                f.write(','.join(row) + '\n')

        print(f"Summary written to: {summary_csv}")

        # Print table to console
        print("\nConfidence metrics:")
        for r in results:
            print(f"\n  Site: {r.get('site_id', 'N/A')}")
            print(f"    Ligand: {r.get('ligand', 'N/A')} ({r.get('ligand_code', 'N/A')})")
            print(f"    Target residues: {r.get('site_residues', 'N/A')}")
            if 'prmsd' in r:
                print(f"    prmsd: {r.get('prmsd', 'N/A')}")
            if 'plddt' in r:
                print(f"    plddt: {r.get('plddt', 'N/A')}")
            if 'plddt_pde' in r:
                print(f"    plddt_pde: {r.get('plddt_pde', 'N/A')}")
            if 'rmsd' in r:
                print(f"    rmsd: {r.get('rmsd', 'N/A')}")
            print(f"    Ligand PDB: {r.get('ligand_pdb', 'N/A')}")
    else:
        print("No results to report.")

    # Print metrics explanation and save to file
    print_metrics_explanation()
    save_metrics_explanation(args.output)


METRICS_EXPLANATION = """
================================================================================
PLACER METRICS EXPLANATION
================================================================================

CONFIDENCE METRICS (use these to assess prediction quality):
----------------------------------------------------------------------------
  prmsd      - Predicted RMSD: RMS of predicted deviations of atomic positions.
               Lower is better. Primary metric for ranking docking results.

  plddt      - Predicted lDDT score (1D track): Confidence score averaged over
               ligand atoms. Higher is better (0-1 scale).

  plddt_pde  - Predicted lDDT score (2D track): Alternative confidence from
               pairwise distance error prediction. Higher is better (0-1 scale).

ACCURACY METRICS (compare to ground truth, if available):
----------------------------------------------------------------------------
  rmsd       - RMSD of ligand atoms between input and predicted structure.
               Measures docking position accuracy. Lower is better.

  kabsch     - Superimposed RMSD after optimal alignment.
               Measures ligand conformation accuracy. Lower is better.

  lddt       - Actual lDDT score between generated and reference structure.
               Higher is better (0-1 scale).

  fape       - Frame Aligned Point Error: All-atom structural loss.
               Lower is better.

================================================================================
INTERPRETATION GUIDELINES (CUTOFFS)
================================================================================

  CONFIDENCE LEVEL    |  prmsd    |  plddt / plddt_pde
  --------------------|-----------|--------------------
  HIGH (trustworthy)  |  < 2.0    |  > 0.8
  MEDIUM (acceptable) |  2.0-4.0  |  0.6-0.8
  LOW (unreliable)    |  > 4.0    |  < 0.6

RECOMMENDATIONS:
  - Use prmsd as the primary metric for ranking docking predictions
  - High confidence: prmsd < 2.0 with plddt > 0.8
  - Acceptable: prmsd < 4.0 with plddt > 0.8 (especially for complex ligands)
  - For best results, generate 50-100 samples and analyze top 10% by prmsd
  - Large/complex ligands may have higher prmsd but can still be valid if
    plddt scores are good

NOTE: These metrics reflect prediction confidence, not biological activity.
      Always validate results with experimental data when possible.

Reference: https://www.biorxiv.org/content/10.1101/2024.09.25.614868v1
================================================================================
"""


def print_metrics_explanation():
    """Print explanation of PLACER metrics and interpretation guidelines."""
    print(METRICS_EXPLANATION)


def save_metrics_explanation(output_dir: str):
    """Save metrics explanation to a file in the output directory."""
    filepath = os.path.join(output_dir, "METRICS_EXPLANATION.txt")
    with open(filepath, 'w') as f:
        f.write(METRICS_EXPLANATION)
    print(f"Metrics explanation saved to: {filepath}")


if __name__ == '__main__':
    main()
