#!/usr/bin/env python3
"""
Script 11: Visualize Priority Mutations on Structure
Downloads the AlphaFold PDB for Gene and generates a PyMOL script to map priority missense mutants.
"""

import sys
from pathlib import Path
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from utils.api_clients import AlphaFoldClient
from utils.constants import UNIPROT_ID, GENE_SYMBOL

def generate_pymol_script(variants_df, pdb_path, output_script):
    """Generates a .pml script to visualize variants in PyMOL."""
    
    # Filter to only the high priority mutations (e.g. missense with GWAS hits or SIFT deleterious)
    # OR just visualize all missense mutations if it's small, otherwise top N
    
    # For publication imagery, we want GWAS variants colored distinctly from others
    gwas_variants = variants_df[variants_df.get('has_gwas', False) == True]
    
    # Also grab frequency info, choosing genome over exome or MAF
    af_col = 'gnomad_af_genome' if 'gnomad_af_genome' in variants_df.columns else 'maf'
    
    # Check for SIFT 'deleterious' matches or PolyPhen 'damaging' matches
    sift_mask = variants_df.get('sift', pd.Series(dtype=str)).astype(str).str.contains('deleterious', na=False)
    polyphen_mask = variants_df.get('polyphen', pd.Series(dtype=str)).astype(str).str.contains('damaging', na=False)
    
    # Common variant threshold (e.g. > 1% MAF)
    common_mask = pd.to_numeric(variants_df.get(af_col, 0), errors='coerce').fillna(0) > 0.01
    
    # Tier 1: GWAS
    # Tier 2: Common Missense
    conseq_mask = variants_df.get('consequence_type', pd.Series(dtype=str)).astype(str).str.contains('missense', na=False)
    common_missense = variants_df[common_mask & conseq_mask]
    
    # Tier 3: Rare Deleterious (filtered so very rare ones are deprioritized compared to common variants)
    # Let's say we only get deleterious variants that are at least > 0.0001 (0.01%) so we ignore ultra-rare noise
    rare_but_real_mask = pd.to_numeric(variants_df.get(af_col, 0), errors='coerce').fillna(0) >= 0.0001
    deleterious_filtered = variants_df[(sift_mask | polyphen_mask) & rare_but_real_mask]
    
    pml = []
    
    # Use absolute path for robust loading regardless of PyMOL's working directory
    # Wrap in quotes to handle any spaces in the directory path
    abs_pdb_path = Path(pdb_path).resolve()
    gene_pymol_name = GENE_SYMBOL.lower()
    
    pml.append(f"load \"{abs_pdb_path}\", {gene_pymol_name}")
    pml.append("hide everything, name H*")
    pml.append("hide everything, hydrogens")
    pml.append("bg_color white")
    pml.append(f"color gray80, {gene_pymol_name}")
    
    # Beautiful cartoon representation
    pml.append(f"show cartoon, {gene_pymol_name}")
    pml.append("set cartoon_transparency, 0.2")
    
    count = 0
    gwas_pos = set()
    
    # Plot GWAS 
    if not gwas_variants.empty and 'protein_start' in gwas_variants.columns:
        gwas_pos = set(gwas_variants['protein_start'].dropna().astype(int))
        for _, row in gwas_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            aa = row.get('amino_acid_change', '')
            name = f"gwas_mut_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color gold, {name}")
            pml.append(f"set sphere_scale, 1.5, {name}")
            
            var_id = row.get('variant_id', f'pos_{pos}')
            pml.append(f"label {name} and name CA, '{var_id} {aa}'")
            count += 1
            
    # Plot Common Missense (Cyan)
    if not common_missense.empty and 'protein_start' in common_missense.columns:
        for _, row in common_missense.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in gwas_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"common_mut_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color cyan, {name}")
            count += 1
            
    # Plot SIFT Deleterious (filtered for ultra-rare noise)
    if not deleterious_filtered.empty and 'protein_start' in deleterious_filtered.columns:
        # Avoid double-plotting if GWAS/Common is also deleterious
        plotted_pos = set(gwas_variants['protein_start'].dropna().astype(int)) if not gwas_variants.empty else set()
        if not common_missense.empty:
            plotted_pos.update(common_missense['protein_start'].dropna().astype(int))
            
        for _, row in deleterious_filtered.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in plotted_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"del_mut_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color red, {name}")
            
            # Optional: label top few deleterious mutations
            # pml.append(f"label {name} and name CA, '{aa}'")
            count += 1
            
    # If we have very few, plot some general missense
    if count < 5 and 'protein_start' in variants_df.columns:
        general_missense = variants_df[variants_df['consequence_type'].str.contains('missense', na=False)]
        plotted_pos = gwas_pos.copy() if not gwas_variants.empty else set()
        if not common_missense.empty:
            plotted_pos.update(common_missense['protein_start'].dropna().astype(int))
        if not deleterious_filtered.empty:
            plotted_pos.update(deleterious_filtered['protein_start'].dropna().astype(int))
            
        added = 0
        for _, row in general_missense.dropna(subset=['protein_start']).iterrows():
            if added >= 10: break
            
            pos = int(float(row['protein_start']))
            if pos in plotted_pos: continue
            
            aa = row.get('amino_acid_change', '')
            name = f"missense_{pos}"
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color deepblue, {name}")
            plotted_pos.add(pos)
            added += 1
            count += 1
            
    pml.append("set label_color, black")
    pml.append("set label_size, 24")
    pml.append("set label_font_id, 7") # bold
    pml.append(f"zoom {gene_pymol_name}")
    pml.append("ray 1200, 1000")
    
    with open(output_script, 'w') as f:
        f.write("\n".join(pml))
        
    print(f"  ✓ PyMOL script generated with {count} mutations to visualize: {output_script}")

def main():
    print("=" * 80)
    print("Structure Fetcher & PyMOL Visualizer")
    print("=" * 80)
    
    structures_dir = Path('structures')
    structures_dir.mkdir(parents=True, exist_ok=True)
    
    pdb_path = structures_dir / f'{GENE_SYMBOL.lower()}_alphafold.pdb'
    
    # 1. Download PDB
    print("\n[1/2] Fetching Gene AlphaFold structure...")
    if not pdb_path.exists():
        af_client = AlphaFoldClient()
        success = af_client.download_structure(UNIPROT_ID, str(pdb_path))
        if success:
            print(f"  ✓ Structure downloaded to {pdb_path}")
        else:
            print(f"  ✗ Failed to download structure for {UNIPROT_ID}")
            return
    else:
        print(f"  ✓ Structure already exists at {pdb_path}")
        
    # 2. Map Priority Mutations
    print("\n[2/2] Mapping priority mutations...")
    master_file = Path('outputs/master_integrated.csv')
    
    if not master_file.exists():
        print("  ✗ Master data file not found. Run previous scripts first.")
        return
        
    master_df = pd.read_csv(master_file, low_memory=False)
    output_script = structures_dir / 'visualize_mutations.pml'
    
    generate_pymol_script(master_df, pdb_path, output_script)
    
    print("\n" + "=" * 80)
    print("VISUALIZATION READY")
    print("=" * 80)
    print("To view interacting variants in PyMOL run:")
    print(f"  pymol {output_script}")
    print("=" * 80)
    print("\n✓ Structure mapping complete!\n")

if __name__ == '__main__':
    main()
