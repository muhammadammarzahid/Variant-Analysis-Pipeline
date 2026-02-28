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
    
    # For publication imagery, we want distinct tiers
    af_col = 'gnomad_af_genome' if 'gnomad_af_genome' in variants_df.columns else 'maf'
    
    # Check for SIFT 'deleterious' matches or PolyPhen 'damaging' matches
    sift_mask = variants_df.get('sift', pd.Series(dtype=str)).astype(str).str.contains('deleterious', na=False)
    polyphen_mask = variants_df.get('polyphen', pd.Series(dtype=str)).astype(str).str.contains('damaging', na=False)
    
    any_deleterious_mask = sift_mask | polyphen_mask
    dual_deleterious_mask = sift_mask & polyphen_mask
    
    # MAF logic
    maf_series = pd.to_numeric(variants_df.get(af_col, 0), errors='coerce').fillna(0)
    common_mask = maf_series > 0.01
    rare_mask = maf_series <= 0.01
    
    # Missense consequence mask
    conseq_mask = variants_df.get('consequence_type', pd.Series(dtype=str)).astype(str).str.contains('missense', na=False)
    
    # Pre-compute Amino Acid Change column so it's exported in the CSV
    if 'amino_acids' in variants_df.columns and 'protein_start' in variants_df.columns:
        def format_aa_change(row):
            aas = str(row.get('amino_acids', ''))
            pos = str(row.get('protein_start', '')).split('.')[0]
            if pd.isna(row.get('protein_start')) or aas == 'nan' or not aas:
                return ''
            if '/' in aas:
                ref, alt = aas.split('/', 1)
                return f"{ref}{pos}{alt}"
            return f"{aas}{pos}"
        
        # Apply only to rows with protein_start
        mask_ps = variants_df['protein_start'].notna()
        variants_df.loc[mask_ps, 'amino_acid_change'] = variants_df.loc[mask_ps].apply(format_aa_change, axis=1)
    else:
        variants_df['amino_acid_change'] = ''
        
    # Tier 1: High-Confidence Clinical/Multi-Omic
    gwas_mask = variants_df.get('has_gwas', pd.Series(dtype=bool)).astype(bool)
    clinvar_mask = variants_df.get('clinical_significance', pd.Series(dtype=str)).astype(str).str.contains('pathogenic', case=False, na=False)
    tier1_mask = gwas_mask | clinvar_mask
    
    tier1_variants = variants_df[tier1_mask]
    
    # Tier 2: Dual-Evidence Deleterious & Rare
    tier2_mask = conseq_mask & dual_deleterious_mask & rare_mask
    tier2_variants = variants_df[tier2_mask & ~tier1_mask]
    
    # Tier 3: Dual-Evidence Deleterious & Common
    tier3_mask = conseq_mask & dual_deleterious_mask & common_mask
    tier3_variants = variants_df[tier3_mask & ~tier1_mask & ~tier2_mask]
    
    # Tier 4: eQTL Missense Variants
    # Check that they have an eQTL association AND the minimum p-value is significant enough
    eqtl_mask = variants_df.get('has_eqtl', pd.Series(dtype=bool)).astype(bool)
    eqtl_pvals = pd.to_numeric(variants_df.get('eqtl_min_pvalue', 1.0), errors='coerce').fillna(1.0)
    significant_eqtl_mask = eqtl_mask & (eqtl_pvals < 1e-5)
    
    tier4_mask = conseq_mask & significant_eqtl_mask
    tier4_variants = variants_df[tier4_mask & ~tier1_mask & ~tier2_mask & ~tier3_mask]
    
    # Tier 5: Common Polymorphisms (Structurally Tolerated)
    tier5_mask = conseq_mask & common_mask & ~any_deleterious_mask
    tier5_variants = variants_df[tier5_mask & ~tier1_mask & ~tier2_mask & ~tier3_mask & ~tier4_mask]
    
    # Store variants we want to export
    export_dfs = []

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
    plotted_pos = set()
    
    # Plot Tier 1: High Confidence / Proven Impact (Red)
    if not tier1_variants.empty and 'protein_start' in tier1_variants.columns:
        for _, row in tier1_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            aa = row.get('amino_acid_change', '')
            
            name = f"t1_mut_{pos}"
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color red, {name}")
            pml.append(f"set sphere_scale, 1.5, {name}")
            
            var_id = row.get('variant_id', f'pos_{pos}')
            pml.append(f"label {name} and name CA, '{var_id} {aa}'")
            count += 1
            plotted_pos.add(pos)
            
        t1_export = tier1_variants.copy()
        t1_export['visualization_tier'] = 'Tier 1: Clinical/GWAS'
        export_dfs.append(t1_export)
            
    # Plot Tier 2: Dual-Evidence Deleterious & Rare (Orange)
    if not tier2_variants.empty and 'protein_start' in tier2_variants.columns:
        for _, row in tier2_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in plotted_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"t2_del_rare_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color orange, {name}")
            pml.append(f"set sphere_scale, 1.3, {name}")
            
            var_id = row.get('variant_id', f'pos_{pos}')
            pml.append(f"label {name} and name CA, '{var_id} {aa}'")
            count += 1
            plotted_pos.add(pos)
            
        t2_export = tier2_variants[~tier2_variants['protein_start'].astype(int).isin(plotted_pos - set(tier2_variants['protein_start'].dropna().astype(int)))].copy()
        if not t2_export.empty:
            t2_export['visualization_tier'] = 'Tier 2: Dual-Evidence Rare'
            export_dfs.append(t2_export)

    # Plot Tier 3: Dual-Evidence Deleterious & Common (Magenta)
    if not tier3_variants.empty and 'protein_start' in tier3_variants.columns:
        for _, row in tier3_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in plotted_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"t3_del_com_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color magenta, {name}")
            pml.append(f"set sphere_scale, 1.2, {name}")
            
            var_id = row.get('variant_id', f'pos_{pos}')
            pml.append(f"label {name} and name CA, '{var_id} {aa}'")
            count += 1
            plotted_pos.add(pos)
            
        t3_export = tier3_variants[~tier3_variants['protein_start'].astype(int).isin(plotted_pos - set(tier3_variants['protein_start'].dropna().astype(int)))].copy()
        if not t3_export.empty:
            t3_export['visualization_tier'] = 'Tier 3: Dual-Evidence Common'
            export_dfs.append(t3_export)

    # Plot Tier 4: eQTL Missense Variants (Yellow)
    if not tier4_variants.empty and 'protein_start' in tier4_variants.columns:
        for _, row in tier4_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in plotted_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"t4_eqtl_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color yellow, {name}")
            pml.append(f"set sphere_scale, 1.1, {name}")
            count += 1
            plotted_pos.add(pos)
            
        t4_export = tier4_variants[~tier4_variants['protein_start'].astype(int).isin(plotted_pos - set(tier4_variants['protein_start'].dropna().astype(int)))].copy()
        if not t4_export.empty:
            t4_export['visualization_tier'] = 'Tier 4: eQTL Missense'
            export_dfs.append(t4_export)
            
    # Plot Tier 5: Common Polymorphisms Tolerated (Cyan)
    if not tier5_variants.empty and 'protein_start' in tier5_variants.columns:
        for _, row in tier5_variants.dropna(subset=['protein_start']).iterrows():
            pos = int(float(row['protein_start']))
            if pos in plotted_pos:
                continue
                
            aa = row.get('amino_acid_change', '')
            name = f"t5_tol_{pos}"
            
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color cyan, {name}")
            pml.append(f"set sphere_scale, 0.8, {name}")
            count += 1
            plotted_pos.add(pos)
            
        t5_export = tier5_variants[~tier5_variants['protein_start'].astype(int).isin(plotted_pos - set(tier5_variants['protein_start'].dropna().astype(int)))].copy()
        if not t5_export.empty:
            t5_export['visualization_tier'] = 'Tier 5: Common Tolerated'
            export_dfs.append(t5_export)
            
    # Fallback if we have very few variants mapped
    if count < 5 and 'protein_start' in variants_df.columns:
        general_missense = variants_df[variants_df['consequence_type'].astype(str).str.contains('missense', na=False)]
        
        added = 0
        fallback_added = []
        for _, row in general_missense.dropna(subset=['protein_start']).iterrows():
            if added >= 10: break
            
            pos = int(float(row['protein_start']))
            if pos in plotted_pos: continue
            
            aa = row.get('amino_acid_change', '')
            name = f"t6_missense_{pos}"
            pml.append(f"select {name}, resi {pos}")
            pml.append(f"show spheres, {name}")
            pml.append(f"color gray50, {name}")
            pml.append(f"set sphere_scale, 0.6, {name}")
            plotted_pos.add(pos)
            added += 1
            count += 1
            fallback_added.append(row)
            
        if fallback_added:
            t6_export = pd.DataFrame(fallback_added)
            t6_export['visualization_tier'] = 'Tier 6: Random Fallback'
            export_dfs.append(t6_export)
            
    pml.append("set label_color, black")
    pml.append("set label_size, 24")
    pml.append("set label_font_id, 7") # bold
    pml.append(f"zoom {gene_pymol_name}")
    pml.append("ray 1200, 1000")
    
    with open(output_script, 'w') as f:
        f.write("\n".join(pml))
        
    print(f"  ✓ PyMOL script generated with {count} mutations to visualize: {output_script}")
    
    # Export the table of prioritized variants
    if export_dfs:
        final_export_df = pd.concat(export_dfs, ignore_index=True)
        # Reorder columns to put tier first
        cols = ['visualization_tier'] + [c for c in final_export_df.columns if c != 'visualization_tier']
        final_export_df = final_export_df[cols]
        export_path = Path('outputs/prioritized_structural_variants.csv')
        export_path.parent.mkdir(parents=True, exist_ok=True)
        final_export_df.to_csv(export_path, index=False)
        print(f"  ✓ Prioritized variants table saved to: {export_path}")
    else:
        print("  ! No variants met criteria for PyMOL structural visualization.")

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
