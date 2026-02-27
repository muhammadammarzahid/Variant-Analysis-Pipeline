#!/usr/bin/env python3
"""
Script 09: Integrate All Data
Combines all data sources into a master dataset and generates summary report.
"""

import json
import sys
from pathlib import Path
import pandas as pd
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))

from utils.constants import DATA_PATHS, GENE_SYMBOL, GENE_ID


def main():
    print("=" * 80)
    print("Data Integration and Summary Generator")
    print("=" * 80)
    
    # Load all data
    print(f"\n[1/5] Loading all datasets...")
    
    # Variants
    variants_file = Path(DATA_PATHS['processed']) / 'variants_with_frequencies.csv'
    if not variants_file.exists():
        variants_file = Path(DATA_PATHS['processed']) / 'variants_basic.csv'
        
    if variants_file.exists():
        variants_df = pd.read_csv(variants_file)
        print(f"  ✓ Variants: {len(variants_df):,}")
    else:
        print("  ✗ Variants file not found")
        return
    
    # eQTL
    eqtl_file = Path('outputs') / 'eqtl_associations.csv'
    if eqtl_file.exists():
        eqtl_df = pd.read_csv(eqtl_file)
        print(f"  ✓ eQTL: {len(eqtl_df):,} associations")
    else:
        eqtl_df = pd.DataFrame()
        print("  ⚠ eQTL file not found")
    
    # GWAS
    gwas_file = Path('outputs') / 'gwas_associations.csv'
    if gwas_file.exists():
        gwas_df = pd.read_csv(gwas_file)
        print(f"  ✓ GWAS: {len(gwas_df):,} associations")
    else:
        gwas_df = pd.DataFrame()
        print("  ⚠ GWAS file not found")
    
    # Annotations
    annot_file = Path('outputs') / 'variants_annotated.csv'
    if annot_file.exists():
        annot_df = pd.read_csv(annot_file)
        print(f"  ✓ Annotations: {len(annot_df):,}")
    else:
        annot_df = pd.DataFrame()
        print("  ⚠ Annotations file not found")
    
    # Domain mapping
    domain_file = Path('outputs') / 'variant_protein_mapping.csv'
    if domain_file.exists():
        domain_df = pd.read_csv(domain_file)
        print(f"  ✓ Domain mapping: {len(domain_df):,}")
    else:
        domain_df = pd.DataFrame()
        print("  ⚠ Domain mapping file not found")
    
    # Merge datasets
    print(f"\n[2/5] Integrating datasets...")
    
    master_df = variants_df.copy()
    
    # Add eQTL info
    if not eqtl_df.empty and 'rsid' in eqtl_df.columns:
        # Summarize by rsID
        eqtl_summary = eqtl_df.dropna(subset=['rsid']).groupby('rsid').agg({
            'tissue': lambda x: '|'.join(x.unique()),
            'pvalue': 'min',
            'variant_id': 'first'  # Keep GTEx variant ID
        }).reset_index()
        eqtl_summary.columns = ['rsid', 'eqtl_tissues', 'eqtl_min_pvalue', 'eqtl_variant_id']
        eqtl_summary['has_eqtl'] = True
        
        # Merge by standard rsID
        master_df = master_df.merge(eqtl_summary, left_on='variant_id', right_on='rsid', how='left')
        master_df['has_eqtl'] = master_df['has_eqtl'].fillna(False)
        master_df = master_df.drop('rsid', axis=1)  # Remove temporary rsid column
    else:
        master_df['has_eqtl'] = False
    
    # Add GWAS info
    if not gwas_df.empty:
        gwas_summary = gwas_df.groupby('variant_id').agg({
            'trait': lambda x: '|'.join(x.unique()) if len(x) > 0 else ''
        }).reset_index()
        gwas_summary.columns = ['variant_id', 'gwas_traits']
        gwas_summary['has_gwas'] = True
        master_df = master_df.merge(gwas_summary, on='variant_id', how='left')
        master_df['has_gwas'] = master_df['has_gwas'].fillna(False)
    else:
        master_df['has_gwas'] = False
        
    # Add annotations
    if not annot_df.empty and 'variant_id' in annot_df.columns:
        annot_sub = annot_df.dropna(subset=['variant_id']).drop_duplicates(subset=['variant_id'])
        cols_to_use = ['variant_id', 'sift', 'polyphen', 'amino_acids', 'codons', 'protein_start', 'protein_end']
        cols_to_use = [c for c in cols_to_use if c in annot_sub.columns]
        master_df = master_df.merge(annot_sub[cols_to_use], on='variant_id', how='left')
        
    # Add domain mapping
    if not domain_df.empty and 'variant_id' in domain_df.columns:
        domain_sub = domain_df.dropna(subset=['variant_id']).drop_duplicates(subset=['variant_id'])
        cols_to_use = ['variant_id', 'domain_affected', 'domain_start', 'domain_end']
        cols_to_use = [c for c in cols_to_use if c in domain_sub.columns]
        master_df = master_df.merge(domain_sub[cols_to_use], on='variant_id', how='left')
    
    print(f"  ✓ Integrated {len(master_df):,} variants")
    
    # Calculate summary statistics
    print(f"\n[3/5] Calculating summary statistics...")
    
    consequence_counts = {}
    for cons in master_df['consequence_type'].dropna():
        if cons and cons.strip():
            for c in cons.split('|'):
                c = c.strip()
                consequence_counts[c] = consequence_counts.get(c, 0) + 1
    
    top_consequences = dict(sorted(consequence_counts.items(), key=lambda x: x[1], reverse=True)[:10])
    
    freq_col = None
    if 'gnomad_af_genome' in master_df.columns:
        freq_col = 'gnomad_af_genome'
    elif 'gnomad_af' in master_df.columns:
        freq_col = 'gnomad_af'
    
    freq_breakdown = {
        'very_rare (<0.1%)': 0,
        'rare (0.1%-1%)': 0,
        'low_frequency (1%-5%)': 0,
        'common (>5%)': 0,
        'unknown/missing': len(master_df)
    }
    
    if freq_col:
        af_values = pd.to_numeric(master_df[freq_col], errors='coerce').fillna(-1)
        freq_breakdown['very_rare (<0.1%)'] = int(((af_values >= 0) & (af_values < 0.001)).sum())
        freq_breakdown['rare (0.1%-1%)'] = int(((af_values >= 0.001) & (af_values < 0.01)).sum())
        freq_breakdown['low_frequency (1%-5%)'] = int(((af_values >= 0.01) & (af_values < 0.05)).sum())
        freq_breakdown['common (>5%)'] = int((af_values >= 0.05).sum())
        freq_breakdown['unknown/missing'] = int((af_values == -1).sum())

    stats = {
        'gene_symbol': GENE_SYMBOL,
        'gene_id': GENE_ID,
        'analysis_date': datetime.now().isoformat(),
        'total_variants': int(len(master_df)),
        'coding_variants': int(master_df['consequence_type'].str.contains('missense|nonsense|frameshift', na=False).sum()),
        'with_eqtl': int(master_df['has_eqtl'].sum()) if 'has_eqtl' in master_df else 0,
        'with_gwas': int(master_df['has_gwas'].sum()) if 'has_gwas' in master_df else 0,
        'clinical_significance': int(master_df[master_df['clinical_significance'] != ''].shape[0]) if 'clinical_significance' in master_df.columns else 0,
        'top_consequence_types': top_consequences,
        'frequency_breakdown': freq_breakdown,
    }
    
    # Save integrated data
    print(f"\n[4/5] Saving integrated dataset...")
    
    output_file = Path('outputs') / 'master_integrated.csv'
    master_df.to_csv(output_file, index=False)
    print(f"  ✓ Master dataset saved to {output_file}")
    
    stats_file = Path('outputs') / 'summary_stats.json'
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"  ✓ Statistics saved to {stats_file}")
    
    # Generate summary report
    print(f"\n[5/5] Generating summary report...")
    
    report_file = Path('outputs') / 'summary_report.txt'
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"Gene MUTATION AND QTL ANALYSIS SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Gene: {GENE_SYMBOL} ({GENE_ID})\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("OVERVIEW\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total Variants Found:              {stats['total_variants']:,}\n")
        f.write(f"Coding Variants:                   {stats['coding_variants']:,}\n")
        f.write(f"Non-coding Variants:               {stats['total_variants'] - stats['coding_variants']:,}\n")
        f.write(f"With Clinical Significance:        {stats['clinical_significance']:,}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("TOP 10 VARIANT CONSEQUENCE TYPES\n")
        f.write("=" * 80 + "\n")
        for cons, count in top_consequences.items():
            f.write(f"{cons:45s} {count:6,}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("VARIANT FREQUENCIES (gnomAD)\n")
        f.write("=" * 80 + "\n")
        f.write(f"Common (>5%):                      {freq_breakdown['common (>5%)']:,}\n")
        f.write(f"Low Frequency (1%-5%):             {freq_breakdown['low_frequency (1%-5%)']:,}\n")
        f.write(f"Rare (0.1%-1%):                    {freq_breakdown['rare (0.1%-1%)']:,}\n")
        f.write(f"Very Rare (<0.1%):                 {freq_breakdown['very_rare (<0.1%)']:,}\n")
        f.write(f"Unknown/Missing:                   {freq_breakdown['unknown/missing']:,}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("REGULATORY EFFECTS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Variants with eQTL:                {stats['with_eqtl']:,}\n")
        if not eqtl_df.empty:
            f.write(f"Total eQTL Associations:           {len(eqtl_df):,}\n")
            f.write(f"Tissues with eQTLs:                {eqtl_df['tissue'].nunique()}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("GWAS ASSOCIATIONS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Variants with GWAS hits:           {stats['with_gwas']:,}\n")
        if not gwas_df.empty:
            f.write(f"Total GWAS Associations:           {len(gwas_df):,}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("KEY FINDINGS\n")
        f.write("=" * 80 + "\n")
        
        # Highlight interesting variants
        multi_evidence = master_df[
            (master_df.get('has_eqtl', False)) & 
            (master_df.get('has_gwas', False))
        ]
        
        f.write(f"Variants with Multiple Evidence:\n")
        f.write(f"  - eQTL + GWAS:                   {len(multi_evidence):,}\n\n")
        
        # Coding variants
        missense = master_df[master_df['consequence_type'].str.contains('missense', na=False)]
        stop_gained = master_df[master_df['consequence_type'].str.contains('stop_gained', na=False)]
        frameshift = master_df[master_df['consequence_type'].str.contains('frameshift', na=False)]
        
        f.write(f"High-Impact Coding Variants:\n")
        f.write(f"  - Missense:                      {len(missense):,}\n")
        f.write(f"  - Stop-gained (nonsense):        {len(stop_gained):,}\n")
        f.write(f"  - Frameshift:                    {len(frameshift):,}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("DATA FILES GENERATED\n")
        f.write("=" * 80 + "\n")
        f.write("- master_integrated.csv:     All variants with integrated annotations\n")
        f.write("- variants_basic.csv:        Basic variant information\n")
        f.write("- eqtl_associations.csv:     Expression QTL data\n")
        f.write("- gwas_associations.csv:     GWAS phenotype associations\n")
        f.write("- domain_mapping.csv:        Coding variants mapped to domains\n")
        f.write("- variants_annotated.csv:    VEP functional annotations\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("NEXT STEPS\n")
        f.write("=" * 80 + "\n")
        f.write("1. Review high-impact coding variants (missense, nonsense, frameshift)\n")
        f.write("2. Investigate variants with both eQTL and GWAS evidence\n")
        f.write("3. Examine tissue-specific eQTL patterns\n")
        f.write("4. Map coding variants to protein domains for structural analysis\n")
        f.write("5. Validate top candidates experimentally\n\n")
        
        f.write("=" * 80 + "\n")
        f.write(f"Report generated by Gene Analysis Pipeline\n")
        f.write("=" * 80 + "\n")
    
    print(f"  ✓ Report saved to {report_file}")
    
    # Print summary to console
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total Variants:           {stats['total_variants']:,}")
    print(f"Coding Variants:          {stats['coding_variants']:,}")
    print(f"With eQTL:                {stats['with_eqtl']:,}")
    print(f"With GWAS:                {stats['with_gwas']:,}")
    print(f"Clinical Significance:    {stats['clinical_significance']:,}")
    print("\nOutputs saved to: outputs/")
    print(f"  - master_integrated.csv")
    print(f"  - summary_report.txt")
    print(f"  - summary_stats.json")
    print("=" * 80)
    print("\n✓ Integration complete!\n")


if __name__ == '__main__':
    main()
