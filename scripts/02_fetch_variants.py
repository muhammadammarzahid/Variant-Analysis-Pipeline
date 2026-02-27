#!/usr/bin/env python3
"""
Script 02: Fetch Gene Variants from Ensembl
Retrieves all known variants in the Gene gene region.
"""

import json
import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import EnsemblClient
from utils.constants import (
    GENE_ID, GENE_SYMBOL, GENE_CHR, 
    GENE_START, GENE_END, DATA_PATHS
)


def main():
    print("=" * 80)
    print("Gene Variant Fetcher")
    print("=" * 80)
    
    # Initialize client
    ensembl = EnsemblClient()
    
    # Fetch variants in Gene region
    print(f"\n[1/2] Fetching variants in region chr{GENE_CHR}:{GENE_START}-{GENE_END}...")
    print(f"  (This may take a few minutes due to API rate limiting)")
    
    variants = ensembl.get_variants_in_region(GENE_CHR, GENE_START, GENE_END)
    
    if not variants:
        print("  ✗ Failed to fetch variants")
        return
    
    print(f"  ✓ Found {len(variants)} variants")
    
    # Process variants
    print(f"\n[2/2] Processing variant information...")
    
    variant_data = []
    
    for v in tqdm(variants, desc="Processing variants"):
        # Extract basic info
        variant_id = v.get('id', '')
        
        # Skip if not a proper variant ID
        if not variant_id or not variant_id.startswith('rs'):
            continue
        
        variant_info = {
            'variant_id': variant_id,
            'chr': v.get('seq_region_name', GENE_CHR),
            'start': v.get('start'),
            'end': v.get('end'),
            'strand': v.get('strand'),
            'alleles': '|'.join(v.get('alleles', [])) if isinstance(v.get('alleles'), list) else str(v.get('alleles', '')),
            'minor_allele': v.get('minor_allele', ''),
            'maf': v.get('minor_allele_freq', None),
            'consequence_type': '|'.join(v.get('consequence_type', [])) if isinstance(v.get('consequence_type'), list) else str(v.get('consequence_type', '')),
            'clinical_significance': '|'.join(v.get('clinical_significance', [])) if isinstance(v.get('clinical_significance'), list) else str(v.get('clinical_significance', '')),
            'source': v.get('source', ''),
        }
        
        variant_data.append(variant_info)
    
    print(f"  ✓ Processed {len(variant_data)} valid variants with rsIDs")
    
    # Save raw JSON
    print(f"\n[Saving] Writing results to disk...")
    raw_dir = Path(DATA_PATHS['raw'])
    raw_dir.mkdir(parents=True, exist_ok=True)
    
    raw_file = raw_dir / 'ensembl_variants.json'
    with open(raw_file, 'w') as f:
        json.dump(variants, f, indent=2)
    print(f"  ✓ Raw data saved to {raw_file}")
    
    # Create DataFrame
    df = pd.DataFrame(variant_data)
    
    # Calculate variant statistics
    stats = {
        'total_variants': int(len(df)),
        'variants_with_maf': int(df['maf'].notna().sum()),
        'variants_with_consequence': int(df['consequence_type'].notna().sum()),
        'variants_with_clinical_sig': int(df[df['clinical_significance'] != ''].shape[0]),
        'unique_consequence_types': len(set([c for cons in df['consequence_type'].dropna() 
                                             for c in cons.split(',') if c])),
    }
    
    # Count consequence types
    all_consequences = []
    for cons in df['consequence_type'].dropna():
        if cons and cons.strip():  # Make sure it's not empty
            all_consequences.extend([c.strip() for c in cons.split('|') if c.strip()])
    consequence_counts = pd.Series(all_consequences).value_counts()
    
    # Save processed CSV
    processed_dir = Path(DATA_PATHS['processed'])
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    csv_file = processed_dir / 'variants_basic.csv'
    df.to_csv(csv_file, index=False)
    print(f"  ✓ Processed data saved to {csv_file}")
    
    # Save statistics
    stats_file = processed_dir / 'variant_stats.json'
    with open(stats_file, 'w') as f:
        json.dump({
            **stats,
            'top_consequences': {k: int(v) for k, v in consequence_counts.head(10).items()}
        }, f, indent=2)
    print(f"  ✓ Statistics saved to {stats_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total Variants Found:     {stats['total_variants']:,}")
    print(f"Variants with MAF:        {stats['variants_with_maf']:,}")
    print(f"With Consequence Type:    {stats['variants_with_consequence']:,}")
    print(f"With Clinical Sig:        {stats['variants_with_clinical_sig']:,}")
    print(f"\nTop 10 Consequence Types:")
    for cons, count in consequence_counts.head(10).items():
        print(f"  {cons:40s} {count:6,}")
    print("=" * 80)
    print("\n✓ Variant fetching complete!\n")


if __name__ == '__main__':
    main()
