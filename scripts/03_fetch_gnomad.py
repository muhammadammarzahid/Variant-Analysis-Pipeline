#!/usr/bin/env python3
"""
Script 03: Fetch gnomAD Population Frequencies
Note: gnomAD API access is limited. This script will attempt to query gnomAD
or fall back to using Ensembl variant annotations which include some frequency data.
"""

import json
import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import GnomADClient, EnsemblClient
from utils.constants import GENE_ID, GENE_SYMBOL, DATA_PATHS


def main():
    print("=" * 80)
    print("gnomAD Population Frequency Fetcher")
    print("=" * 80)
    
    # Load variants from previous step
    print(f"\n[1/3] Loading variants from previous step...")
    variants_file = Path(DATA_PATHS['processed']) / 'variants_basic.csv'
    
    if not variants_file.exists():
        print(f"  ✗ File not found: {variants_file}")
        print("  Please run 02_fetch_variants.py first")
        return
    
    df = pd.read_csv(variants_file)
    print(f"  ✓ Loaded {len(df)} variants")
    
    from utils.api_clients import MyVariantClient
    
    # Try MyVariant API for batch processing
    print(f"\n[2/3] Querying MyVariant.info for {len(df)} variants...")
    myvariant = MyVariantClient()
    
    # Extract rsIDs where available (variant_id should be rsID)
    rsids_to_query = df['variant_id'].dropna().tolist()
    rsids_to_query = [vid for vid in rsids_to_query if str(vid).startswith('rs')]
    
    # Also we might have some variants without rsID but we can only query MyVariant by rsID easily.
    # The script currently uses 'variant_id' column for queries.
    
    print(f"  ✓ Found {len(rsids_to_query)} variants with rsIDs to query")
    
    myvariant_results = myvariant.get_variants_bulk(rsids_to_query)
    
    freq_data = []
    
    if myvariant_results:
        print(f"  ✓ Received data for {len(myvariant_results)} variants")
        for res in myvariant_results:
            variant_id = res.get('query')
            if not variant_id:
                continue
                
            genome = res.get('gnomad_genome', {})
            exome = res.get('gnomad_exome', {})
            
            # gnomAD genome
            genome_af = genome.get('af', {})
            genome_ac = genome.get('ac', {})
            genome_an = genome.get('an', {})
            genome_hom = genome.get('hom', {})
            
            # gnomAD exome
            exome_af = exome.get('af', {})
            exome_ac = exome.get('ac', {})
            exome_an = exome.get('an', {})
            exome_hom = exome.get('hom', {})
            
            freq_info = {
                'variant_id': variant_id,
                'gnomad_af_genome': genome_af.get('af') if genome_af else None,
                'gnomad_ac_genome': genome_ac.get('ac') if genome_ac else None,
                'gnomad_an_genome': genome_an.get('an') if genome_an else None,
                'gnomad_hom_genome': genome_hom.get('hom') if genome_hom else None,
                'gnomad_af_exome': exome_af.get('af') if exome_af else None,
                'gnomad_ac_exome': exome_ac.get('ac') if exome_ac else None,
                'gnomad_an_exome': exome_an.get('an') if exome_an else None,
                'gnomad_hom_exome': exome_hom.get('hom') if exome_hom else None,
            }
            freq_data.append(freq_info)
            
        print(f"  ✓ Processed frequencies for {len(freq_data)} variants")
    else:
        print("  ⚠ MyVariant.info API did not return data or failed.")
    
    # Create DataFrame with frequency data
    freq_df = pd.DataFrame(freq_data)
    
    # Merge with original variants
    merged_df = df.merge(freq_df, on='variant_id', how='left')
    
    # Calculate statistics
    print(f"\n[Statistics] Calculating frequency statistics...")
    stats = {
        'total_variants': len(merged_df),
        'with_genome_freq': merged_df['gnomad_af_genome'].notna().sum(),
        'with_exome_freq': merged_df['gnomad_af_exome'].notna().sum(),
        'common_variants_genome': (merged_df['gnomad_af_genome'] > 0.01).sum() if 'gnomad_af_genome' in merged_df else 0,
        'rare_variants_genome': ((merged_df['gnomad_af_genome'] > 0) & (merged_df['gnomad_af_genome'] <= 0.01)).sum() if 'gnomad_af_genome' in merged_df else 0,
    }
    
    # Save results
    print(f"\n[Saving] Writing results to disk...")
    processed_dir = Path(DATA_PATHS['processed'])
    
    output_file = processed_dir / 'variants_with_frequencies.csv'
    merged_df.to_csv(output_file, index=False)
    print(f"  ✓ Data saved to {output_file}")
    
    stats_file = processed_dir / 'frequency_stats.json'
    with open(stats_file, 'w') as f:
        json_stats = {k: int(v) if hasattr(v, 'item') else v for k, v in stats.items()}
        json.dump(json_stats, f, indent=2)
    print(f"  ✓ Statistics saved to {stats_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total Variants:           {stats['total_variants']:,}")
    print(f"With Genome Frequency:    {stats['with_genome_freq']:,}")
    print(f"With Exome Frequency:     {stats['with_exome_freq']:,}")
    print(f"Common (AF > 1%):         {stats['common_variants_genome']:,}")
    print(f"Rare (AF ≤ 1%):           {stats['rare_variants_genome']:,}")
    print("=" * 80)
    print("\n✓ Frequency data collection complete!\n")
    print("Note: Full gnomAD data may require direct API access or bulk downloads.")
    print("Consider using the gnomAD browser for comprehensive frequency data.")


if __name__ == '__main__':
    main()
