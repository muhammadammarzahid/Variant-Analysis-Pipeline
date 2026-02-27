#!/usr/bin/env python3
"""
Script 06: Fetch GWAS Associations
Retrieves GWAS associations for Gene variants.
"""

import json
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import GWASCatalogClient
from utils.constants import GENE_SYMBOL, GENE_CHR, GENE_START, GENE_END, DATA_PATHS


def main():
    print("=" * 80)
    print("GWAS Catalog Fetcher")
    print("=" * 80)
    
    gwas = GWASCatalogClient()
    
    # Method 1: Query by gene name (returns nearby variants)
    print(f"\n[1/3] Querying GWAS Catalog by gene name: {GENE_SYMBOL}...")
    data_gene = gwas.get_associations_by_gene(GENE_SYMBOL)
    
    associations = []
    if data_gene and '_embedded' in data_gene:
        snps = data_gene['_embedded'].get('singleNucleotidePolymorphisms', [])
        print(f"  ✓ Found {len(snps)} SNPs near Gene (may be outside gene boundaries)")
        
        for snp in snps:
            snp_id = snp.get('rsId', '')
            chr_name = snp.get('locations', [{}])[0].get('chromosomeName') if snp.get('locations') else None
            chr_pos = snp.get('locations', [{}])[0].get('chromosomePosition') if snp.get('locations') else None
            
            # Fetch detailed association data for this SNP
            assoc_data = gwas.get_associations_by_snp(snp_id)
            
            if assoc_data and '_embedded' in assoc_data:
                assoc_list = assoc_data['_embedded'].get('associations', [])
                for assoc in assoc_list:
                    trait_name = 'Unknown'
                    
                    # Try to get trait from embedded efoTraits first
                    if 'efoTraits' in assoc and len(assoc['efoTraits']) > 0:
                        trait_name = assoc['efoTraits'][0].get('trait', 'Unknown')
                    # Otherwise try to get from betaUnit/description
                    elif 'betaUnit' in assoc and assoc['betaUnit']:
                        trait_name = f"Quantitative trait ({assoc['betaUnit']})"
                    elif 'description' in assoc and assoc['description']:
                        trait_name = assoc['description']
                    
                    pvalue = assoc.get('pvalue')
                    study_name = 'GWAS Catalog'
                    if 'study' in assoc and 'title' in assoc['study']:
                        study_name = assoc['study']['title'][:100]  # Truncate long titles
                    
                    associations.append({
                        'variant_id': snp_id,
                        'chr': chr_name,
                        'pos': chr_pos,
                        'trait': trait_name,
                        'pvalue': pvalue,
                        'study': study_name,
                        'location': 'near_gene'
                    })
            else:
                # Fallback if no detailed associations found
                associations.append({
                    'variant_id': snp_id,
                    'chr': chr_name,
                    'pos': chr_pos,
                    'trait': 'Multiple traits - see GWAS Catalog',
                    'pvalue': None,
                    'study': 'GWAS Catalog',
                    'location': 'near_gene'
                })
    else:
        print("  ⚠ No GWAS associations found by gene name")
    
    # Method 2: Query by genomic region (exact gene boundaries)
    print(f"\n[2/3] Querying GWAS Catalog by region: chr{GENE_CHR}:{GENE_START}-{GENE_END}...")
    data_region = gwas.get_associations_by_region(GENE_CHR, GENE_START, GENE_END)
    
    if data_region and '_embedded' in data_region:
        snps = data_region['_embedded'].get('singleNucleotidePolymorphisms', [])
        print(f"  ✓ Found {len(snps)} SNPs within Gene gene boundaries")
        
        for snp in snps:
            snp_id = snp.get('rsId', '')
            chr_name = snp.get('locations', [{}])[0].get('chromosomeName') if snp.get('locations') else None
            chr_pos = snp.get('locations', [{}])[0].get('chromosomePosition') if snp.get('locations') else None
            
            # Fetch detailed association data for this SNP
            assoc_data = gwas.get_associations_by_snp(snp_id)
            
            if assoc_data and '_embedded' in assoc_data:
                assoc_list = assoc_data['_embedded'].get('associations', [])
                for assoc in assoc_list:
                    trait_name = 'Unknown'
                    
                    # Try to get trait from embedded efoTraits first
                    if 'efoTraits' in assoc and len(assoc['efoTraits']) > 0:
                        trait_name = assoc['efoTraits'][0].get('trait', 'Unknown')
                    # Otherwise try to get from betaUnit/description
                    elif 'betaUnit' in assoc and assoc['betaUnit']:
                        trait_name = f"Quantitative trait ({assoc['betaUnit']})"
                    elif 'description' in assoc and assoc['description']:
                        trait_name = assoc['description']
                    
                    pvalue = assoc.get('pvalue')
                    study_name = 'GWAS Catalog'
                    if 'study' in assoc and 'title' in assoc['study']:
                        study_name = assoc['study']['title'][:100]  # Truncate long titles
                    
                    associations.append({
                        'variant_id': snp_id,
                        'chr': chr_name,
                        'pos': chr_pos,
                        'trait': trait_name,
                        'pvalue': pvalue,
                        'study': study_name,
                        'location': 'within_gene'
                    })
            else:
                # Fallback if no detailed associations found
                associations.append({
                    'variant_id': snp_id,
                    'chr': chr_name,
                    'pos': chr_pos,
                    'trait': 'Multiple traits - see GWAS Catalog',
                    'pvalue': None,
                    'study': 'GWAS Catalog',
                    'location': 'within_gene'
                })
    else:
        print("  ⚠ No GWAS associations found within gene boundaries")
    
    # Save results
    print(f"\n[3/3] Saving results...")
    
    output_file = Path('outputs') / 'gwas_associations.csv'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    df = pd.DataFrame(associations)
    
    # Remove duplicates but keep all info
    df = df.drop_duplicates(subset=['variant_id', 'pos'])
    
    df.to_csv(output_file, index=False)
    print(f"  ✓ GWAS data saved to {output_file}")
    
    # Statistics
    within_gene = df[df['location'] == 'within_gene']
    near_gene = df[df['location'] == 'near_gene']
    
    stats = {
        'total_associations': len(associations),
        'unique_variants': int(df['variant_id'].nunique()) if len(df) > 0 else 0,
        'within_gene': int(len(within_gene)),
        'near_gene': int(len(near_gene)),
    }
    
    stats_file = Path(DATA_PATHS['processed']) / 'gwas_stats.json'
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total Associations:       {stats['total_associations']:,}")
    print(f"Unique Variants:          {stats['unique_variants']:,}")
    print(f"  Within gene:            {stats['within_gene']:,}")
    print(f"  Near gene:              {stats['near_gene']:,}")
    print("=" * 80)
    print("\n✓ GWAS fetching complete!\n")


if __name__ == '__main__':
    main()
