#!/usr/bin/env python3
"""
Script 04: Fetch eQTL Data from GTEx
Retrieves expression QTL associations for Gene across tissues.
"""

import json
import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import GTExClient
from utils.constants import GENE_ID, GTEX_TISSUES, DATA_PATHS


def main():
    print("=" * 80)
    print("GTEx eQTL Fetcher")
    print("=" * 80)
    
    gtex = GTExClient()
    
    print(f"\n[1/2] Querying GTEx for {GENE_ID} across {len(GTEX_TISSUES)} tissues...")
    print("  (This may take 5-10 minutes due to API rate limiting)")
    
    all_eqtls = []
    tissues_with_eqtls = 0
    
    for tissue in tqdm(GTEX_TISSUES, desc="Querying tissues"):
        data = gtex.get_eqtls_for_gene(GENE_ID, tissue)
        
        if data and 'data' in data:
            eqtls = data['data']
            if eqtls and len(eqtls) > 0:
                tissues_with_eqtls += 1
                for eqtl in eqtls:
                    all_eqtls.append({
                        'variant_id': eqtl.get('variantId', ''),  # GTEx standard format
                        'rsid': eqtl.get('snpId', ''),            # rsID for proper merging
                        'gene_id': eqtl.get('geneId', ''),
                        'tissue': tissue,
                        'pvalue': eqtl.get('pValue'),
                        'nes': eqtl.get('nes'),  # normalized effect size
                    })
    
    print(f"\n  ✓ Found eQTLs in {tissues_with_eqtls} tissues")
    print(f"  ✓ Total eQTL associations: {len(all_eqtls)}")
    
    # Save results
    print(f"\n[2/2] Saving results...")
    
    if all_eqtls:
        df = pd.DataFrame(all_eqtls)
        
        output_dir = Path(DATA_PATHS['raw'])
        output_dir.mkdir(parents=True, exist_ok=True)
        
        csv_file = Path('outputs') / 'eqtl_associations.csv'
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv_file, index=False)
        print(f"  ✓ eQTL data saved to {csv_file}")
        
        # Statistics
        stats = {
            'total_eqtls': int(len(df)),
            'unique_variants': int(df['variant_id'].nunique()),
            'tissues_with_eqtls': int(tissues_with_eqtls),
            'significant_eqtls_5e8': int((df['pvalue'] < 5e-8).sum()),
            'top_tissues': {k: int(v) for k, v in df['tissue'].value_counts().head(10).to_dict().items()},
        }
        
        stats_file = Path(DATA_PATHS['processed']) / 'eqtl_stats.json'
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(f"Total eQTL Associations:  {stats['total_eqtls']:,}")
        print(f"Unique Variants:          {stats['unique_variants']:,}")
        print(f"Tissues with eQTLs:       {stats['tissues_with_eqtls']}/{len(GTEX_TISSUES)}")
        print(f"Genome-wide Significant:  {stats['significant_eqtls_5e8']:,}")
        print("=" * 80)
    else:
        print("  ⚠ No eQTL data found")
    
    print("\n✓ eQTL fetching complete!\n")


if __name__ == '__main__':
    main()
