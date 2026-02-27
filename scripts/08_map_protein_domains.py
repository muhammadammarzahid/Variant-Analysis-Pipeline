#!/usr/bin/env python3
"""
Script 08: Map Variants to Protein Domains
Maps coding variants to protein domains and structural features.
"""

import json
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from utils.constants import DATA_PATHS

def main():
    print("=" * 80)
    print("Protein Domain Mapper")
    print("=" * 80)
    
    # Load gene info with protein domains
    print(f"\n[1/3] Loading protein domain information...")
    gene_info_file = Path(DATA_PATHS['raw']) / 'gene_info.json'
    
    if not gene_info_file.exists():
        print("  ✗ Gene info not found. Run 01_fetch_gene_info.py first")
        return
    
    with open(gene_info_file, 'r') as f:
        gene_data = json.load(f)
    
    protein_info = gene_data.get('protein', {})
    features = protein_info.get('features', [])
    domains = [f for f in features if f.get('type') in ['Domain', 'Region']]
    
    print(f"  ✓ Found {len(domains)} domains/regions")
    for d in domains[:5]:
        loc = d.get('location', {})
        start = loc.get('start', {}).get('value', '?')
        end = loc.get('end', {}).get('value', '?')
        desc = d.get('description', 'Unknown')
        print(f"    - {desc} ({start}-{end})")
    
    # Load variants
    print(f"\n[2/3] Loading variants...")
    annot_file = Path('outputs') / 'variants_annotated.csv'
    
    if not annot_file.exists():
        print("  ✗ Annotated variants file not found. Run 07_annotate_vep.py first")
        return
    
    df = pd.read_csv(annot_file)
    
    # Filter for coding variants that have protein positions
    coding_df = df.dropna(subset=['protein_start'])
    
    print(f"  ✓ Found {len(coding_df)} coding variants with protein mapping")
    
    # Map to domains
    print(f"\n[3/3] Creating domain mapping...")
    
    domain_mapping = []
    
    for _, row in coding_df.iterrows():
        try:
            prot_start = int(float(row['protein_start']))
            prot_end = int(float(row['protein_end'])) if not pd.isna(row['protein_end']) else prot_start
        except ValueError:
            continue
            
        # Find affected domains
        affected_domains = []
        domain_starts = []
        domain_ends = []
        
        for d in domains:
            loc = d.get('location', {})
            d_start_val = loc.get('start', {}).get('value')
            d_end_val = loc.get('end', {}).get('value')
            
            if d_start_val and d_end_val:
                try:
                    d_start = int(d_start_val)
                    d_end = int(d_end_val)
                    if (prot_start <= d_end and prot_end >= d_start):
                        affected_domains.append(d.get('description', 'Unknown'))
                        domain_starts.append(str(d_start))
                        domain_ends.append(str(d_end))
                except ValueError:
                    continue
        
        domain_mapping.append({
            'variant_id': row.get('variant_id'),
            'input_coords': row.get('input_coords'),
            'consequence': row.get('most_severe_consequence'),
            'protein_position': prot_start,
            'amino_acid_change': row.get('amino_acids'),
            'codon_change': row.get('codons'),
            'domain_affected': '|'.join(affected_domains) if affected_domains else 'None',
            'domain_start': '|'.join(domain_starts) if domain_starts else None,
            'domain_end': '|'.join(domain_ends) if domain_ends else None,
        })
    
    mapping_df = pd.DataFrame(domain_mapping)
    
    # Save results
    output_file = Path('outputs') / 'variant_protein_mapping.csv'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    mapping_df.to_csv(output_file, index=False)
    print(f"  ✓ Domain mapping saved to {output_file}")
    
    # Group by domain for summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Mapped Coding Variants:    {len(mapping_df):,}")
    print(f"Known Protein Domains:    {len(domains)}\n")
    
    print(f"Variants by Domain:")
    domain_counts = mapping_df['domain_affected'].value_counts()
    for domain, count in domain_counts.items():
        print(f"  {domain}: {count}")
    print("=" * 80)
    print("\n✓ Domain mapping complete!\n")

if __name__ == '__main__':
    main()
