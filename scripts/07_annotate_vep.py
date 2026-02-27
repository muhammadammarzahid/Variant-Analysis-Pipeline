#!/usr/bin/env python3
"""
Script 07: Annotate Variants with Ensembl VEP
Adds functional predictions and consequence annotations.
"""

import json
import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import EnsemblClient
from utils.constants import DATA_PATHS, CANONICAL_TRANSCRIPT


def main():
    print("=" * 80)
    print("Ensembl VEP Annotator")
    print("=" * 80)
    
    # Load variants
    print(f"\n[1/3] Loading variants...")
    variants_file = Path(DATA_PATHS['processed']) / 'variants_basic.csv'
    
    if not variants_file.exists():
        print(f"  ✗ File not found: {variants_file}")
        return
    
    df = pd.read_csv(variants_file)
    print(f"  ✓ Loaded {len(df)} variants")
    
    # Filter for coding variants only (for detailed annotation)
    coding_variants = df[df['consequence_type'].str.contains('missense|nonsense|frameshift|splice|stop_gained|start_lost', na=False)]
    
    sample_df = coding_variants.copy()
    
    print(f"  Annotating {len(sample_df)} coding variants")
    
    # Annotate with VEP
    print(f"\n[2/3] Querying Ensembl VEP...")
    print("  (This may take 10-20 minutes for large variant sets)")
    
    ensembl = EnsemblClient()
    
    annotations = []
    batch_size = 200
    
    for i in tqdm(range(0, len(sample_df), batch_size), desc="VEP annotation"):
        batch = sample_df.iloc[i:i+batch_size]
        
        # Format for VEP
        vep_input = []
        input_to_id = {}
        for _, row in batch.iterrows():
            # Format: chr start end allele strand
            alleles = row['alleles'].replace('|', '/') if pd.notna(row['alleles']) else 'A/G'
            inp_str = f"{row['chr']} {row['start']} {row['end']} {alleles} +"
            input_to_id[inp_str] = row.get('variant_id')
            vep_input.append({
                'chr': str(row['chr']),
                'start': int(row['start']),
                'end': int(row['end']),
                'allele': alleles
            })
        
        # Query VEP
        vep_result = ensembl.annotate_variants_vep(vep_input)
        
        if vep_result:
            for result in vep_result:
                # Extract most severe consequence
                most_severe = result.get('most_severe_consequence', '')
                
                # Extract transcript consequences
                transcripts = result.get('transcript_consequences', [])
                canonical_transcript = None
                
                for tc in transcripts:
                    if tc.get('transcript_id') == CANONICAL_TRANSCRIPT:
                        canonical_transcript = tc
                        break
                
                if not canonical_transcript and transcripts:
                    canonical_transcript = transcripts[0]
                
                # Extract predictions
                sift = None
                polyphen = None
                if canonical_transcript:
                    sift = canonical_transcript.get('sift_prediction')
                    polyphen = canonical_transcript.get('polyphen_prediction')
                
                inp = result.get('input', '')
                annotations.append({
                    'variant_id': input_to_id.get(inp),
                    'input_coords': inp,
                    'most_severe_consequence': most_severe,
                    'sift': sift,
                    'polyphen': polyphen,
                    'amino_acids': canonical_transcript.get('amino_acids') if canonical_transcript else None,
                    'codons': canonical_transcript.get('codons') if canonical_transcript else None,
                    'protein_start': canonical_transcript.get('protein_start') if canonical_transcript else None,
                    'protein_end': canonical_transcript.get('protein_end') if canonical_transcript else None,
                })
    
    print(f"  ✓ Annotated {len(annotations)} variants")
    
    # Save results
    print(f"\n[3/3] Saving annotations...")
    
    output_file = Path('outputs') / 'variants_annotated.csv'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    annot_df = pd.DataFrame(annotations)
    annot_df.to_csv(output_file, index=False)
    print(f"  ✓ Annotations saved to {output_file}")
    
    # Statistics
    stats = {
        'total_annotated': len(annot_df),
        'with_sift': annot_df['sift'].notna().sum(),
        'sift_deleterious': annot_df['sift'].str.contains('deleterious', na=False).sum() if 'sift' in annot_df else 0,
        'with_polyphen': annot_df['polyphen'].notna().sum(),
        'polyphen_damaging': annot_df['polyphen'].str.contains('damaging', na=False).sum() if 'polyphen' in annot_df else 0,
    }
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Variants Annotated:       {stats['total_annotated']:,}")
    print(f"With SIFT:                {stats['with_sift']:,}")
    print(f"  SIFT Deleterious:       {stats['sift_deleterious']:,}")
    print(f"With PolyPhen:            {stats['with_polyphen']:,}")
    print(f"  PolyPhen Damaging:      {stats['polyphen_damaging']:,}")
    print("=" * 80)
    print("\n✓ VEP annotation complete!\n")


if __name__ == '__main__':
    main()
