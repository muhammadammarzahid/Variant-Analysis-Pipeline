#!/usr/bin/env python3
"""
Script 01: Fetch Gene Gene and Protein Information
Retrieves basic gene metadata from Ensembl and UniProt.
"""

import json
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.api_clients import EnsemblClient, UniProtClient, AlphaFoldClient
from utils.constants import (
    GENE_ID, UNIPROT_ID, GENE_SYMBOL,
    GENE_CHR, GENE_START, GENE_END, DATA_PATHS
)


def main():
    print("=" * 80)
    print("Gene Gene and Protein Information Fetcher")
    print("=" * 80)
    
    # Initialize clients
    ensembl = EnsemblClient()
    uniprot = UniProtClient()
    alphafold = AlphaFoldClient()
    
    results = {
        'gene': None,
        'transcripts': None,
        'protein': None,
        'structure': None
    }
    
    # 1. Fetch gene information from Ensembl
    print(f"\n[1/4] Fetching gene information for {GENE_SYMBOL} ({GENE_ID})...")
    gene_info = ensembl.get_gene_info(GENE_ID)
    if gene_info:
        results['gene'] = gene_info
        print(f"  ✓ Gene: {gene_info.get('display_name')}")
        print(f"  ✓ Location: chr{gene_info.get('seq_region_name')}:{gene_info.get('start')}-{gene_info.get('end')}")
        print(f"  ✓ Strand: {'+' if gene_info.get('strand') == 1 else '-'}")
        print(f"  ✓ Biotype: {gene_info.get('biotype')}")
    else:
        print("  ✗ Failed to fetch gene information")
    
    # 2. Fetch transcript information
    print(f"\n[2/4] Fetching transcript information...")
    transcripts = ensembl.get_transcripts(GENE_ID)
    if transcripts:
        results['transcripts'] = transcripts
        print(f"  ✓ Found {len(transcripts)} transcripts")
        for t in transcripts[:3]:  # Show first 3
            print(f"    - {t.get('id')}: {t.get('biotype')} ({t.get('length')} bp)")
        if len(transcripts) > 3:
            print(f"    ... and {len(transcripts) - 3} more")
    else:
        print("  ✗ Failed to fetch transcript information")
    
    # 3. Fetch protein information from UniProt
    print(f"\n[3/4] Fetching protein information from UniProt ({UNIPROT_ID})...")
    protein_info = uniprot.get_protein_info(UNIPROT_ID)
    if protein_info:
        results['protein'] = protein_info
        
        # Extract key information
        seq_length = protein_info.get('sequence', {}).get('length', 'Unknown')
        print(f"  ✓ Protein: {protein_info.get('uniProtkbId')}")
        print(f"  ✓ Length: {seq_length} amino acids")
        
        # Extract domains
        features = protein_info.get('features', [])
        domains = [f for f in features if f.get('type') in ['Domain', 'Region']]
        if domains:
            print(f"  ✓ Found {len(domains)} domains/regions:")
            for domain in domains[:5]:
                location = domain.get('location', {})
                start = location.get('start', {}).get('value', '?')
                end = location.get('end', {}).get('value', '?')
                desc = domain.get('description', 'Unknown')
                print(f"    - {desc} ({start}-{end})")
        
        # Extract function
        comments = protein_info.get('comments', [])
        function_comments = [c for c in comments if c.get('commentType') == 'FUNCTION']
        if function_comments:
            print(f"  ✓ Function:")
            for fc in function_comments[:1]:  # Show first function comment
                texts = fc.get('texts', [])
                if texts:
                    func_text = texts[0].get('value', '')
                    # Truncate long text
                    if len(func_text) > 200:
                        func_text = func_text[:197] + '...'
                    print(f"    {func_text}")
    else:
        print("  ✗ Failed to fetch protein information")
    
    # 4. Download AlphaFold structure
    print(f"\n[4/4] Downloading AlphaFold structure...")
    structure_path = Path(DATA_PATHS['raw']) / f"AF-{UNIPROT_ID}-F1-model_v4.pdb"
    success = alphafold.download_structure(UNIPROT_ID, str(structure_path))
    if success:
        print(f"  ✓ Structure saved to {structure_path}")
        results['structure'] = str(structure_path)
        
        # Extract pLDDT scores
        plddt_scores = alphafold.get_plddt_scores(str(structure_path))
        if plddt_scores:
            avg_plddt = sum(plddt_scores.values()) / len(plddt_scores)
            print(f"  ✓ Average pLDDT confidence: {avg_plddt:.1f}")
            results['avg_plddt'] = avg_plddt
    else:
        print("  ✗ Failed to download structure")
    
    # Save raw results
    print(f"\n[Saving] Writing results to disk...")
    raw_dir = Path(DATA_PATHS['raw'])
    raw_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = raw_dir / 'gene_info.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"  ✓ Raw data saved to {output_file}")
    
    # Create summary file
    summary = {
        'gene_id': GENE_ID,
        'gene_symbol': GENE_SYMBOL,
        'chromosome': gene_info.get('seq_region_name') if gene_info else GENE_CHR,
        'start': gene_info.get('start') if gene_info else GENE_START,
        'end': gene_info.get('end') if gene_info else GENE_END,
        'strand': gene_info.get('strand') if gene_info else None,
        'gene_length': (gene_info.get('end') - gene_info.get('start') + 1) if gene_info else None,
        'num_transcripts': len(transcripts) if transcripts else 0,
        'canonical_transcript': gene_info.get('canonical_transcript') if gene_info else None,
        'uniprot_id': UNIPROT_ID,
        'protein_length': protein_info.get('sequence', {}).get('length') if protein_info else None,
        'num_domains': len([f for f in protein_info.get('features', []) if f.get('type') in ['Domain', 'Region']]) if protein_info else 0,
        'alphafold_structure': str(structure_path) if success else None,
        'avg_plddt': results.get('avg_plddt'),
    }
    
    processed_dir = Path(DATA_PATHS['processed'])
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    summary_file = processed_dir / 'gene_metadata.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  ✓ Summary saved to {summary_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Gene Symbol:          {summary['gene_symbol']}")
    print(f"Gene ID:              {summary['gene_id']}")
    print(f"Location:             chr{summary['chromosome']}:{summary['start']}-{summary['end']}")
    print(f"Gene Length:          {summary['gene_length']:,} bp" if summary['gene_length'] else "Gene Length:          Unknown")
    print(f"Transcripts:          {summary['num_transcripts']}")
    print(f"Protein Length:       {summary['protein_length']} aa" if summary['protein_length'] else "Protein Length:       Unknown")
    print(f"Protein Domains:      {summary['num_domains']}")
    print(f"AlphaFold Structure:  {'Available' if summary['alphafold_structure'] else 'Not downloaded'}")
    if summary.get('avg_plddt'):
        print(f"Avg. Structure Conf:  {summary['avg_plddt']:.1f}")
    print("=" * 80)
    print("\n✓ Gene information collection complete!\n")


if __name__ == '__main__':
    main()
