#!/usr/bin/env python3
"""
Master Pipeline Runner

Usage:
    python run_pipeline.py --gene <GENE_SYMBOL>

This script fetches the specific genomic and UniProt metrics for the target gene,
updates the universal configuration file (`data/processed/current_gene.json`), and 
sequentially executes all data-fetching and integrating scripts (01 - 11).
"""

import argparse
import json
import os
import subprocess
import requests
from pathlib import Path
import time

def fetch_gene_data(symbol: str):
    print(f"Fetching configuration data for gene: {symbol}...")
    
    # 1. Ensembl Lookup
    ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?expand=1"
    headers = {"Content-Type": "application/json"}
    
    resp = requests.get(ensembl_url, headers=headers)
    if resp.status_code != 200:
        raise ValueError(f"Failed to fetch gene {symbol} from Ensembl: {resp.status_code}")
        
    gene_info = resp.json()
    gene_id = gene_info.get("id")
    chrom = gene_info.get("seq_region_name")
    start = gene_info.get("start")
    end = gene_info.get("end")
    
    # Find canonical transcript
    canonical = None
    transcripts = gene_info.get("Transcript", [])
    for t in transcripts:
        if t.get("is_canonical") == 1:
            canonical = t.get("id")
            break
            
    if not canonical and transcripts:
        # Fallback to longest coding
        coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
        if coding:
            canonical = max(coding, key=lambda x: x.get("end", 0) - x.get("start", 0)).get("id")
        else:
            canonical = transcripts[0].get("id")

    # 2. UniProt Lookup
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{symbol} AND organism_id:9606&format=json"
    u_resp = requests.get(url)
    if u_resp.status_code != 200:
        raise ValueError(f"Failed to fetch {symbol} from UniProt: {u_resp.status_code}")
        
    u_data = u_resp.json()
    results = u_data.get("results", [])
    if not results:
        raise ValueError(f"No UniProt Entry found for gene {symbol}")

    uniprot_entry = results[0]
    uniprot_id = uniprot_entry.get("primaryAccession")
    sequence = uniprot_entry.get("sequence", {})
    protein_length = sequence.get("length", 0)

    # 3. Assemble and write config
    config = {
        "symbol": symbol.upper(),
        "id": gene_id,
        "chr": str(chrom),
        "start": start,
        "end": end,
        "uniprot_id": uniprot_id,
        "canonical_transcript": canonical,
        "protein_length": protein_length
    }
    
    config_dir = Path("data/processed")
    config_dir.mkdir(parents=True, exist_ok=True)
    with open(config_dir / "current_gene.json", "w") as f:
        json.dump(config, f, indent=4)
        
    print(f"✓ Configuration defined perfectly! Targeting {symbol} ({gene_id}) | UniProt: {uniprot_id}")
    return config

def main():
    parser = argparse.ArgumentParser(description="Run the full variant mutation analysis pipeline on a target gene.")
    parser.add_argument("--gene", required=True, help="Official Gene Symbol (e.g. SESN2, TP53)")
    args = parser.parse_args()

    # Create root output dirs if they are destroyed
    for d in ["data/raw", "data/processed", "data/cache", "outputs/figures", "structures"]:
        Path(d).mkdir(parents=True, exist_ok=True)

    # Initialize generic target configs based on provided symbol
    fetch_gene_data(args.gene.upper())
    
    # Target Scripts exactly in sequence
    scripts = [
        "scripts/01_fetch_gene_info.py",
        "scripts/02_fetch_variants.py",
        "scripts/03_fetch_gnomad.py",
        "scripts/04_fetch_eqtl.py",
        "scripts/05_fetch_pqtl.py",
        "scripts/06_fetch_gwas.py",
        "scripts/07_annotate_vep.py",
        "scripts/08_map_protein_domains.py",
        "scripts/09_integrate_all.py",
        "scripts/10_generate_figures.py",
        "scripts/11_visualize_structure.py"
    ]
    
    print("\n" + "="*80)
    print(f"STARTING PIPELINE EXECUTION FOR {args.gene.upper()}")
    print("="*80 + "\n")
    
    for script_file in scripts:
        print(f"\n🚀 Running {script_file}...")
        
        # Execute script and pipe outputs live
        result = subprocess.run(
            ["pixi", "run", "python", script_file], 
            text=True
        )
        
        if result.returncode != 0:
            print(f"\n❌ FAILED at {script_file}. Pipeline aborted.")
            return
            
    print("\n" + "="*80)
    print(f"PIPELINE COMPLETE - 100% DONE FOR {args.gene.upper()}")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
