#!/usr/bin/env python3
"""
Script 10: Generate Publication-Ready Figures
Creates various plots for variant analysis.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import GENE_SYMBOL, PROTEIN_LENGTH

def set_aesthetics():
    """Set publication-ready aesthetics for plots."""
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 300,
    })

def plot_domain_map(variants_df, domains):
    """Plot variants along the protein domain map."""
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Protein backbone
    max_len = PROTEIN_LENGTH if PROTEIN_LENGTH > 0 else 500
    ax.plot([0, max_len], [0, 0], color='gray', linewidth=4, zorder=1)
    
    # Domains
    colors = ['#4A90E2', '#50E3C2', '#F5A623', '#D0021B']
    for i, d in enumerate(domains):
        try:
            start = int(d.get('location', {}).get('start', {}).get('value'))
            end = int(d.get('location', {}).get('end', {}).get('value'))
            desc = d.get('description', f'Domain {i+1}')
            
            # Draw domain box
            rect = plt.Rectangle((start, -0.4), end - start, 0.8, 
                                 facecolor=colors[i % len(colors)], 
                                 alpha=0.6, edgecolor='black', zorder=2)
            ax.add_patch(rect)
            
            # Add label
            ax.text((start + end)/2, 0.5, desc, ha='center', va='bottom', 
                    rotation=45 if len(desc) > 15 else 0,
                    fontsize=9)
        except (ValueError, TypeError):
            continue
            
    # Variants subset (e.g., deleterious missense or GWAS)
    if 'protein_position' in variants_df.columns:
        high_impact = variants_df[variants_df['sift'] == 'deleterious']
        gwas = variants_df[variants_df.get('has_gwas', False) == True]
        
        # Plot deleterious
        if not high_impact.empty:
            ax.scatter(high_impact['protein_position'], np.zeros(len(high_impact)) - 0.6, 
                       color='red', marker='v', s=50, label='Deleterious (SIFT)', zorder=3)
            
        # Plot GWAS
        if not gwas.empty and 'protein_position' in gwas.columns:
            gwas_pos = gwas['protein_position'].dropna()
            if not gwas_pos.empty:
                ax.scatter(gwas_pos, np.zeros(len(gwas_pos)) + 0.6, 
                           color='gold', marker='*', s=150, edgecolor='black', 
                           label='GWAS Associated', zorder=4)

    ax.set_xlim(-10, max_len + 10)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([])
    ax.set_xlabel('Amino Acid Position')
    ax.set_title(f'{GENE_SYMBOL} Protein Domain Map & High-Impact Variants')
    ax.legend(loc='upper right')
    
    sns.despine(left=True)
    plt.tight_layout()
    return fig

def plot_eqtl_volcano(eqtl_df):
    """Plot volcano plot for eQTL associations."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if len(eqtl_df) > 0 and 'pvalue' in eqtl_df.columns:
        # Avoid log(0)
        min_pval = eqtl_df[eqtl_df['pvalue'] > 0]['pvalue'].min() 
        eqtl_df['log_pval'] = -np.log10(eqtl_df['pvalue'].clip(lower=min_pval))
        
        # We don't have effect size (NES/beta) in standard summary usually, simulate spreading for visual if missing
        if 'nes' in eqtl_df.columns:
            x_col = 'nes'
            x_label = 'Normalized Effect Size'
        else:
            # Fake spread by using position or random jitter if NES is missing, but better to check
            print("Warning: 'nes' (effect size) not in eQTL data. Plotting by position instead of effect size.")
            # For a proper volcano, let's just use placeholder if NES missing
            x_col = 'pos' if 'pos' in eqtl_df.columns else 'variant_id'
            x_label = 'Position / Variant'
            if x_col == 'variant_id':
                eqtl_df['x_val'] = range(len(eqtl_df))
                x_col = 'x_val'

        sig_threshold = -np.log10(5e-8)
        
        # Plot non-significant
        mask_ns = eqtl_df['log_pval'] < sig_threshold
        ax.scatter(eqtl_df.loc[mask_ns, x_col], eqtl_df.loc[mask_ns, 'log_pval'], 
                   alpha=0.5, color='gray', s=20, label='n.s.')
        
        # Plot significant
        mask_sig = eqtl_df['log_pval'] >= sig_threshold
        ax.scatter(eqtl_df.loc[mask_sig, x_col], eqtl_df.loc[mask_sig, 'log_pval'], 
                   alpha=0.7, color='#D0021B', s=40, label='p < 5e-8')
        
        ax.axhline(sig_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_xlabel(x_label)
        ax.set_ylabel('-log10(P-value)')
        ax.set_title(f'{GENE_SYMBOL} eQTL Associations')
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'No eQTL data available', ha='center', va='center')
        
    sns.despine()
    plt.tight_layout()
    return fig

def plot_af_vs_impact(master_df):
    """Scatter plot of allele frequency vs prediction score."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    has_freq = 'gnomad_af' in master_df.columns or 'gnomad_af_genome' in master_df.columns
    has_impact = 'sift' in master_df.columns
    
    if has_freq and has_impact:
        freq_col = 'gnomad_af' if 'gnomad_af' in master_df.columns else 'gnomad_af_genome'
        
        df_plot = master_df.dropna(subset=[freq_col, 'sift']).copy()
        
        if len(df_plot) > 0:
            # Map SIFT to numeric for y-axis (sift score is usually 0-1 where 0 is deleterious)
            # Since we only have categorical ('deleterious', 'tolerated') in some outputs, handle that
            sift_mapping = {'deleterious': 1, 'tolerated': 0, 'deleterious_low_confidence': 0.8, 'tolerated_low_confidence': 0.2}
            df_plot['impact_score'] = df_plot['sift'].map(sift_mapping).fillna(-0.1)
            
            # Scatter with jitter on Y
            df_plot['y_jitter'] = df_plot['impact_score'] + np.random.normal(0, 0.05, len(df_plot))
            
            # X as log AF
            df_plot['log_af'] = np.log10(df_plot[freq_col].astype(float) + 1e-6)
            
            colors = df_plot['impact_score'].apply(lambda x: '#D0021B' if x >= 0.8 else '#4A90E2')
            
            ax.scatter(df_plot['log_af'], df_plot['y_jitter'], c=colors, alpha=0.6, s=30, edgecolor='w')
            
            ax.set_xlabel('log10(Allele Frequency)')
            ax.set_ylabel('SIFT Impact Category (Jittered)')
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['Tolerated', 'Deleterious'])
            ax.set_title(f'Allele Frequency vs. Predicted Impact')
        else:
            ax.text(0.5, 0.5, 'No overlapping frequency and impact data', ha='center', va='center')
    else:
        ax.text(0.5, 0.5, 'Missing frequency or impact data', ha='center', va='center')
        
    sns.despine()
    plt.tight_layout()
    return fig

def main():
    print("=" * 80)
    print("Generating Publication-Ready Figures")
    print("=" * 80)
    
    set_aesthetics()
    
    output_dir = Path('outputs/figures')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Load Data
    master_file = Path('outputs/master_integrated.csv')
    eqtl_file = Path('outputs/eqtl_associations.csv')
    gene_info_file = Path('data/raw/gene_info.json')
    
    master_df = pd.read_csv(master_file, low_memory=False) if master_file.exists() else pd.DataFrame()
    eqtl_df = pd.read_csv(eqtl_file) if eqtl_file.exists() else pd.DataFrame()
    
    domains = []
    if gene_info_file.exists():
        with open(gene_info_file, 'r') as f:
            gene_data = json.load(f)
            features = gene_data.get('protein', {}).get('features', [])
            domains = [f for f in features if f.get('type') in ['Domain', 'Region']]
            
    # 2. Generate Plots
    
    # Plot 1: Domain Map
    if not master_df.empty and domains:
        print("  Generating Domain Map...")
        fig1 = plot_domain_map(master_df, domains)
        fig1.savefig(output_dir / '1_protein_domain_map.pdf')
        fig1.savefig(output_dir / '1_protein_domain_map.png', dpi=300)
        plt.close(fig1)
        
    # Plot 2: eQTL Volcano
    if not eqtl_df.empty:
        print("  Generating eQTL Volcano Plot...")
        fig2 = plot_eqtl_volcano(eqtl_df)
        fig2.savefig(output_dir / '2_eqtl_volcano.pdf')
        fig2.savefig(output_dir / '2_eqtl_volcano.png', dpi=300)
        plt.close(fig2)
        
    # Plot 3: Frequency vs Impact
    if not master_df.empty:
        print("  Generating Frequency vs Impact Plot...")
        fig3 = plot_af_vs_impact(master_df)
        fig3.savefig(output_dir / '3_frequency_vs_impact.pdf')
        fig3.savefig(output_dir / '3_frequency_vs_impact.png', dpi=300)
        plt.close(fig3)
        
    print(f"\n✓ Figures saved to {output_dir}")

if __name__ == '__main__':
    main()
