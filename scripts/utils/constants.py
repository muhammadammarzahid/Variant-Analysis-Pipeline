"""Constants for Gene analysis."""

import json
from pathlib import Path

# Load dynamic gene configuration
current_gene_file = Path('data/processed/current_gene.json')
if current_gene_file.exists():
    with open(current_gene_file, 'r') as f:
        gene_data = json.load(f)
else:
    # Fallback to SESN2 defaults if running out of the box without the runner script
    gene_data = {
        'symbol': 'SESN2',
        'id': 'ENSG00000130766',
        'chr': '1',
        'start': 28259473,
        'end': 28282491,
        'uniprot_id': 'P58004',
        'canonical_transcript': 'ENST00000253063',
        'protein_length': 480
    }

# Gene configuration
GENE_ID = gene_data.get('id', '')
UNIPROT_ID = gene_data.get('uniprot_id', '')
GENE_SYMBOL = gene_data.get('symbol', 'UNKNOWN')
GENE_CHR = str(gene_data.get('chr', ''))
GENE_START = int(gene_data.get('start', 0))
GENE_END = int(gene_data.get('end', 0))
CANONICAL_TRANSCRIPT = gene_data.get('canonical_transcript', '')
PROTEIN_LENGTH = int(gene_data.get('protein_length', 0))

# API Endpoints
API_URLS = {
    'ensembl': 'https://rest.ensembl.org',
    'gnomad': 'https://gnomad.broadinstitute.org/api',
    'gtex': 'https://gtexportal.org/api/v2',
    'gwas_catalog': 'https://www.ebi.ac.uk/gwas/rest/api',
    'opengwas': 'https://gwas.mrcieu.ac.uk/api',
    'phenoscanner': 'http://www.phenoscanner.medschl.cam.ac.uk/api',
    'uniprot': 'https://rest.uniprot.org',
    'alphafold': 'https://alphafold.ebi.ac.uk',
    'opentargets': 'https://api.genetics.opentargets.org/graphql',
    'myvariant': 'https://myvariant.info/v1',
}

# GTEx Tissues (all 54)
GTEX_TISSUES = [
    'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum',
    'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial',
    'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum',
    'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus',
    'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1',
    'Brain_Substantia_nigra', 'Breast_Mammary_Tissue',
    'Cells_Cultured_fibroblasts', 'Cells_EBV-transformed_lymphocytes',
    'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction',
    'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Fallopian_Tube',
    'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Kidney_Cortex',
    'Kidney_Medulla', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal',
    'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate',
    'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg',
    'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis',
    'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'
]

# Variant consequence severity
CONSEQUENCE_SEVERITY = {
    'high': [
        'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
        'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost'
    ],
    'moderate': [
        'inframe_insertion', 'inframe_deletion', 'missense_variant',
        'protein_altering_variant'
    ],
    'low': [
        'splice_region_variant', 'incomplete_terminal_codon_variant',
        'start_retained_variant', 'stop_retained_variant', 'synonymous_variant'
    ],
    'modifier': [
        '5_prime_UTR_variant', '3_prime_UTR_variant', 'intron_variant',
        'intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
        'non_coding_transcript_variant', 'regulatory_region_variant'
    ]
}

# Data file paths
DATA_PATHS = {
    'raw': 'data/raw',
    'processed': 'data/processed',
    'cache': 'data/cache',
    'outputs': 'outputs'
}
