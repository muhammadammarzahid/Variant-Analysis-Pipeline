# Variant Mutation and QTL Analysis Pipeline

Comprehensive dynamic analysis pipeline for gathering and analyzing mutations, eQTLs, pQTLs, and GWAS associations for *any* target human gene (e.g., SESN2, PTEN).

## Overview

This pipeline systematically collects and integrates data from multiple free APIs to analyze:
- All known genetic variants in the target gene region
- Expression QTLs (eQTLs) from GTEx
- Protein QTLs (pQTLs) from available databases
- GWAS associations and phenotype links
- Functional predictions and structural impacts (SIFT / PolyPhen)
- Protein domain mapping for coding variants
- Automated 3D Structural PyMOL Mapping from AlphaFold

## Project Structure

```
Sesn2_mutations/
├── scripts/                      # Analysis scripts
│   ├── 01_fetch_gene_info.py     # Gene/protein metadata
│   ├── 02_fetch_variants.py      # Ensembl variants
│   ├── 03_fetch_gnomad.py        # Population frequencies
│   ├── 04_fetch_eqtl.py          # GTEx eQTL data
│   ├── 05_fetch_pqtl.py          # pQTL associations
│   ├── 06_fetch_gwas.py          # GWAS catalog
│   ├── 07_annotate_vep.py        # VEP annotations
│   ├── 08_map_protein_domains.py # Domain mapping
│   ├── 09_integrate_all.py       # Integration & summary
│   └── utils/                    # Utility modules
├── data/                         # Raw and processed data
│   ├── raw/                      # Raw API responses
│   ├── processed/                # Cleaned datasets
│   └── cache/                    # API response cache
├── outputs/                      # Final results
│   ├── master_integrated.csv     # All variants with annotations
│   ├── eqtl_associations.csv     # eQTL results
│   ├── gwas_associations.csv     # GWAS hits
│   ├── domain_mapping.csv        # Variant-domain mapping
│   ├── variants_annotated.csv    # VEP annotations
│   └── summary_report.txt        # Human-readable summary
├── config.yaml                   # Configuration file
├── pixi.toml                     # Package management
└── README.md                     # This file
```

## Installation

This project uses [pixi](https://prefix.dev/) for dependency management.

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Install dependencies
pixi install
```

## Usage

### Run Entire Pipeline Dynamically

The pipeline is managed by a master orchestrator (`run_pipeline.py`) that accepts any gene symbol. It automatically queries Ensembl and UniProt for dynamic coordinates, canonical transcripts, and IDs, then sequentially executes all scripts (`01`-`11`).

```bash
pixi run python run_pipeline.py --gene PTEN
```

This seamlessly executes all data collection, annotation, integration, and graphical rendering steps.

### Run Individual Steps

```bash
# 1. Fetch gene information
pixi run python scripts/01_fetch_gene_info.py

# 2. Fetch variants
pixi run python scripts/02_fetch_variants.py

# 3. Fetch population frequencies (via MyVariant.info)
pixi run python scripts/03_fetch_gnomad.py

# 4. Fetch eQTL data (takes 10-20 minutes)
pixi run python scripts/04_fetch_eqtl.py

# 5. Fetch pQTL data (placeholder - manual curation recommended)
pixi run python scripts/05_fetch_pqtl.py

# 6. Fetch GWAS associations
pixi run python scripts/06_fetch_gwas.py

# 7. Annotate with VEP (takes 10-20 minutes)
pixi run python scripts/07_annotate_vep.py

# 8. Map to protein domains
pixi run python scripts/08_map_protein_domains.py

# 9. Integrate all data
pixi run python scripts/09_integrate_all.py
```

## Data Sources

| Source | Data Type | API Endpoint |
|--------|-----------|--------------|
| Ensembl REST API | Variants, Gene Info, VEP | https://rest.ensembl.org |
| UniProt | Protein Domains, Features | https://rest.uniprot.org |
| GTEx Portal | eQTL Associations | https://gtexportal.org/api/v2 |
| GWAS Catalog | GWAS Associations | https://www.ebi.ac.uk/gwas/rest/api |
| MyVariant.info | gnomAD Frequencies | https://myvariant.info/v1/variant |
| AlphaFold DB | Predicted Structure | https://alphafold.ebi.ac.uk |

## Output Files

### Main Output: master_integrated.csv

Contains all variants with integrated annotations:
- Variant ID, position, alleles
- Consequence types
- Population frequencies (if available)
- eQTL associations (tissue, p-value)
- GWAS associations (traits, p-value)
- Clinical significance
- Domain information

### Additional Outputs

1. **eqtl_associations.csv**: Tissue-specific expression QTLs
2. **gwas_associations.csv**: GWAS phenotype associations
3. **domain_mapping.csv**: Coding variants mapped to protein domains
4. **variants_annotated.csv**: VEP functional predictions (SIFT, PolyPhen)
5. **summary_report.txt**: Human-readable analysis summary
6. **summary_stats.json**: Summary statistics in JSON format



## Limitations

1. **gnomAD / MyVariant**: Relies on third-party Elasticsearch aggregation; novel or highly rare variants might occasionally drop.
2. **pQTL**: Limited free API access; manual curation from literature recommended
3. **eQTL**: GTEx API may have rate limits; full analysis takes time
4. **VEP**: Sample subset annotated due to API rate limits
5. **AlphaFold**: Structure download may fail (URL changes); manual download available



## Configuration

The pipeline's exact biological thresholds and API behaviors are fully exposed and controlled via `config.yaml`:

### Biological Thresholds
```yaml
thresholds:
  gnomad_af_rare: 0.01        # Allele Frequency (AF) defining "Rare" variants
  gnomad_af_common: 0.05      # AF bounding "Common" variants
  eqtl_pvalue: 0.00001        # Minimum p-value to associate a regulatory eQTL 
  gwas_pvalue: 0.00000005     # Standard genome-wide GWAS significance (< 5x10^-8)
  sift_deleterious: 0.05      # SIFT score ceiling indicating deleterious mutations
  polyphen_damaging: 0.45     # PolyPhen score floor establishing damaging impact
```

### API Behavior and Reliability
```yaml
api:
  rate_limit: 15              # Global requests per second to avoid API ban
  retry_attempts: 3           # Exponential backoff auto-recovery max attempts
  retry_delay: 2              # Base wait seconds during an HTTP 429 Rate Limit
  cache_responses: true       # Prevents duplicate external calls to speed up reruns
```

### Visual & Output Options
```yaml
output:
  include_non_coding: true        # Track 3' UTR, 5' UTR, and Intron associations
  min_confidence_domains: 50      # PyMOL AlphaFold pLDDT visualization baseline
  max_variants_display: 10        # Limit maximum priority points plotted structurally
```

*(Note: The `sesn2` gene mapping block at the top of the raw `config.yaml` is kept only as a legacy fallback. It is automatically overridden by `run_pipeline.py` during execution).*

## Troubleshooting

### API Rate Limits
- The pipeline implements automatic rate limiting and retry logic
- Responses are cached to avoid re-querying
- If you encounter 429 errors, wait a few minutes and retry

### Missing Data
- Some APIs may be temporarily unavailable
- Scripts will continue with available data
- Check individual output files for completeness

### Memory Issues
- Large variant sets may require substantial memory
- Consider filtering variants by frequency or consequence type
- Process in batches if needed

## Citation

If you use this pipeline, please cite the data sources:
- Ensembl: https://www.ensembl.org
- GTEx: https://gtexportal.org
- GWAS Catalog: https://www.ebi.ac.uk/gwas
- UniProt: https://www.uniprot.org
- MyVariant: https://myvariant.info

## License

This project is open-source and licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues with the pipeline, please refer to the individual API documentation links above.

---

**Analysis Date**: Generated dynamically
**Pipeline Version**: 1.0
**Target Gene**: Variable (Provided at runtime via `run_pipeline.py`)
