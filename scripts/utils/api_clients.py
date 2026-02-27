"""API client classes for fetching data from various sources."""

import requests
import time
import json
import os
from pathlib import Path
from typing import Optional, Dict, Any, List
import yaml
from tqdm import tqdm

from .constants import API_URLS, DATA_PATHS


class BaseAPIClient:
    """Base class for API clients with rate limiting and caching."""
    
    def __init__(self, base_url: str, rate_limit: int = 15, cache_dir: str = "data/cache"):
        self.base_url = base_url
        self.rate_limit = rate_limit  # requests per second
        self.min_interval = 1.0 / rate_limit
        self.last_request_time = 0
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'SESN2-Analysis/1.0'})
        
        # Load config
        with open('config.yaml', 'r') as f:
            self.config = yaml.safe_load(f)
    
    def _rate_limit_wait(self):
        """Wait to respect rate limiting."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()
    
    def _get_cache_path(self, cache_key: str) -> Path:
        """Get cache file path for a given key."""
        safe_key = cache_key.replace('/', '_').replace(':', '_')
        return self.cache_dir / f"{safe_key}.json"
    
    def _load_from_cache(self, cache_key: str) -> Optional[Dict]:
        """Load data from cache if available."""
        if not self.config['api']['cache_responses']:
            return None
        
        cache_path = self._get_cache_path(cache_key)
        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Cache read error: {e}")
        return None
    
    def _save_to_cache(self, cache_key: str, data: Dict):
        """Save data to cache."""
        if not self.config['api']['cache_responses']:
            return
        
        cache_path = self._get_cache_path(cache_key)
        try:
            with open(cache_path, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            print(f"Cache write error: {e}")
    
    def get(self, endpoint: str, params: Optional[Dict] = None, use_cache: bool = True) -> Optional[Dict]:
        """Make GET request with rate limiting and retry logic."""
        cache_key = f"{endpoint}_{str(params)}"
        
        # Try cache first
        if use_cache:
            cached = self._load_from_cache(cache_key)
            if cached:
                return cached
        
        url = f"{self.base_url}{endpoint}"
        retry_attempts = self.config['api']['retry_attempts']
        retry_delay = self.config['api']['retry_delay']
        
        for attempt in range(retry_attempts):
            try:
                self._rate_limit_wait()
                response = self.session.get(url, params=params, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    self._save_to_cache(cache_key, data)
                    return data
                elif response.status_code == 429:  # Rate limited
                    wait_time = retry_delay * (2 ** attempt)
                    print(f"Rate limited, waiting {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"Error {response.status_code}: {response.text}")
                    return None
                    
            except requests.exceptions.RequestException as e:
                print(f"Request failed (attempt {attempt + 1}/{retry_attempts}): {e}")
                if attempt < retry_attempts - 1:
                    time.sleep(retry_delay * (2 ** attempt))
        
        return None
    
    def post(self, endpoint: str, data: Dict, params: Optional[Dict] = None, use_cache: bool = True) -> Optional[Dict]:
        """Make POST request with rate limiting and retry logic."""
        cache_key = f"post_{endpoint}_{hash(str(data))}_{str(params)}"
        
        # Try cache first
        if use_cache:
            cached = self._load_from_cache(cache_key)
            if cached:
                return cached
        
        url = f"{self.base_url}{endpoint}"
        retry_attempts = self.config['api']['retry_attempts']
        retry_delay = self.config['api']['retry_delay']
        
        for attempt in range(retry_attempts):
            try:
                self._rate_limit_wait()
                response = self.session.post(url, json=data, params=params, timeout=60)
                
                if response.status_code == 200:
                    result = response.json()
                    self._save_to_cache(cache_key, result)
                    return result
                elif response.status_code == 429:
                    wait_time = retry_delay * (2 ** attempt)
                    print(f"Rate limited, waiting {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    print(f"Error {response.status_code}: {response.text}")
                    return None
                    
            except requests.exceptions.RequestException as e:
                print(f"Request failed (attempt {attempt + 1}/{retry_attempts}): {e}")
                if attempt < retry_attempts - 1:
                    time.sleep(retry_delay * (2 ** attempt))
        
        return None


class EnsemblClient(BaseAPIClient):
    """Client for Ensembl REST API."""
    
    def __init__(self):
        super().__init__(API_URLS['ensembl'])
        self.session.headers.update({'Content-Type': 'application/json'})
    
    def get_gene_info(self, gene_id: str) -> Optional[Dict]:
        """Get gene information."""
        return self.get(f"/lookup/id/{gene_id}")
    
    def get_transcripts(self, gene_id: str) -> Optional[List[Dict]]:
        """Get all transcripts for a gene."""
        data = self.get(f"/lookup/id/{gene_id}", params={'expand': '1'})
        if data and 'Transcript' in data:
            return data['Transcript']
        return None
    
    def get_variants_in_region(self, chr: str, start: int, end: int) -> Optional[List[Dict]]:
        """Get all variants in a genomic region."""
        endpoint = f"/overlap/region/human/{chr}:{start}-{end}"
        return self.get(endpoint, params={'feature': 'variation'})
    
    def get_variant_info(self, variant_id: str) -> Optional[Dict]:
        """Get information about a specific variant."""
        return self.get(f"/variation/human/{variant_id}")
    
    def annotate_variants_vep(self, variants: List[Dict]) -> Optional[List[Dict]]:
        """Annotate variants using VEP. Variants should be in format: {'chr': '1', 'start': 12345, 'allele': 'A/G'}"""
        # VEP accepts HGVS or region format
        vep_input = []
        for v in variants:
            # Format: chr start end allele strand
            vep_input.append(f"{v['chr']} {v['start']} {v['end']} {v['allele']} +")
        
        endpoint = "/vep/human/region"
        params = {"protein": 1, "domains": 1}
        data = {"variants": vep_input}
        return self.post(endpoint, data, params=params)


class GnomADClient:
    """Client for gnomAD GraphQL API."""
    
    def __init__(self):
        self.base_url = API_URLS['gnomad']
        self.session = requests.Session()
        
        # Load config
        with open('config.yaml', 'r') as f:
            self.config = yaml.safe_load(f)
    
    def query_region(self, gene_id: str, chr: str, start: int, end: int) -> Optional[Dict]:
        """Query gnomAD for variants in a region."""
        query = """
        query GnomadVariants($geneId: String!, $datasetId: DatasetId!) {
          gene(gene_id: $geneId, reference_genome: GRCh38) {
            variants(dataset: $datasetId) {
              variant_id
              pos
              ref
              alt
              genome {
                ac
                an
                af
                homozygote_count
                filters
              }
              exome {
                ac
                an
                af
                homozygote_count
                filters
              }
            }
          }
        }
        """
        
        variables = {
            "geneId": gene_id,
            "datasetId": "gnomad_r4"
        }
        
        try:
            response = self.session.post(
                self.base_url,
                json={"query": query, "variables": variables},
                timeout=60
            )
            if response.status_code == 200:
                return response.json()
            else:
                print(f"gnomAD query failed: {response.status_code}")
                return None
        except Exception as e:
            print(f"gnomAD query error: {e}")
            return None


class GTExClient(BaseAPIClient):
    """Client for GTEx Portal API."""
    
    def __init__(self):
        super().__init__(API_URLS['gtex'], rate_limit=10)  # Be conservative
    
    def get_eqtls_for_gene(self, gene_id: str, tissue: str) -> Optional[Dict]:
        """Get eQTL associations for a gene in a specific tissue.
        
        Fixed implementation:
        - Uses GET with query parameters (not POST)
        - Uses correct gene version .4 for GTEx v8 (v26 gencode)
        - Uses correct endpoint /association/singleTissueEqtl
        - Uses correct parameter name 'gencodeId'
        """
        # GTEx v8 uses gencode v26 which has Gene as version .4
        if gene_id == 'ENSG00000130766' and '.' not in gene_id:
            gene_id = 'ENSG00000130766.4'
        elif '.' not in gene_id:
            gene_id = f"{gene_id}.4"  # Default to .4 for GTEx v8
        
        endpoint = "/association/singleTissueEqtl"
        params = {
            'gencodeId': gene_id,
            'tissueSiteDetailId': tissue,
            'datasetId': 'gtex_v8',
            'pageSize': 250  # Get up to 250 eQTLs per tissue
        }
        return self.get(endpoint, params=params)


class GWASCatalogClient(BaseAPIClient):
    """Client for GWAS Catalog API."""
    
    def __init__(self):
        super().__init__(API_URLS['gwas_catalog'], rate_limit=10)
    
    def get_associations_by_gene(self, gene_name: str) -> Optional[Dict]:
        """Get GWAS associations for a gene (returns nearby variants)."""
        endpoint = "/singleNucleotidePolymorphisms/search/findByGene"
        params = {'geneName': gene_name}
        return self.get(endpoint, params=params)
    
    def get_associations_by_region(self, chrom: str, start: int, end: int) -> Optional[Dict]:
        """Get GWAS associations for a genomic region.
        
        Args:
            chrom: Chromosome (e.g., "1", "X")
            start: Start position (bp)
            end: End position (bp)
        """
        endpoint = "/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"
        params = {
            'chrom': chrom,
            'bpStart': start,
            'bpEnd': end
        }
        return self.get(endpoint, params=params)
    
    def get_associations_by_snp(self, rs_id: str) -> Optional[Dict]:
        """Get all GWAS associations for a specific SNP.
        
        Args:
            rs_id: rsID of the variant (e.g., "rs123456")
        
        Returns:
            Dictionary containing association details including traits and p-values
        """
        endpoint = f"/singleNucleotidePolymorphisms/{rs_id}/associations"
        return self.get(endpoint)


class UniProtClient(BaseAPIClient):
    """Client for UniProt REST API."""
    
    def __init__(self):
        super().__init__(API_URLS['uniprot'], rate_limit=10)
    
    def get_protein_info(self, uniprot_id: str) -> Optional[Dict]:
        """Get protein information including domains."""
        endpoint = f"/uniprotkb/{uniprot_id}.json"
        return self.get(endpoint)
    
    def get_protein_features(self, uniprot_id: str) -> Optional[Dict]:
        """Get detailed protein features."""
        return self.get_protein_info(uniprot_id)


class AlphaFoldClient:
    """Client for AlphaFold Database."""
    
    def __init__(self):
        self.base_url = API_URLS['alphafold']
    
    def download_structure(self, uniprot_id: str, output_path: str) -> bool:
        """Download AlphaFold predicted structure."""
        api_url = f"{self.base_url}/api/prediction/{uniprot_id}"
        try:
            api_resp = requests.get(api_url, timeout=30)
            if api_resp.status_code == 200 and len(api_resp.json()) > 0:
                pdb_url = api_resp.json()[-1].get('pdbUrl')
                if not pdb_url:
                    print(f"No PDB URL found in AlphaFold API for {uniprot_id}")
                    return False
                    
                response = requests.get(pdb_url, timeout=30)
                if response.status_code == 200:
                    with open(output_path, 'w') as f:
                        f.write(response.text)
                    return True
                else:
                    print(f"AlphaFold PDB download failed: {response.status_code}")
                    return False
            else:
                print(f"AlphaFold API prediction lookup failed: {api_resp.status_code}")
                return False
        except Exception as e:
            print(f"AlphaFold download error: {e}")
            return False
    
    def get_plddt_scores(self, pdb_file: str) -> Dict[int, float]:
        """Extract pLDDT confidence scores from PDB file."""
        scores = {}
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        # pLDDT is stored in the B-factor column
                        residue_num = int(line[22:26].strip())
                        plddt = float(line[60:66].strip())
                        scores[residue_num] = plddt
        except Exception as e:
            print(f"Error parsing PDB file: {e}")
        return scores


class MyVariantClient(BaseAPIClient):
    """Client for MyVariant.info REST API."""
    
    def __init__(self):
        super().__init__(API_URLS['myvariant'], rate_limit=5)
    
    def get_variants_bulk(self, rsids: List[str]) -> Optional[List[Dict]]:
        """Get variants in bulk from MyVariant.info."""
        endpoint = "/variant"
        
        # Split into batches of 500 (MyVariant limit for POST is 1000, 500 is safer)
        results = []
        batch_size = 500
        for i in range(0, len(rsids), batch_size):
            batch = rsids[i:i+batch_size]
            data = {
                "ids": ",".join(batch),
                "fields": "gnomad_exome,gnomad_genome"
            }
            res = self.post(endpoint, data=data, use_cache=True)
            if res:
                results.extend(res)
            time.sleep(0.5)
        
        return results if results else None
