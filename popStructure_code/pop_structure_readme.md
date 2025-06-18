# Population Structure Analysis Pipeline

## Overview
This bash script performs comprehensive population structure analysis on genomic VCF data using multiple approaches including PCA and ADMIXTURE clustering.

## Usage
```bash
bash pop_structure_all_samples.sh input.vcf.gz metadata [imiss] [lmiss] [hwp] [k] [threads]
```

### Parameters
- `input.vcf.gz`: Input VCF file (required)
- `metadata`: Sample metadata file (required)
- `imiss`: Individual missingness threshold (default: 0.5)
- `lmiss`: Locus missingness threshold (default: 0.5)
- `hwp`: Hardy-Weinberg equilibrium p-value cutoff (default: 0.01)
- `k`: Maximum number of clusters for ADMIXTURE (default: 5)
- `threads`: Number of processing threads (default: 4)

## Analysis Steps

### 1. VCF Quality Filtering
- Applies basic quality filters (minimum GQ=15, DP=10)
- Removes indels and keeps only biallelic sites
- Filters samples and loci based on missingness thresholds
- Applies minor allele frequency filter (MAF > 0.05)

### 2. Chromosome Naming (for Plink)
- Checks if chromosomes are numeric
- Adds "scaffold_" prefix to numeric chromosomes for compatibility

### 3. Linkage Disequilibrium Pruning
- Uses PLINK to identify sites in linkage disequilibrium
- Removes correlated sites (window size: 50, step: 5, r² > 0.5)
- Applies Hardy-Weinberg equilibrium filtering

### 4. Principal Component Analysis (PCA)
- **PLINK PCA**: Standard PCA implementation with R plotting
- **VCF2PCACluster**: Alternative PCA method using Bradburd module

### 5. Population Structure Analysis
- **ADMIXTURE**: Tests K=2 to specified maximum K with cross-validation
- Generates admixture plots and cross-validation results
- *(FastStructure section commented out but available)*

## Output Structure
```
output/pop_struct_all/
├── filtered VCF files
├── plink_results/          # PLINK PCA results
├── plink_pca_results/      # PCA plots
├── vcf2pca_results/        # Alternative PCA results
└── admixture_results/      # ADMIXTURE clustering results
```

## Dependencies
- bcftools, vcftools, htslib
- PLINK
- R with tidyverse
- VCF2PCACluster (Bradburd module)
- ADMIXTURE (Bradburd module)

## Output Files
- Filtered VCF: `pop_all.filtered.vcf.gz`
- PCA results and plots
- ADMIXTURE Q-matrices and cross-validation results
- Log file: `pop_structure_all_samples.log`