# VCF Quality Control Pipeline

A comprehensive pipeline for evaluating VCF file quality from ipyrad assemblies, assessing sample missingness, sequencing depth, and site quality.

## Overview

This pipeline evaluates VCF files by:
- Counting samples and variants
- Calculating individual and site missingness rates
- Measuring mean sequencing depth per sample
- Identifying invariant sites
- Generating comprehensive plots and flagging problematic samples

## Files

- **`vcf_qc.sh`** - Main bash script that orchestrates the entire QC pipeline
- **`vcf_qc.R`** - R script for statistical analysis and visualization (called by bash script)
- **`check.invarSites.sh`** - Script to identify invariant sites (called by bash script)

## Usage

```bash
bash code/vcf_qc.sh input.vcf.gz
```

**Requirements:**
- Modules: `Bioinformatics`, `bcftools`, `vcftools`, `Rtidyverse`
- Input: Compressed VCF file (`.vcf.gz`)

## Directory Structure

```
project/
├── code/
│   ├── vcf_qc.sh           # Main pipeline script
│   ├── vcf_qc.R            # R analysis script
│   └── check.invarSites.sh # Invariant sites checker
├── input.vcf.gz            # Your VCF file
└── output/
    └── vcf_qc/
        ├── vcf_qc.log                    # Pipeline log
        ├── vcf_qc_r.log                  # R script log
        ├── vcf_qc_summary_plots.pdf      # Main visualization
        ├── samples_missing_50percent.csv # Samples >50% missing
        ├── samples_missing_75percent.csv # Samples >75% missing
        ├── samples_missing_80percent.csv # Samples >80% missing
        ├── samples_low_depth.csv         # Samples <5x depth
        ├── indiv_miss.imiss             # Individual missingness
        ├── indiv_depth.idepth           # Individual depth
        └── site_miss.lmiss              # Site missingness
```

## Key Outputs

### Visualization
- **`vcf_qc_summary_plots.pdf`** - Six-panel figure showing:
  - Individual missingness distribution
  - SNP missingness distribution  
  - Individual depth distribution
  - Per-sample missingness plot
  - Top 10 chromosomes with highest SNP missingness
  - Per-sample depth plot

### Sample Flagging
- **`samples_missing_[50|75|80]percent.csv`** - Lists samples exceeding missingness thresholds
- **`samples_low_depth.csv`** - Lists samples with <5x mean sequencing depth

## Quality Control Thresholds

- **Missingness**: 50%, 75%, 80% cutoffs for sample flagging
- **Depth**: <5x mean depth cutoff for low-coverage samples
- **Visual reference**: 50% missingness line on individual plots

## Workflow

1. **`vcf_qc.sh`** loads required modules and sets up output directory
2. Runs `bcftools` and `vcftools` to extract basic statistics
3. Calls **`check.invarSites.sh`** to identify invariant (monomorphic) sites
4. Executes **`vcf_qc.R`** for summary stats and visualization

Use the flagged sample lists to inform filtering decisions for downstream phylogenetic or population genetic analyses.