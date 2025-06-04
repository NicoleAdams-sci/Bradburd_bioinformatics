#!/bin/bash
#################### ipyrad VCF QC ####################
# Evaluate the quality of the VCF from ipyrad: individual sample and site missingness, depth, check for invariant sites

# Usage: code/vcf_qc.sh input.vcf.gz

# Load modules
module load Bioinformatics bcftools vcftools Rtidyverse

# Define inputs
VCF=$1

THREADS=4

# Define output
mkdir -p output/vcf_qc
OUTDIR="output/vcf_qc"

# Write running output to screen and log file
exec > >(tee $OUTDIR/vcf_qc.log) 2>&1  


# Count the number of samples and variants in the ipyrad VCF
echo "--- No. of samples in $VCF ---"
bcftools query -l $VCF | wc -l

echo "--- No. of sites in $VCF ---"
bcftools view -H $VCF --threads $THREADS | wc -l


# Get individual missingness using vcftools
vcftools --gzvcf $VCF --missing-indv --out $OUTDIR/indiv_miss

# Get mean individual mean depth using vcftools
vcftools --gzvcf $VCF --depth --out $OUTDIR/indiv_depth

# Get site missingness using vcftools
vcftools --gzvcf $VCF --missing-site --out $OUTDIR/site_miss

# Check for invariant (monomorphic) sites
echo "--- Checking for invariant sites in $VCF ---"
bash code/check.invarSites.sh $VCF

# Change into output directory
cd $OUTDIR

# Plot missingness and depth distributions in R
Rscript ../../code/vcf_qc.R indiv_miss site_miss indiv_depth


# Get version information for modules used
module list