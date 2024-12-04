#!/bin/bash
#SBATCH --job-name assoc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
##SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output assoc.%j.out
##SBATCH --error assoc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=username@institution.edu

# Script to do GEA on low cov data using ANGSD -doAsso. Input files are 1. a beagle file 2. a sample file that includes environmental data (e.g. precip for each sample lat/long) and covariate info (e.g. PCs from PCAngsd to control for population structure)
# for sample file format see https://www.chg.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html

## submit this script like: sbatch ~/pman/scripts/doAssoc_rBay.sbatch pman.noLow.dwnSamp.noX.beagle.gprobs_2000_08.noUn.beagle.gz sampleFile.txt $OUTNAM

date
time

module load Bioinformatics angsd

IMP_BEAG=$1
SAMPFILE=$2
OUTNAM=$4


BEAGDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp"
OUTDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp/gea"
FAI="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna_noUn.fai"
# Make sure the .fai file has no chromosomes with '_' ** this is only for imputed beagle

# Run association with structure covariants
# doAsso 4 = latent genotype model (will give effect sizes, depending on sample size)
# doAsso 2 = Score Test (will not give effect sizes)
# Pvalue 1 = prints P-value instead of a likelihood ratio (0 = no -default)
# whichPhe = comma separated list of "phenotypes" (environ variables) to test from sampleFile.txt
# whichCov = comma separated list of covariates from sampleFile.txt
# doMaf 4 = Freq from genotype probabilities

angsd -beagle $BEAGDIR/$IMP_BEAG -fai $FAI \
-beagle $BEAGF \
-P ${SLURM_NTASKS_PER_NODE} \
-fai $FAI -doAsso 2 -doMaf 4 -Pvalue 1 -minMaf 0.1 \
-sampleFile sampleFile.txt \
-whichPhe ENV1,ENV2 \
-whichCov Lib,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
-out $OUTDIR/HETH.contemp.male.PCcovar