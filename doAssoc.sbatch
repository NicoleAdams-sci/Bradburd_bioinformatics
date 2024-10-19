#!/bin/bash
#SBATCH --job-name assoc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
##SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output assoc.%j.out
##SBATCH --error assoc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=username@institution.edu

# Script to do GEA on low cov data using ANGSD -doAsso. Input files are 1. a beagle file 2. environmental data file (e.g. precip for each sample lat/long) 3. covariate file (e.g. PCs from PCAngsd to control for population structure)

## submit this script like: for ENV in `ls gsb.*2000_08.txt | grep -v "PRISM"`; do OUTNAM=$(echo $ENV| cut -f2 -d"."); echo $ENV; sbatch ~/pman/scripts/doAssoc.sbatch pman.noLow.dwnSamp.noX.beagle.gprobs_2000_08.noUn.beagle.gz gsb.genPCs.imput.txt $ENV $OUTNAM; done

date
time

module load Bioinformatics angsd

IMP_BEAG=$1
COV=$2
ENV=$3
OUTNAM=$4


BEAGDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp"
OUTDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp/gea"
FAI="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna_noUn.fai"
# Make sure the .fai file has no chromosomes with '_' ** this is only for imputed beagle

# Run association with structure covariants
# doAsso 4 = latent genotype model
# Pvalue 1 = prints P-value instead of a likelihood ratio (0 = no -default)
# cov = covariate file
# doMaf 4 = Freq from genotype probabilities

angsd -beagle $BEAGDIR/$IMP_BEAG -fai $FAI \
-P 5 -doMaf 4 -doAsso 4 -Pvalue 1 -yQuant $ENV -cov $COV -out $OUTDIR/$OUTNAM\.assoc -nThreads ${SLURM_NTASKS_PER_NODE}