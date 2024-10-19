#!/bin/bash
#SBATCH --job-name annot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard
#SBATCH --output annot.%j.out
##SBATCH --error annot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=username@institution.edu

# Script to get annotations using LD-annot after reformatting GFF and a VCF 
# Need a VCF of your beagle file to do this, run angsd with "-doBcf 1" to get it
# Modified from Rachel Bay's annot.sbatch for HETH

# Download and configure LD-annot prior to running script
#wget https://github.com/ArnaudDroitLab/LD-annot/archive/refs/heads/master.zip
#unzip master.zip
#cd LD-annot-master/
#chmod +x LD-annot0.4.py
#unzip plink_linux_x86_64_20190304.zip


## submit this script like: sbatch ~/pman/scripts/annot.sbatch

date
time

module load Bioinformatics plink/1.9 bcftools/1.12-g4b275e vcftools/0.1.15 # LD-annot can only use plink 1.9

REFDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/"
DIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp/gea/doAssoc_impute/annot/LD-annot-master"
VCF="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp/pman.noLow.dwnSamp.vcf.gz"

###Note: vcf and gff need to be edited so that chr and sample names do not include special characters
cat $VCF | \
sed 's/_H.merge.bam//g' | \
sed 's/_M.merge.bam//g' | \
sed 's/_H.merge.dwn.bam//g' | \
sed 's/_M.merge.dwn.bam//g' | \
sed -Ee 's/NC\_(\S*)\.1/chr\1/g' | \
sed -Ee 's/NW\_(\S*)\.1/scaff\1/g' | \
sed -Ee 's/NC\_(\S*)\.2/chr\1/g' > $DIR/pman.noLow.dwnSamp.fixednames.vcf

zcat $REFDIR/GCF_003704035.1_HU_Pman_2.1.3_genomic.gff.gz | \
sed -Ee 's/NC\_(\S*)\.1/chr\1/g' | \
sed -Ee 's/NW\_(\S*)\.1/scaff\1/g' | \
sed -Ee 's/NC\_(\S*)\.2/chr\1/g' > $DIR/fixedchrom.gff

vcftools --vcf $DIR/pman.noLow.dwnSamp.fixednames.vcf --012 --out $DIR/pman.noLow.dwnSamp.fixednames

cd $DIR

python3 LD-annot0.4.py pman.noLow.dwnSamp.fixednames.vcf fixedchrom.gff pman.noLow.dwnSamp.fixednames.012.pos gene 0.9 pman.noLow.dwnSamp_0.9