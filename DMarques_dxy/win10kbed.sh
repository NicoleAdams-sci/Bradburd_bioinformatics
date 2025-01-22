#!/bin/bash

# Run this script like this: bash win10kbed.sh GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.fai

module load Bioinformatics bedtools2

FAI=$1

# create bed file from reference index
cat $FAI | cut -f1,2,3 > $FAI.bed

# create 10kb windows per chr
bedtools makewindows -g $FAI.bed -w 10000 > $FAI.10kb.bed

# make list of chr names
less $FAI.bed | cut -f1 |head -24 | tr '\n' ' ' > $FAI.autoNames.txt
