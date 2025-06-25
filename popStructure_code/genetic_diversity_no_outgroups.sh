#!/bin/bash
#################### Genetic diversity analyses - no outgroups ####################
# This bash script filters a genomic VCF first globally then within genetic groups identified using general population structure analyses (i.e. pop_structure_all_samples.sh) then performs genetic diversity calculations and analyses

# Usage: bash code/genetic_diversity_no_outgroups.sh input.vcf.gz metadata geneticgroups <imiss> <lmiss> <hwp> <threads>

# Load modules
module load Bioinformatics bcftools vcftools htslib plink Rtidyverse

# Read in named and positional inputs
VCF=$1
META=$2		# Metadata file; should have column headers
GENGRPS=$3		# File of genetic groups assignments; should have column headers
shift 3  # Remove first two positional args

# Set defaults
IMISS=0.5	# individual missingness, default 0.5 (50%)
LMISS=0.5	# locus missingness, default 0.5 (50%)
HWP=0.000001	# pvalue for HWE cutoff filtering, default 1e-6
THREADS=4

# Read in named inputs
for arg in "$@"; do
    case $arg in
        imiss=*) IMISS="${arg#*=}" ;;
        lmiss=*) LMISS="${arg#*=}" ;;
        hwp=*) HWP="${arg#*=}" ;;
        threads=*) THREADS="${arg#*=}" ;;
        *) echo "Unknown parameter: $arg" ;;
    esac
done


# Define outputs
mkdir -p output/gen_div
OUTDIR="output/gen_div"

OUTNAM="gen_div"

META_NOHEAD=${META%.*}\_nohead.txt
tail -n +2 $META > $META_NOHEAD


# Write running output to screen and log file
exec > >(tee $OUTDIR/genetic_diversity_no_outgroups.log) 2>&1

# Convert individual missingness threshold from decimal to percent to use in file names
IMISS2=$(echo "$IMISS * 100" | bc | cut -f1 -d".")
LMISS2=$(echo "$LMISS * 100" | bc | cut -f1 -d".")

echo "Using $IMISS2% individual missingness and $LMISS2% locus missingness as cutoffs"


########## Filter VCF for population structure analyses ##########
# # Basic quality filters
# ## minimum genotype quality -- can adjust; minimum genotype depth -- can adjust, eep only biallelic alleles
# vcftools --gzvcf $VCF \
# 		--minGQ 15 \
# 		--minDP 10 \
#         --remove-indels \
#         --min-alleles 2 \
#         --max-alleles 2 \
#         --recode \
#         --recode-INFO-all \
#         --stdout | bgzip -c > $OUTDIR/$OUTNAM.vcf.gz

cd $OUTDIR

# # Calculate individual sample missingness
# vcftools --gzvcf $OUTNAM.vcf.gz --missing-indv --out $OUTNAM
# 
# # Create a list of individuals that pass the missingness threshold
# cat $OUTNAM.imiss | awk -v imiss=$IMISS '$5 < imiss' | awk '{print $1}' > $OUTNAM.$IMISS2.txt
# 
# # Filter VCF using the list of individuals that pass the missingness threshold
# # And filter based on locus missingness
# vcftools --gzvcf $OUTNAM.vcf.gz \
#   --keep $OUTNAM.$IMISS2.txt \
#   --max-missing $LMISS \
#   --recode \
#   --recode-INFO-all \
#   --stdout | bgzip -c > $OUTNAM.i$IMISS2.l$LMISS2.vcf.gz


########## Population-specific MAF and HWE filtering ##########
# Get list of unique genetic groups (assuming column 2 has the group assignments)
GENETIC_GROUPS=$(tail -n +2 ../../$GENGRPS | cut -f2 | sort | uniq)

# Create sample lists for each genetic group
for GROUP in $GENETIC_GROUPS; do
    echo "Creating sample list for genetic group: $GROUP"
    tail -n +2 ../../$GENGRPS | awk -v grp="$GROUP" '$2 == grp {print $1}' > ${OUTNAM}_${GROUP}_samples.txt
    echo "Group $GROUP has $(wc -l < ${OUTNAM}_${GROUP}_samples.txt) samples"
done

# Filter each genetic group separately for MAF and HWE
# FILTERED_VCFS=""
# for GROUP in $GENETIC_GROUPS; do
#     echo "Processing genetic group: $GROUP"
#     
#     # Check if sample list is not empty
#     if [ -s ${OUTNAM}_${GROUP}_samples.txt ]; then
#         # Filter this genetic group with MAF and HWE
#         vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.vcf.gz \
#           --keep ${OUTNAM}_${GROUP}_samples.txt \
#           --maf 0.01 \
#           --hwe $HWP \
#           --recode \
#           --recode-INFO-all \
#           --stdout | bgzip -c > ${OUTNAM}_${GROUP}_i${IMISS2}.l${LMISS2}.maf05.hwe.vcf.gz
#         
#         # Add to list of filtered VCFs
#         FILTERED_VCFS="$FILTERED_VCFS ${OUTNAM}_${GROUP}_i${IMISS2}.l${LMISS2}.maf05.hwe.vcf.gz"
#         
#         echo "Completed filtering for group $GROUP"
#     else
#         echo "Warning: No samples found for genetic group $GROUP"
#     fi
# done


# Merge the population-filtered VCFs back together
# echo "Merging population-filtered VCFs..."
# if [ -n "$FILTERED_VCFS" ]; then
#     # Index all VCFs first
#     for vcf in $FILTERED_VCFS; do
#         tabix $vcf
#     done
#     
#     # Merge using bcftools
#     bcftools merge $FILTERED_VCFS -O z -o $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz
#     
#     # Index the merged file
#     tabix $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz
#     
#     echo "Population-specific filtering complete. Final VCF: $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz"
#     
#     # Optional: Clean up intermediate files
#     # rm $FILTERED_VCFS
#     # rm ${OUTNAM}_*_samples.txt
#     
# else
#     echo "Error: No filtered VCFs were created"
#     exit 1
# fi


########## Calculate individual heterozygosity and inbreeding coefficient via VCFtools ##########
# echo "Calculating pairwise Fst between all genetic groups..."
# vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
# 	--het \
# 	--out $OUTNAM
# 	
# echo "VCFtools heterozygosity calculations completed!"


########## Calculate Weir and Cockerham's (1984) Fst via VCFtools ##########
## Fst between genetic groups
echo "Calculating pairwise Fst between all genetic groups..."

# Convert GENETIC_GROUPS to array for easier handling
GROUPS_ARRAY=($GENETIC_GROUPS)
NUM_GROUPS=${#GROUPS_ARRAY[@]}

echo "Found $NUM_GROUPS genetic groups for Fst calculations: ${GROUPS_ARRAY[@]}"

# # Calculate pairwise Fst between all group combinations
# for (( i=0; i<$NUM_GROUPS; i++ )); do
#     for (( j=i+1; j<$NUM_GROUPS; j++ )); do
#         GROUP1=${GROUPS_ARRAY[$i]}
#         GROUP2=${GROUPS_ARRAY[$j]}
#         
#         echo "Calculating Fst between $GROUP1 and $GROUP2..."
#         
#         if [ -s ${OUTNAM}_${GROUP1}_samples.txt ] && [ -s ${OUTNAM}_${GROUP2}_samples.txt ]; then
#             # Run vcftools and capture all output
#             vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
#                      --weir-fst-pop ${OUTNAM}_${GROUP1}_samples.txt \
#                      --weir-fst-pop ${OUTNAM}_${GROUP2}_samples.txt \
#                      --out ${OUTNAM}_fst_${GROUP1}_vs_${GROUP2} \
#                      > ${OUTNAM}_fst_${GROUP1}_vs_${GROUP2}.log 2>&1
#             
#             # Extract Fst value directly from log and save to summary
#             if [ -f ${OUTNAM}_fst_${GROUP1}_vs_${GROUP2}.log ]; then
#                 FST_VALUE=$(grep "weighted Fst estimate" ${OUTNAM}_fst_${GROUP1}_vs_${GROUP2}.log | awk '{print $NF}')
#                 echo -e "${GROUP1}\t${GROUP2}\t${FST_VALUE}" >> ${OUTNAM}_fst_summary.txt
#                 echo "Fst between $GROUP1 and $GROUP2: $FST_VALUE"
#             fi
#         fi
#     done
# done
# 
# # Calculate global Fst across all populations simultaneously
# echo "Calculating global Fst across all populations..."
# FST_POP_FLAGS=""
# for GROUP in $GENETIC_GROUPS; do
#     if [ -s ${OUTNAM}_${GROUP}_samples.txt ]; then
#         FST_POP_FLAGS="$FST_POP_FLAGS --weir-fst-pop ${OUTNAM}_${GROUP}_samples.txt"
#     fi
# done
# 
# if [ -n "$FST_POP_FLAGS" ]; then
#     vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
#              $FST_POP_FLAGS \
#              --out ${OUTNAM}_fst_all_populations \
#              > ${OUTNAM}_fst_all_populations.log 2>&1
#     
#     echo "Global Fst calculation complete: ${OUTNAM}_fst_all_populations.weir.fst"
#     
#     # Extract and display the overall Fst value
#     if [ -f ${OUTNAM}_fst_all_populations.log ]; then
#         GLOBAL_FST=$(grep "weighted Fst estimate" ${OUTNAM}_fst_all_populations.log | awk '{print $NF}')
#         if [ -n "$GLOBAL_FST" ]; then
#             echo "Global weighted Fst across all populations: $GLOBAL_FST"
#             echo -e "All_Populations\tGlobal\t$GLOBAL_FST" >> ${OUTNAM}_fst_summary.txt
#         else
#             echo "Warning: Could not extract global Fst value from log file"
# fi
# 
# echo "All Fst calculations completed!"


########## Calculate nucleotide diversity via pixy *** UNDER CONSTRUCTION *** ##########
# Check for invariant (monomorphic) sites in starting VCF
echo "--- Checking for invariant sites in $VCF ---"
bash ../../code/check.invarSites.sh $VCF

#invariant + variant file for pixy

vcftools --vcf populations.all.vcf --remove lowDP.indv --recode --recode-INFO-all --out populations_no_lowDp
vcftools --vcf populations_no_lowDp.recode.vcf --remove ~/shovel-bugs/cryptic_species.txt --recode --stdout | bgzip -c > initial_filtered_indv.vcf.gz 

#separate the invariants from variants

vcftools --gzvcf initial_filtered_indv.vcf.gz --remove-indels --max-alleles 1 --recode --stdout | bgzip -c > invariant_sites.vcf.gz #this keeps invariants

vcftools --gzvcf invariant_sites.vcf.gz --max-missing 0.75 --minDP 20 --recode --stdout | bgzip -c > test_missing_invariant_sites.vcf.gz #this keeps fewer invariants

# create a filtered VCF containing only variant sites
vcftools --gzvcf initial_filtered_indv.vcf.gz --mac 1 --remove-indels --max-missing 0.75 --maf 0.05 --min-meanDP 20  --recode --stdout | bgzip -c > variant_sites.vcf.gz

# index both vcfs using tabix
tabix test_missing_invariant_sites.vcf.gz
tabix variant_sites.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat --allow-overlaps test_missing_invariant_sites.vcf.gz variant_sites.vcf.gz -O z -o final_filtered_sites.vcf.gz

bcftools view --header-only final_filtered_sites.vcf.gz | grep "##contig" | cut -f3 -d "=" | sed 's/>/''/g' | awk '{print $0, " chr" int((NR-1)/10000+1)}'  >chrlist.txt

bcftools annotate --rename-chrs chrlist.txt final_filtered_sites.vcf.gz -Oz -o final_renamed.vcf.gz

bcftools view -H final_renamed.vcf.gz > rename.body.vcf

zgrep -v "^#" final_renamed.vcf.gz |wc -l
echo {1..10603559} | tr ' ' '\n' > contpos.txt

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1; next} {if (FNR in a) $2=a[FNR]; print $0}' contpos.txt rename.body.vcf > rename2.body.vcf

bcftools sort final_renamed.vcf.gz -Oz -o final_sorted.vcf.gz

bcftools view -h final_sorted.vcf.gz > sort.head.txt

cat sort.head.txt rename2.body.vcf | bgzip > final_merged.vcf.gz

tabix final_merged.vcf.gz

##pixy analysis

pixy --stats fst \
--vcf final_merged.vcf.gz \
--populations new_pop_file.txt \
--window_size 10603559 \
--n_cores 4

pixy --stats pi \
--vcf final_merged.vcf.gz \
--populations new_indv_file.txt \
--window_size 10603559 \
--n_cores 4

