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

# Set defaults (for input parameters)
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


# More defaults (non-inputs)
DP=10			# minimum depth cutoff per locus
MAF=0.01		# minor allele freq cutoff

# Define outputs
mkdir -p output/gen_div
OUTDIR="output/gen_div"
# Create FST subdirectory
mkdir -p $OUTDIR/fst

OUTNAM="gen_div"

META_NOHEAD=${META%.*}\_nohead.txt
tail -n +2 $META > $META_NOHEAD


# Write running output to screen and log file
exec > >(tee $OUTDIR/genetic_diversity_no_outgroups.log) 2>&1

# Convert individual missingness threshold from decimal to percent to use in file names
IMISS2=$(echo "$IMISS * 100" | bc | cut -f1 -d".")
LMISS2=$(echo "$LMISS * 100" | bc | cut -f1 -d".")

echo "Using $IMISS2% individual missingness and $LMISS2% locus missingness as cutoffs"


# ########## Filter VCF for population structure analyses ##########
# # Basic quality filters
# ## minimum genotype quality -- can adjust; minimum genotype depth -- can adjust, eep only biallelic alleles
# vcftools --gzvcf $VCF \
# 		--minGQ 15 \
# 		--minDP $DP \
#         --remove-indels \
#         --min-alleles 2 \
#         --max-alleles 2 \
#         --recode \
#         --recode-INFO-all \
#         --stdout | bgzip -c > $OUTDIR/$OUTNAM.vcf.gz

cd $OUTDIR

# # # Calculate individual sample missingness
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
# 

########## Population-specific MAF and HWE filtering ##########
# Get list of unique genetic groups (assuming column 2 has the group assignments)
GENETIC_GROUPS=$(tail -n +2 ../../$GENGRPS | cut -f2 | sort | uniq)

# # Create sample lists for each genetic group
# for GROUP in $GENETIC_GROUPS; do
#     echo "Creating sample list for genetic group: $GROUP"
#     tail -n +2 ../../$GENGRPS | awk -v grp="$GROUP" '$2 == grp {print $1}' > ${OUTNAM}_${GROUP}_samples.txt
#     echo "Group $GROUP has $(wc -l < ${OUTNAM}_${GROUP}_samples.txt) samples"
# done
# 
# # Filter each genetic group separately for MAF and HWE
# FILTERED_VCFS=""
# for GROUP in $GENETIC_GROUPS; do
#     echo "Processing genetic group: $GROUP"
#     
#     # Check if sample list is not empty
#     if [ -s ${OUTNAM}_${GROUP}_samples.txt ]; then
#         # Filter this genetic group with MAF and HWE
#         vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.vcf.gz \
#           --keep ${OUTNAM}_${GROUP}_samples.txt \
#           --maf $MAF \
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
# 
# 
# # Merge the population-filtered VCFs back together
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
# 

########## Calculate individual heterozygosity and inbreeding coefficient via VCFtools ##########
# echo "Calculating pairwise Fst between all genetic groups..."
# vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
# 	--het \
# 	--out $OUTNAM
# 	
# echo "VCFtools heterozygosity calculations completed!"

# Function to calculate pairwise FST between two populations
calculate_pairwise_fst() {
    local pop1=$1
    local pop2=$2
    local outnam=$3
    local vcf_file=$4
    local sample_suffix=$5
    local output_dir=$6
    
    echo "Calculating Fst between $pop1 and $pop2..."
    
    if [ -s ${outnam}_${pop1}_${sample_suffix}_samples.txt ] && [ -s ${outnam}_${pop2}_${sample_suffix}_samples.txt ]; then
        # Run vcftools and capture all output
        vcftools --gzvcf $vcf_file \
                 --weir-fst-pop ${outnam}_${pop1}_${sample_suffix}_samples.txt \
                 --weir-fst-pop ${outnam}_${pop2}_${sample_suffix}_samples.txt \
                 --out ${output_dir}/${outnam}_fst_${sample_suffix}_${pop1}_vs_${pop2} \
                 > ${output_dir}/${outnam}_fst_${sample_suffix}_${pop1}_vs_${pop2}.log 2>&1
        
        # Extract Fst value and save to summary
        if [ -f ${output_dir}/${outnam}_fst_${sample_suffix}_${pop1}_vs_${pop2}.log ]; then
            FST_VALUE=$(grep "weighted Fst estimate" ${output_dir}/${outnam}_fst_${sample_suffix}_${pop1}_vs_${pop2}.log | awk '{print $NF}')
            echo -e "${pop1}\t${pop2}\t${FST_VALUE}" >> ${output_dir}/${outnam}_fst_${sample_suffix}_summary.txt
            echo "Fst between $pop1 and $pop2: $FST_VALUE"
        fi
    else
        echo "Warning: Sample files not found for $pop1 or $pop2"
    fi
}

########## Calculate Weir and Cockerham's (1984) Fst via VCFtools ##########
## Fst between POPULATIONS
echo "Calculating pairwise Fst between populations from metadata..."

# Get list of unique populations from metadata (assuming column with population info)
# Adjust the column number based on your metadata structure
POPULATIONS=$(tail -n +2 ../../$META | cut -f2 | sort | uniq)  # Change f3 to correct column

# Create sample lists for each population
for POP in $POPULATIONS; do
    echo "Creating sample list for population: $POP"
    tail -n +2 ../../$META | awk -v pop="$POP" '$2 == pop {print $1}' > ${OUTNAM}_${POP}_pop_samples.txt  # Change $2 to correct column
    echo "Population $POP has $(wc -l < ${OUTNAM}_${POP}_pop_samples.txt) samples"
done

# Convert POPULATIONS to array for easier handling
POPS_ARRAY=($POPULATIONS)
NUM_POPS=${#POPS_ARRAY[@]}

echo "Found $NUM_POPS populations for Fst calculations: ${POPS_ARRAY[@]}"

# Initialize summary file
echo -e "Population1\tPopulation2\tFst" > fst/${OUTNAM}_fst_pop_summary.txt

# Calculate pairwise Fst between all population combinations in parallel
echo "Starting parallel FST calculations for populations..."
for (( i=0; i<$NUM_POPS; i++ )); do
    for (( j=i+1; j<$NUM_POPS; j++ )); do
        POP1=${POPS_ARRAY[$i]}
        POP2=${POPS_ARRAY[$j]}
        
        if [ ! -f fst/${OUTNAM}_fst_pop_${POP1}_vs_${POP2}.weir.fst ]; then
            # Run in background with job control
            (
                calculate_pairwise_fst "$POP1" "$POP2" "$OUTNAM" "$OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz" "pop" "fst"
            ) &
            
            # Limit number of parallel jobs
            while [ $(jobs -r | wc -l) -ge $THREADS ]; do
                sleep 1
            done
        else
            echo "Fst file for $POP1 vs $POP2 already exists, skipping calculation..."
            # Extract and add to summary if not already there
            if [ -f fst/${OUTNAM}_fst_pop_${POP1}_vs_${POP2}.log ]; then
                FST_VALUE=$(grep "weighted Fst estimate" fst/${OUTNAM}_fst_pop_${POP1}_vs_${POP2}.log | awk '{print $NF}')
                # Check if this pair is already in the summary file
                if ! grep -q "${POP1}.*${POP2}" fst/${OUTNAM}_fst_pop_summary.txt 2>/dev/null; then
                    echo -e "${POP1}\t${POP2}\t${FST_VALUE}" >> fst/${OUTNAM}_fst_pop_summary.txt
                fi
            fi
        fi
    done
done

# Wait for all background jobs to complete
wait

echo "All population pairwise FST calculations completed!"

# Calculate global Fst across all populations simultaneously
echo "Calculating global Fst across all populations (metadata-defined)..."
FST_POP_FLAGS=""
for POP in $POPULATIONS; do
    if [ -s ${OUTNAM}_${POP}_pop_samples.txt ]; then
        FST_POP_FLAGS="$FST_POP_FLAGS --weir-fst-pop ${OUTNAM}_${POP}_pop_samples.txt"
    fi
done

if [ -n "$FST_POP_FLAGS" ]; then
    vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
             $FST_POP_FLAGS \
             --out fst/${OUTNAM}_fst_all_metadata_populations \
             > fst/${OUTNAM}_fst_all_metadata_populations.log 2>&1
    
    echo "Global Fst calculation complete: fst/${OUTNAM}_fst_all_metadata_populations.weir.fst"
    
    # Extract and display the overall Fst value
    if [ -f fst/${OUTNAM}_fst_all_metadata_populations.log ]; then
        GLOBAL_FST=$(grep "weighted Fst estimate" fst/${OUTNAM}_fst_all_metadata_populations.log | awk '{print $NF}')
        if [ -n "$GLOBAL_FST" ]; then
            echo "Global weighted Fst across all metadata populations: $GLOBAL_FST"
            echo -e "All_Metadata_Populations\tGlobal\t$GLOBAL_FST" >> fst/${OUTNAM}_fst_pop_summary.txt
        else
            echo "Warning: Could not extract global Fst value from log file"
        fi
    fi
fi

echo "All population Fst calculations completed!"

## Fst between GENETIC GROUPS
echo "Calculating pairwise Fst between all genetic groups..."

# Create sample lists for each genetic group
for GROUP in $GENETIC_GROUPS; do
    echo "Creating sample list for genetic group: $GROUP"
    tail -n +2 ../../$GENGRPS | awk -v grp="$GROUP" '$2 == grp {print $1}' > ${OUTNAM}_${GROUP}_samples.txt
    echo "Group $GROUP has $(wc -l < ${OUTNAM}_${GROUP}_samples.txt) samples"
done

# Convert GENETIC_GROUPS to array for easier handling
GROUPS_ARRAY=($GENETIC_GROUPS)
NUM_GROUPS=${#GROUPS_ARRAY[@]}

echo "Found $NUM_GROUPS genetic groups for Fst calculations: ${GROUPS_ARRAY[@]}"

# Initialize summary file
echo -e "Group1\tGroup2\tFst" > fst/${OUTNAM}_fst_genetic_group_summary.txt

# Calculate pairwise Fst between all group combinations in parallel
echo "Starting parallel FST calculations for genetic groups..."
for (( i=0; i<$NUM_GROUPS; i++ )); do
    for (( j=i+1; j<$NUM_GROUPS; j++ )); do
        GROUP1=${GROUPS_ARRAY[$i]}
        GROUP2=${GROUPS_ARRAY[$j]}
        
        # Run in background with job control
        (
            calculate_pairwise_fst "$GROUP1" "$GROUP2" "$OUTNAM" "$OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz" "group" "fst"
        ) &
        
        # Limit number of parallel jobs
        while [ $(jobs -r | wc -l) -ge $THREADS ]; do
            sleep 1
        done
    done
done

# Wait for all background jobs to complete
wait

echo "All genetic group pairwise FST calculations completed!"

# Calculate global Fst across all groups simultaneously
echo "Calculating global Fst across all genetic groups..."
FST_GROUP_FLAGS=""
for GROUP in $GENETIC_GROUPS; do
    if [ -s ${OUTNAM}_${GROUP}_samples.txt ]; then
        FST_GROUP_FLAGS="$FST_GROUP_FLAGS --weir-fst-pop ${OUTNAM}_${GROUP}_samples.txt"
    fi
done

if [ -n "$FST_GROUP_FLAGS" ]; then
    vcftools --gzvcf $OUTNAM.i$IMISS2.l$LMISS2.maf05.hwe.merged.vcf.gz \
             $FST_GROUP_FLAGS \
             --out fst/${OUTNAM}_fst_all_genetic_groups \
             > fst/${OUTNAM}_fst_all_genetic_groups.log 2>&1
    
    echo "Global Fst calculation complete: fst/${OUTNAM}_fst_all_genetic_groups.weir.fst"
    
    # Extract and display the overall Fst value
    if [ -f fst/${OUTNAM}_fst_all_genetic_groups.log ]; then
        GLOBAL_FST=$(grep "weighted Fst estimate" fst/${OUTNAM}_fst_all_genetic_groups.log | awk '{print $NF}')
        if [ -n "$GLOBAL_FST" ]; then
            echo "Global weighted Fst across all genetic groups: $GLOBAL_FST"
            echo -e "All_Genetic_Groups\tGlobal\t$GLOBAL_FST" >> fst/${OUTNAM}_fst_genetic_group_summary.txt
        else
            echo "Warning: Could not extract global Fst value from log file"
        fi
    fi
fi

echo "All Fst calculations completed!"

########## Calculate nucleotide diversity via pixy (Bradburd module) ##########
# # Check for invariant (monomorphic) sites in starting VCF
# echo "--- Checking for invariant sites in $VCF ---"
# bash ../../code/check.invarSites.sh $VCF
# 
# echo "Creating combined VCF for pixy analysis..."
# 
# # Step 1: Filter individuals with high missingness and create base filtered VCF
# echo "Filtering individuals with high missingness..."
# vcftools --gzvcf ../../$VCF \
#     --keep $OUTNAM.$IMISS2.txt \
#     --recode --recode-INFO-all --stdout 2> initial_filter_4pixy.log | bgzip -c > initial_filtered_4pixy.vcf.gz
# 
# # Step 2: Create invariant sites VCF in one pipeline
# echo "Processing invariant sites..."
# vcftools --gzvcf initial_filtered_4pixy.vcf.gz \
#     --remove-indels \
#     --max-alleles 1 \
#     --max-missing $LMISS \
#     --minDP $DP \
#     --recode --stdout 2> invariant_sites.log | bgzip -c > invariant_sites_4pixy_filtered.vcf.gz
# 
# # Step 3: Create variant sites VCF in one pipeline
# echo "Processing variant sites..."
# vcftools --gzvcf initial_filtered_4pixy.vcf.gz \
#     --mac 1 \
#     --remove-indels \
#     --max-missing $LMISS \
#     --minDP $DP \
#     --maf $MAF \
#     --recode --stdout 2> variant_sites.log | bgzip -c > variant_sites_4pixy.vcf.gz
# 
# # Step 4: Index and combine VCFs
# echo "Indexing and combining VCFs..."
# tabix invariant_sites_4pixy_filtered.vcf.gz
# tabix variant_sites_4pixy.vcf.gz
# 
# bcftools concat --allow-overlaps invariant_sites_4pixy_filtered.vcf.gz variant_sites_4pixy.vcf.gz -O z -o $OUTNAM\_final_4pixy.vcf.gz
# 
# # Step 5: Rename chromosomes and prepare final VCF
# echo "Renaming chromosomes and preparing final VCF..."
# # Create chromosome renaming file
# bcftools view --header-only $OUTNAM\_final_4pixy.vcf.gz | grep "##contig" | cut -f3 -d "=" | sed 's/>/''/g' | awk '{print $0, " chr" int((NR-1)/10000+1)}' > chrlist.txt
# 
# # Rename chromosomes (this creates a properly renamed VCF with correct header)
# bcftools annotate --rename-chrs chrlist.txt $OUTNAM\_final_4pixy.vcf.gz -Oz -o $OUTNAM\_final_4pixy_renamed.vcf.gz
# 
# # Count variants for position replacement
# VARIANT_COUNT=$(bcftools view -H $OUTNAM\_final_4pixy_renamed.vcf.gz | wc -l)
# echo "Total variants: $VARIANT_COUNT"
# 
# # Create position list
# echo $(seq 1 $VARIANT_COUNT) | tr ' ' '\n' > contpos.txt
# 
# # Extract body and replace positions
# echo "Replacing positions in VCF body..."
# bcftools view -H $OUTNAM\_final_4pixy_renamed.vcf.gz | \
# awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1; next} {if (FNR in a) $2=a[FNR]; print $0}' contpos.txt - | \
# bgzip > processed_body.vcf.gz
# 
# # Sort the renamed VCF first to ensure canonical order
# echo "Sorting renamed VCF..."
# bcftools sort $OUTNAM\_final_4pixy_renamed.vcf.gz -Oz -o final_sorted.vcf.gz
# 
# # Extract header from the sorted, renamed VCF (this ensures header matches content)
# bcftools view -h final_sorted.vcf.gz > final_header.txt
# 
# # Combine correct header with processed body
# cat final_header.txt <(zcat processed_body.vcf.gz) | bgzip > final_merged.vcf.gz
# tabix final_merged.vcf.gz
# 
# # Step 6: Prepare pixy input files
# echo "Preparing pixy input files..."
# # Create individuals file directly from VCF
# bcftools query -l final_merged.vcf.gz | awk 'BEGIN{OFS=FS="\t"} {print $1, $1}' > $OUTNAM\_individuals_4pixy.txt
# 
# # Create populations file by joining with metadata
# awk 'NR==FNR{ids[$1]=1; next} $1 in ids' $OUTNAM\_individuals_4pixy.txt $META_NOHEAD > $OUTNAM\_pops_4pixy.txt
# 
# # Step 7: Run pixy analyses
# echo "Running pixy analyses..."
# module purge; module use /nfs/turbo/lsa-bradburd/shared/Lmod/; module load pixy
# 
# echo "Calculating Fst with pixy..."
# pixy --stats fst \
#     --vcf final_merged.vcf.gz \
#     --populations $OUTNAM\_pops_4pixy.txt \
#     --window_size $VARIANT_COUNT \
#     --n_cores $THREADS
# 
# echo "Calculating nucleotide diversity (pi) with pixy..."
# pixy --stats pi \
#     --vcf final_merged.vcf.gz \
#     --populations $OUTNAM\_individuals_4pixy.txt \
#     --window_size $VARIANT_COUNT \
#     --n_cores $THREADS
# 
# # Optional: Clean up intermediate files (uncomment to enable)
# # echo "Cleaning up intermediate files..."
# # rm initial_filtered_4pixy.vcf.gz invariant_sites_4pixy_filtered.vcf.gz variant_sites_4pixy.vcf.gz
# # rm $OUTNAM\_final_4pixy.vcf.gz $OUTNAM\_final_4pixy_renamed.vcf.gz final_sorted.vcf.gz
# # rm processed_body.vcf.gz final_header.txt chrlist.txt contpos.txt
# 
# echo "Pixy analysis completed!"