#!/bin/bash
#################### Pop structure analysis - all samples ####################
# This bash script filters (globally) a genomic VCF then performs population structure analyses

# Usage: bash code/pop_structure_all_samples.sh input.vcf.gz metadata <imiss> <lmiss> <hwp> <k> <threads>

# Load modules
module load Bioinformatics bcftools vcftools htslib plink Rtidyverse

# Read in named and positional inputs
VCF=$1
META=$2		# Should have column headers
shift 2  # Remove first two positional args

# Set defaults
IMISS=0.5	# individual missingness, default 0.5 (50%)
LMISS=0.5	# locus missingness, default 0.5 (50%)
HWP=0.000001	# pvalue for HWE cutoff filtering, default 1e-6
K=5			# Maximum value of clusters (K) want to test with admixture, default 5 clusters
THREADS=4

# Read in named inputs
for arg in "$@"; do
    case $arg in
        imiss=*) IMISS="${arg#*=}" ;;
        lmiss=*) LMISS="${arg#*=}" ;;
        hwp=*) HWP="${arg#*=}" ;;
        k=*) K="${arg#*=}" ;;
        threads=*) THREADS="${arg#*=}" ;;
        *) echo "Unknown parameter: $arg" ;;
    esac
done


# Define outputs
mkdir -p output/pop_struct_all
OUTDIR="output/pop_struct_all"

OUTNAM="pop_all"

META_NOHEAD=${META%.*}\_nohead.txt
tail -n +2 $META > $META_NOHEAD


# Write running output to screen and log file
exec > >(tee $OUTDIR/pop_structure_all_samples.log) 2>&1

# Convert individual missingness threshold from decimal to percent to use in file names
IMISS2=$(echo "$IMISS * 100" | bc | cut -f1 -d".")
LMISS2=$(echo "$LMISS * 100" | bc | cut -f1 -d".")

echo "Using $IMISS2% individual missingness and $LMISS2% locus missingness as cutoffs"


########## Filter VCF for population structure analyses ##########
# Basic quality filters
## minimum genotype quality -- can adjust; minimum genotype depth -- can adjust, eep only biallelic alleles
vcftools --gzvcf $VCF \
		--minGQ 15 \
		--minDP 10 \
        --remove-indels \
        --min-alleles 2 \
        --max-alleles 2 \
        --recode \
        --recode-INFO-all \
        --stdout | bgzip -c > $OUTDIR/$OUTNAM.vcf.gz

cd $OUTDIR

# # Calculate individual sample missingness
vcftools --gzvcf $OUTNAM.vcf.gz --missing-indv --out $OUTNAM

# Create a list of individuals that pass the missingness threshold
cat $OUTNAM.imiss | awk -v imiss=$IMISS '$5 < imiss' | awk '{print $1}' > $OUTNAM.$IMISS2.txt

# Filter VCF using the list of individuals that pass the missingness threshold
# And filter based on locus missingness
## minor allele frequency cutoff, locus missingness cutoff, minor allele frequency cutoff
vcftools --gzvcf $OUTNAM.vcf.gz \
  --keep $OUTNAM.$IMISS2.txt \
  --max-missing $LMISS \
  --maf 0.05 \
  --recode \
  --recode-INFO-all \
  --stdout | bgzip -c > $OUTNAM.i$IMISS2.l$LMISS2.maf05.vcf.gz


# Check if first chromosome is numeric
echo "Checking if chromosomes are numeric..."
FIRST_CHROM=$(bcftools view -H $OUTNAM.i$IMISS2.l$LMISS2.maf05.vcf.gz | head -1 | cut -f1)

if [[ $FIRST_CHROM =~ ^[0-9]+$ ]]; then
    echo "Adding prefix to numeric chromosomes..."
    # More direct sed approach
    zcat $OUTNAM.i$IMISS2.l$LMISS2.maf05.vcf.gz | sed 's/^[^#]/scaffold_&/' | bgzip -c > $OUTNAM.i$IMISS2.l$LMISS2.maf05.scaff.vcf.gz
    VCF_TO_USE="$OUTNAM.i$IMISS2.l$LMISS2.maf05.scaff.vcf.gz"
else
    VCF_TO_USE="$OUTNAM.i$IMISS2.l$LMISS2.maf05.vcf.gz"
fi


# Filter for independent sites
## Identify list of sites in/out of LD using Plink
mkdir -p plink_results
cd plink_results

plink --vcf ../$VCF_TO_USE --allow-extra-chr --make-bed --set-all-var-ids @_# --threads $THREADS --out $OUTNAM.i$IMISS2.l$LMISS2.maf05

plink --bfile $OUTNAM.i$IMISS2.l$LMISS2.maf05 --allow-extra-chr --indep-pairwise 50 5 0.5 --out $OUTNAM.i$IMISS2.l$LMISS2.maf05_ld


# echo "--- reformat file of LD site to rm for $VCF ---"
sed 's/scaffold_\([^_]*\)_\([^_]*\)_.*/scaffold_\1\t\2/' $OUTNAM.i$IMISS2.l$LMISS2.maf05_ld.prune.out > $OUTNAM.i$IMISS2.l$LMISS2.maf05_ld.prune.reformat.out


echo "--- rm LD sites from $VCF ---" and filter on HWE
cd ..

vcftools --gzvcf $VCF_TO_USE \
  --exclude-positions plink_results/$OUTNAM.i$IMISS2.l$LMISS2.maf05_ld.prune.reformat.out \
  --hwe $HWP \
  --recode \
	--recode-INFO-all \
	--stdout | bgzip -c > $OUTNAM.filtered.vcf.gz


# Index vcf
tabix $OUTNAM.filtered.vcf.gz

# # Count the number of samples and variants in the filtered VCF
echo "--- No. of samples in $OUTNAM.filtered.vcf.gz ---"
bcftools query -l $OUTNAM.filtered.vcf.gz | wc -l

echo "--- No. of sites in $OUTNAM.filtered.vcf.gz ---"
bcftools view -H $OUTNAM.filtered.vcf.gz --threads $THREADS | wc -l



########## PCA using Plink ##########
echo "--- Running PCA ---"
plink --vcf $OUTNAM.filtered.vcf.gz --allow-extra-chr --make-bed --set-all-var-ids @_# --export ped --threads $THREADS --pca --out plink_results/$OUTNAM.filtered

# Plot PCA using R
mkdir -p plink_pca_results
Rscript ../../code/plot_pca.R plink_results/$OUTNAM.filtered.eigenval plink_results/$OUTNAM.filtered.eigenvec ../../$META pca


########## PCA using VCF2PCACluster (Bradburd module) ##########
module use /nfs/turbo/lsa-bradburd/shared/Lmod/
module load VCF2PCACluster
mkdir -p vcf2pca_results
VCF2PCACluster -InVCF $OUTNAM.filtered.vcf.gz -OutPut vcf2pca_results/pca -InSampleGroup ../../$META_NOHEAD


########## Structure analysis with Admixture (Bradburd module) ##########
module use /nfs/turbo/lsa-bradburd/shared/Lmod/
module load admixture

echo "--- Run K=2 to $K with cross-validation ---"
mkdir -p admixture_results

# copy plink files to admix directory to run
cp plink_results/$OUTNAM.filtered.{bed,bim,fam,map} admixture_results/

cd admixture_results

# fix chr format bc admix only accepts integer chr names - replace first column with zeros
for file in $OUTNAM.filtered.bim $OUTNAM.filtered.map; do
    if [ -f "$file" ]; then
        echo "Fixing chromosome format in $file"
        awk '{$1="0"; print $0}' "$file" > "${file}.tmp"
        mv "$file" "${file}_og"  # backup original
        mv "${file}.tmp" "$file"
    fi
done

# run admixture for every K with bootstrapping (-B#, can add and change bootstrap value *bootstrap not working atm*)
for i in $(seq 2 $K); do
    admixture --cv $OUTNAM.filtered.bed $i -j$THREADS| tee log${i}.out
done

# Consolidate cross validation 
grep -h CV log*.out > cross_validation.txt

# Plot admixture results using R
Rscript ../../../code/plot_admixture.R . $OUTNAM cross_validation.txt $OUTNAM.filtered.fam


########## Structure analysis using fastStructure (Bradburd module) ##########
# stay in admixture results bc fastStructure also needs the reformatted bim file

# load fastStructure module via structure_threader
module purge
module use /nfs/turbo/lsa-bradburd/shared/Lmod/
module load structure_threader/1.0

# Make indfile input for fastStructure from Plink fam file
cat $OUTNAM.filtered.fam | awk '{print $2}' > $OUTNAM.ind

# Run fastStructure for every K with 3 replicates each (although fs might ignore -R)
structure_threader run -K $K -R 3 -i $OUTNAM.filtered.fam --ind $OUTNAM.ind -o ../fastStructure_results/ -t $THREADS -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure


########## Admixture events analysis with TreeMix ##########
#mkdir -p $OUTDIR/treemix_results
#cd $OUTDIR/treemix_results
