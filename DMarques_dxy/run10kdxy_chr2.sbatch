#!/bin/bash
#SBATCH --job-name run10kdxy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=1G
#SBATCH --time=24:00:00
##SBATCH --account=test
#SBATCH --partition=standard
#SBATCH --output run10kdxy.%j.out
##SBATCH --error run10kdxy.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nicoleadams.sci@gmail.com
#SBATCH --array=1-25

# set reference and bed files directory
REF="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna"
FAI="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.fai"
BEDDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/reference/10kb_beds"
BAMDIR="/nfs/turbo/lsa-bradburd/NicoleAdams/pman/bwaMap/mergedBams"
OUTDIR="nfs/turbo/lsa-bradburd/NicoleAdams/pman/genoL.downSamp/dxy/sfs_chr"
# Get index number for chromosomes
i=${SLURM_ARRAY_TASK_ID}
d=$(cut -f1 GCF_003704035.1_HU_Pman_2.1.3_genomic.10kb.bed | sort | uniq | grep "NC") # only the chr, not the scaffolds
c=$(echo $d | cut -f $i -d " ")
# Population names
a=$1
b=$2
aNam=$3
bNam=$4
# deinfe bed and output names
input_bed="$BEDDIR/${c}.10kb.bed"
output_sfs="$OUTDIR/${aNam}.${bNam}.chr${c}.10kwin.sfs"
# Angsd filters
FILT="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 25 -C 50 -baq 1"
# Choose SAF file based on autosome or sex chromosome
# Process each line in the input BED file
while read -r chrom start end; do
  region_file="$OUTDIR/${i}.rf"
  echo "${chrom}:${start}-${end}" > "$region_file"
  
    # Output files for angsd
  output_a="${OUTDIR}/tmp.${a}.${i}"
  output_b="${OUTDIR}/tmp.${b}.${i}"

  # Run ANGSD for the first set of data, if the saf file does not exist
  if [ ! -f "${output_a}.saf.idx" ]; then
    angsd -P "${SLURM_NTASKS_PER_NODE}" -i "${a}" -out "${output_a}" \
          -doSaf 1 -GL 2 -rf "${region_file}" "${FILT}" \
          -fai "${FAI}" -anc "${REF}"
  fi

  # Run ANGSD for the second set of data, if the saf file does not exist
  if [ ! -f "${output_b}.saf.idx" ]; then
    angsd -P "${SLURM_NTASKS_PER_NODE}" -i "${b}" -out "${output_b}" \
          -doSaf 1 -GL 2 -rf "${region_file}" "${FILT}" \
          -fai "${FAI}" -anc "${REF}"
  fi

  # Calculate the site frequency spectrum, if both saf files exist
  if [ -f "${output_a}.saf.idx" ] && [ -f "${output_b}.saf.idx" ]; then
    result=$(realSFS "${output_a}.saf.idx" "${output_b}.saf.idx")
    # Output the chromosome, start, end positions, and the result
    echo "${chrom} ${start} ${end} ${result}" >> "$output_sfs"
  else
    echo "Error: SAF index files missing for region ${chrom}:${start}-${end}" >&2
  fi

  # Cleanup temporary files
  rm -f "${region_file}" "${output_a}.saf.pos.gz" "${output_b}.saf.pos.gz" "${output_a}.saf.gz" "${output_b}.saf.gz"

done < "$input_bed"
