#!/bin/bash
#SBATCH --job-name=dxy-array
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --mem=2G
#SBATCH --account=bradburd0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH  --mail-user=nicolead@umich.edu
#SBATCH --array=1-21                  # Array jobs, number of chrs

# Modified from David Marques' run10kdxy.pbs  
# Submit this job like this: run10kdxy.sh pop1 pop2 win10k.bed

# Population names
a=$1
b=$2

# windowed bed file
BED=$3

# Define programs and constants
ANGSDIR=/home/akimmitt/angsd.94
WINDIR=/home/akimmitt/.cargo/bin
ANC=flinflon/$REF
FAI=reference.fai #

# Get index number for chromosomes
i=$SLURM_ARRAY_TASK_ID
d=$(cat GCF_003704035.1_HU_Pman_2.1.3_genomic.fna.fai_autoNames.txt) #chr names, I think
c=$(echo $d | cut -f $i -d " ")

# Choose SAF file based on autosome or sex chromosome
cat "${SLURM_SUBMIT_DIR}/../../regions/windows.chr${c}.10k.bed" | while read -r chr start end; do
    region="${chr}:${start}-${end}"
    
    echo "$region" > "${SLURM_SUBMIT_DIR}/${i}.rf"

    angsd -b "${SLURM_SUBMIT_DIR}/../bamlists/${a}.bamlist" \
          -out "${SLURM_SUBMIT_DIR}/tmp.${a}.${i}" \
          -doSaf 1 -GL 2 -rf "${SLURM_SUBMIT_DIR}/${i}.rf" \
          -minQ 30 -minMapQ 30 \
          -fai $FAI \
          -anc $REF

    angsd -b "${SLURM_SUBMIT_DIR}/../bamlists/${b}.bamlist" \
          -out "${SLURM_SUBMIT_DIR}/tmp.${b}.${i}" \
          -doSaf 1 -GL 2 -rf "${SLURM_SUBMIT_DIR}/${i}.rf" \
          -minQ 30 -minMapQ 30 \
          -fai $FAI \
          -anc $ANC

    echo -n "$chr $start $end "
    
    $WINDIR/winsfs "${SLURM_SUBMIT_DIR}/tmp.${a}.${i}.saf.idx" \
            "${SLURM_SUBMIT_DIR}/tmp.${b}.${i}.saf.idx"
    
done | awk '{if(NF>4){print $NF}}' > "${SLURM_SUBMIT_DIR}/${a}.${b}.chr${c}.10kwin.sfs"

#rm $SLURM_SUBMIT_DIR"/"$i".rf" $SLURM_SUBMIT_DIR"/"tmp.*.$i.*