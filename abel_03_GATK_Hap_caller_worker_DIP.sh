#!/bin/bash

# Job name:
#SBATCH --job-name=HC

# Project:
# dont change
#SBATCH --account=uio       

# Wall clock limit:
# not too high
#SBATCH --time=10:00:00  

# Max memory usage per core (MB) (NB usually 16cores*4G = 64G):
#SBATCH --mem-per-cpu=30G

## Set up job environment      
## never change, it is sourcing an script
source /cluster/bin/jobsetup

## loading modules
module load gatk 

## Copy input files to the work directory:
# ONE BAM file named AA016_a.realign.bam (and AA016_a.realign.bai index) are in dir named "AA016_a.1" (the.1 should vary from 1 ... X across the dirs)
# NB. to do so use command: n=1; for dir in AA*; do mv "$dir" "$(printf "$dir.$n")";   n=$((n+1)); done
BAMDIR=~/aa/GATK/tetraploids_june_2016/diploids/AA?????.$TASK_ID
REF=~/aa/JIC_reference

cp $BAMDIR/*realign.* $SCRATCH
cp $REF/*.* $SCRATCH

# copy output files
chkfile "*vcf*"

# do some work
# Run HaplotypeCaller on each DIPLOID sample to produce genomicVCF
cd $SCRATCH

bamfile=(*.realign.bam)
BAM_file_base=${bamfile%.realign.bam}

java -Xmx30g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller \
  -R alygenomes.fasta \
  -I $BAM_file_base".realign.bam" \
  -ERC GVCF \
  -variant_index_type LINEAR \
  -variant_index_parameter 128000 \
  -ploidy 2 \
  -o $BAM_file_base".raw.g.vcf.gz"



