#!/bin/bash
# Project: 
#SBATCH --account=uio
# Wall clock limit:
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=3G
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

# load bwa and samtools
module load bwa
module load samtools

# specify input files
# R1 and R2 files named xxxx are in dir named "AA016_a.1" (the.1 should vary from 1 ... X across the dirs)
# NB. to do so use command: n=1; for dir in AA*; do mv "$dir" "$(printf "$dir.$n")";   n=$((n+1)); done
FASTQDIR=~/aa/GATK/tetraploids_june_2016/trimmed/AA?????.$TASK_ID
REF=~/aa/JIC_reference/alygenomes.fasta

cp -r $FASTQDIR $SCRATCH
cp $REF $SCRATCH
cp $REF.* $SCRATCH

cd $SCRATCH
chkfile "*.bam*" "*.bai"

NAME=$(echo $FASTQDIR | cut -f1 -d.)

## If you want the task ID appended to the file name 
#OUTFILE=$NAME.sam.$TASK_ID
## If not
OUTFILE=$NAME  # if it is a bam file, samtools will automatically add *.bam

# mapping separately paired files and the orphaned reads (reads without a mate)
bwa mem $REF $FASTQDIR/*_trm_run2_R1* $FASTQDIR/*_trm_run2_R2* | samtools view -bS - | samtools sort - $OUTFILE"_paired_run2"
bwa mem $REF $FASTQDIR/*_unpaired_run2_R1*  | samtools view -bS - | samtools sort - $OUTFILE"_unpairedR1_run2"
bwa mem $REF $FASTQDIR/*_unpaired_run2_R2*  | samtools view -bS - | samtools sort - $OUTFILE"_unpairedR2_run2"


