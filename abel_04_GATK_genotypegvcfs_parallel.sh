#!/bin/bash

# Job name:
#SBATCH --job-name=genotypegvcfs

# Project:
# dont change
#SBATCH --account=uio       

# Wall clock limit:
# not too high
#SBATCH --time=20:00:00  

# Max memory usage per core (MB) (NB usually 16cores*4G = 64G):
#SBATCH --mem-per-cpu=3G

# Number of cores:
#SBATCH --cpus-per-task=16

## Set up job environment      
## never change, it is sourcing an script
source /cluster/bin/jobsetup

## loading modules
module load gatk 

## Copy input files to the work directory:
cp $SUBMITDIR/*.raw.g.vcf* $SCRATCH
cp ~/aa/JIC_reference/*.* $SCRATCH

## Mark outfiles for automatic copying to $SUBMITDIR:
## copy all output files
chkfile "*raw.vcf*"

## Do some work:
cd $SCRATCH

#first prepare list of all samples
ls *.raw.g.vcf.gz | sed 's/^/-V /' | sed 's/$/ /' > join.vcf.samplelist.txt
SAMPLELIST=$(<join.vcf.samplelist.txt)
echo "Joint genotyping samples $SAMPLELIST"

#then run Genotype GVCFs
# NB check that (i) -Xmx has cpus-per-task*mem-per-cpu (e.g. 16*3=48) and (ii) -nt = cpus-per-task
java -Xmx48g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R alygenomes.fasta \
$SAMPLELIST \
-nt 16 \
--includeNonVariantSites \
-o joint.raw.vcf.gz

