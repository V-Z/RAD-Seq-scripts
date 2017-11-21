#!/bin/bash

# Job name:
#SBATCH --job-name=selectvariants

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
cp $SUBMITDIR/joint.raw.vcf.gz $SCRATCH
cp ~/aa/JIC_reference/*.* $SCRATCH
cp ~/aa/JIC_reference/Blacklist.Excess.Het.List.Genes.With.5orMore.2PopsMin.AllHet.Sites.GATK.intervals $SCRATCH
cp ~/aa/JIC_reference/Blacklist_allscafs_ExcSites1.6x_maxRD.GATK.intervals $SCRATCH

## Mark outfiles for automatic copying to $SUBMITDIR:
## copy all output files
chkfile "*raw.SNP*"

## Do some work:
cd $SCRATCH

# hardfilter SNPs from this joint file using recommended parameters for SNPS + exclude SNPs from set of putatively paralogou loci and with excess depth (from Harvard)
java -Xmx30g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R alygenomes.fasta \
   -V joint.raw.vcf.gz \
   -selectType SNP \
   -o joint.raw.SNP.vcf.gz \
   -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 \
   --excludeIntervals Blacklist.Excess.Het.List.Genes.With.5orMore.2PopsMin.AllHet.Sites.GATK.intervals \
   --excludeIntervals Blacklist_allscafs_ExcSites1.6x_maxRD.GATK.intervals \
   --excludeNonVariants





