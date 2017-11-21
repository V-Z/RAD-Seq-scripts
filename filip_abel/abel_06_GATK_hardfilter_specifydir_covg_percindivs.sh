#!/bin/bash

REFSEQ="/usit/abel/u1/filipko/aa/JIC_reference/alygenomes.fasta"
EXCLUDEPARALOGS="/usit/abel/u1/filipko/aa/JIC_reference/excess_heterozygotes_AAgenomes_harvard.intervals"

# I have joint.raw.vcf.gz file in a specified folder, e.g., named jointgenotyping_tet_run12
# I have to specify the name of this dir when launching this script, then min coverage to account for a SNP and min. proportion of indivs (e.g. 0.9 for 90%) that contain that SNP in the whole dataset defined by given coverage

# usage: bash abel_05_GATK_hardfilter_specifydir_covg_percindivs.sh jointgenotyping_tet_run12 4 0.9

AnalysisDir=$1
MinCoverage=$2
MinPercIndivs=$3

echo "$AnalysisDir"
echo "MinCoverage: $MinCoverage"
echo "MinPercIndivs: $MinPercIndivs"
cd $AnalysisDir

: <<'EOF'

# hardfilter SNPs from this joint file using recommended parameters for SNPS + exclude SNPs from set of putatively paralogou loci (from Harvard)
echo "Extracting non-paralogous SNPs from chr 1-8 joint vcf..."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.raw.vcf.gz \
   -selectType SNP \
   -o joint.raw.SNP.vcf.gz \
   -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 \
   --excludeIntervals $EXCLUDEPARALOGS \
   --excludeNonVariants

EOF

#create a filtering criterion
echo "Hard Filtering SNPs from joint vcf..."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T VariantFiltration \
   -R $REFSEQ \
   -V joint.raw.SNP.vcf.gz \
   -o joint.hardfilt.SNP.vcf.gz \
   --filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5" \
   --filterName "AA_filter1"
# now possible to view with IGV - the no-PASSSing variants will be marked by lighter colours

# extract SNPs (both bi- and multiallelic) passing the AA_filter1 criterion
echo "Extracting only passing SNPs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.vcf.gz \
   -o joint.hardfilt.SNP.PASS.vcf.gz \
   --excludeFiltered

# extract BIALLELIC SNPs
echo "Extracting only BIALLELIC SNPs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.vcf.gz \
   -o joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   --excludeFiltered \
   --restrictAllelesTo BIALLELIC

# this tool works until 3.4 ... does not work on multiallelic SNPs produced by 3.5, just biallelic
# emit only sites in which x% of samples have sufficient coverage (e.g. 50% at least 8x coverage)
echo "Creating intervals for loci covered min $MinCoverage-times in defined $MinPercIndivs of indivs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.3.0/GenomeAnalysisTK.jar -T CoveredByNSamplesSites \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   -out dp$MinCoverage.perc$MinPercIndivs.BIALL.intervals \
   --minCoverage $MinCoverage \
   --percentageOfSamples $MinPercIndivs

# select variants based on this interval list (NB variants with < defined coverage will be still present in VCF)
echo "Selecting variants based on presence in $MinPercIndivs of indivs"
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   -o joint.hardfilt.SNP.PASS.BIALL.dp$MinCoverage.perc$MinPercIndivs.vcf.gz \
   -L dp$MinCoverage.perc$MinPercIndivs.BIALL.intervals

#EOF

echo Finished!!!







