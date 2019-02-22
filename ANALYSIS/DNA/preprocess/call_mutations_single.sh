#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 40G
#SBATCH -t 0-40:00

cd ../../../PREPROCESS/DNA/
mkdir -p mutect-results

module load centos6/0.0.1-fasrc01
module load java/1.7.0_60-fasrc01

# run MuTect if output file does not already exist 

if [ ! \( -f mutect-results/$(basename $TUMOR_BAM .bam)\.vcf \) ] ; then
  java -Xmx2g -jar mutect-src/mutect/target/mutect-1.1.7.jar \
    --analysis_type MuTect \
    --reference_sequence annotation/GATK_bundle_b37/human_g1k_v37.fasta \
    --dbsnp annotation/GATK_bundle_b37/dbsnp_138.b37.vcf \
    --input_file:normal ../../DATA/DNA/$NORMAL_BAM \
    --input_file:tumor ../../DATA/DNA/$TUMOR_BAM \
    --out mutect-results/$(basename $TUMOR_BAM .bam)\_call_stats.txt \
    --coverage_file mutect-results/$(basename $TUMOR_BAM .bam)\_coverage.wig.txt \
    -vcf mutect-results/$(basename $TUMOR_BAM .bam)\.vcf 
fi

# can't use the COSMIC file, or get a java error....
#--cosmic annotation/COSMIC/b37_cosmic_v54_120711.vcf \

# for testing short regions use, e.g.
# --intervals 1:10000000-10100000 \

# on output, extract PASS mutations in separate vcf file
if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_muts.vcf ]; then
  grep -v REJECT mutect-results/$(basename $TUMOR_BAM .bam)\.vcf  > mutect-results/$(basename $TUMOR_BAM .bam)\_muts.vcf 
fi

##################### run vcf2maf
# first install VEP (https://github.com/Ensembl/ensembl-vep) - install perl module DBI
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf

## VEP installed - need to install vcf2maf now https://github.com/mskcc/vcf2maf

# convert vcf to maf
if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf ]; then
  module load tabix
  perl mskcc-vcf2maf-*/vcf2maf.pl \
    --input-vcf mutect-results/$(basename $TUMOR_BAM .bam)\_muts.vcf \
    --output-maf mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf \
    --tumor-id $(basename $TUMOR_BAM .bam) \
    --normal-id $(basename $NORMAL_BAM .bam) \
    --ref-fasta $HOME/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
fi

# extract HLA muts
if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf ]; then
  head -2 mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf > mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
  grep HLA mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf >> mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
fi

