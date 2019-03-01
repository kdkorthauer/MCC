#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 40G
#SBATCH -t 0-40:00

cd ../../../PREPROCESS/DNA/
RESDIR=mutect-results
mkdir -p $RESDIR


module load centos6/0.0.1-fasrc01
module load java/1.7.0_60-fasrc01

# run MuTect if output file does not already exist 

if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \) ] ; then
  java -Xmx2g -jar mutect-src/mutect/target/mutect-1.1.7.jar \
    --analysis_type MuTect \
    --reference_sequence annotation/GATK_bundle_b37/human_g1k_v37.fasta \
    --dbsnp annotation/GATK_bundle_b37/dbsnp_138.b37.vcf \
    --input_file:normal ../../DATA/DNA/$NORMAL_BAM \
    --input_file:tumor ../../DATA/DNA/$TUMOR_BAM \
    --out $RESDIR/$(basename $TUMOR_BAM .bam)\_call_stats.txt \
    --coverage_file $RESDIR/$(basename $TUMOR_BAM .bam)\_coverage.wig.txt \
    -vcf $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf 
fi

# for testing short regions use, e.g.
# --intervals 1:10000000-10100000 \

# on output, extract PASS mutations in separate vcf file
if [ ! -f $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.vcf ]; then
  grep -v REJECT $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf  > $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.vcf 
fi


##################### run vcf2maf
# first install VEP (https://github.com/Ensembl/ensembl-vep) - install perl module DBI
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf

## VEP installed - need to install vcf2maf now https://github.com/mskcc/vcf2maf

# account for different N GT column name in 5369
TVAR=$(grep DFCI5 $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.vcf);

if [ ! -z "$TVAR" ]; then
  NORMAL_BAM2=${NORMAL_BAM/DFCI-5/DFCI5};
else
  NORMAL_BAM2=$NORMAL_BAM;
fi

# convert vcf to maf
if [ ! -f $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.maf ]; then
  module load tabix
  perl mskcc-vcf2maf-*/vcf2maf.pl \
    --input-vcf $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.vcf \
    --output-maf $RESDIR/$(basename $TUMOR_BAM .bam)\_muts.maf \
    --tumor-id $(basename $TUMOR_BAM .bam) \
    --normal-id $(basename $NORMAL_BAM .bam) \
    --vcf-normal-id $(basename $NORMAL_BAM2 .bam) \
    --ref-fasta $HOME/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
fi

# extract HLA muts
#if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf ]; then
#  head -2 mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf > mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#  grep HLA mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf >> mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#fi
