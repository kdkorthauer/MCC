#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 40G
#SBATCH -t 0-40:00

cd ../../../PREPROCESS/DNA/
RESDIR=mutect2-results
mkdir -p $RESDIR

module load gatk/4.0.2.1-fasrc01
module load jdk/1.8.0_45-fasrc01
module load samtools
  # gatk Mutect2 --help


################
# build PoN
################


# account for different N GT column name in 5369
if echo "$NORMAL_BAM" | grep -q "5369";then
  NORMAL_BAM2=${NORMAL_BAM/DFCI-5/DFCI5};
else
  NORMAL_BAM2=$NORMAL_BAM;
fi

# create somatic panel of normals
if [ \( -f $RESDIR/DFCI-5367-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5368-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5369-N-01.vcf \) ] ; then
if [ ! \( -f $RESDIR/pon.vcf.gz \) ] ; then
 gatk CreateSomaticPanelOfNormals \
   -vcfs $RESDIR/DFCI-5367-N-01.vcf \
   -vcfs $RESDIR/DFCI-5368-N-01.vcf \
   -vcfs $RESDIR/DFCI-5368-N-01.vcf \
   -O $RESDIR/pon.vcf.gz
fi
fi
fi
fi

# run MuTect if output file does not already exist 
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \) ] ; then
if [ \( -f $RESDIR/pon.vcf.gz \) ] ; then

if [ ! -z "$NORMAL_BAM" ]; then
  # matched normal
  gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I ../../DATA/DNA/$TUMOR_BAM \
     -I ../../DATA/DNA/$NORMAL_BAM \
     -tumor ${TUMOR_BAM%.*} \
     -normal ${NORMAL_BAM2%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --af-of-alleles-not-in-resource 0.00001 \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
else
  # no matched normal
  gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I ../../DATA/DNA/$TUMOR_BAM \
     -tumor ${TUMOR_BAM%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --af-of-alleles-not-in-resource 0.00000005 \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
fi

fi
fi

# filter variants
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf \) ] ; then
  gatk FilterMutectCalls \
    -V $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \
    -O $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf
fi

# extract passing variants
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)_pass\.vcf \) ] ; then
   grep "#\|PASS" $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf > $RESDIR/$(basename $TUMOR_BAM .bam)_pass\.vcf
fi


##################### run vcf2maf
# first install VEP (https://github.com/Ensembl/ensembl-vep) - install perl module DBI
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf

## VEP installed - need to install vcf2maf now https://github.com/mskcc/vcf2maf


# convert vcf to maf
if [ ! -f $RESDIR/$(basename $TUMOR_BAM .bam)\.maf ]; then
  module load tabix
  perl mskcc-vcf2maf-*/vcf2maf.pl \
    --input-vcf $RESDIR/$(basename $TUMOR_BAM .bam)\_pass.vcf \
    --output-maf $RESDIR/$(basename $TUMOR_BAM .bam)\.maf \
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