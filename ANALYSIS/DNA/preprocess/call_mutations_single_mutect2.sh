#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 40G
#SBATCH -t 0-40:00

cd ../../../PREPROCESS/DNA/
RESDIR=mutect2-gatk4-results
mkdir -p $RESDIR
GATK=/n/irizarryfs01_backed_up/kkorthauer/softwareTools/gatk

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
if [ \( -f $RESDIR/DFCI-5350-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5351-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5367-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5368-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5369-N-01.vcf \) ] ; then
if [ ! \( -f $RESDIR/pon.vcf.gz \) ] ; then
 # use gatk3 to create PON - new version req splitting over all contigs 
 # see if there is a solution after 4.1.1 is released (should be late march 2019)
 # https://github.com/broadinstitute/gatk/pull/5675
 module load gatk
 gatk CreateSomaticPanelOfNormals \
   -vcfs $RESDIR/DFCI-5367-N-01.vcf \
   -vcfs $RESDIR/DFCI-5368-N-01.vcf \
   -vcfs $RESDIR/DFCI-5368-N-01.vcf \
   -O $RESDIR/pon.vcf.gz
fi
fi
fi
fi
fi
fi

# run MuTect if output file does not already exist 
# no af filter needed: https://gatkforums.broadinstitute.org/gatk/discussion/comment/57078
# need -genotype-germline-sites flag to run with PureCN (https://www.bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/PureCN.pdf)
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \) ] ; then
if [ \( -f $RESDIR/pon.vcf.gz \) ] ; then

if [ ! -z "$NORMAL_BAM" ]; then
  # matched normal
  $GATK/gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I ../../DATA/DNA/$TUMOR_BAM \
     -I ../../DATA/DNA/$NORMAL_BAM \
     -tumor ${TUMOR_BAM%.*} \
     -normal ${NORMAL_BAM2%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     --genotype-germline-sites \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
else
  # no matched normal
  $GATK/gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I ../../DATA/DNA/$TUMOR_BAM \
     -tumor ${TUMOR_BAM%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     --genotype-germline-sites \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
fi

fi
fi

# filter variants
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf \) ] ; then
  $GATK/gatk FilterMutectCalls \
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
    --ref-fasta $HOME/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --vep-forks 1
fi

# extract HLA muts
#if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf ]; then
#  head -2 mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf > mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#  grep HLA mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf >> mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#fi
