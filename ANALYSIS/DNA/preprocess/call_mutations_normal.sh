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

################
# build PoN
################

# account for different N GT column name in 5369
if echo "$NORMAL_BAM" | grep -q "5369";then
  NORMAL_BAM2=${NORMAL_BAM/DFCI-5/DFCI5};
else
  NORMAL_BAM2=$NORMAL_BAM;
fi

# first call variants on each normal sample
if [ ! \( -f $RESDIR/$(basename $NORMAL_BAM .bam)\.vcf \) ] ; then
  $GATK/gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I ../../DATA/DNA/$NORMAL_BAM \
     -tumor ${NORMAL_BAM2%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     -O $RESDIR/$(basename $NORMAL_BAM .bam)\.vcf
fi
