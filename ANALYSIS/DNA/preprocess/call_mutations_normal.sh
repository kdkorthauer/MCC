#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 40G
#SBATCH -t 0-40:00

cd ../../../PREPROCESS/DNA/
RESDIR=mutect2-gatk4.1.2.0-results-xengsort
mkdir -p $RESDIR
#GATK=/n/irizarryfs01_backed_up/kkorthauer/softwareTools/gatk
module load gatk
GATK=""
#module load jdk/1.8.0_45-fasrc01

################
# build PoN
################

# account for different N GT column name in 5369
if echo "$NORMAL_BAM" | grep -q "5369";then
  NORMAL_BAM2=${NORMAL_BAM/DFCI-5/DFCI5};
else
  NORMAL_BAM2=$NORMAL_BAM;
fi

ulimit -c unlimited

# first call variants on each normal sample
if [ ! \( -f /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$NORMAL_BAM\.bai \) ] ; then
  samtools index /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$NORMAL_BAM
fi

if [ ! \( -f $RESDIR/$(basename $NORMAL_BAM .bam)\.vcf \) ] ; then
  gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$NORMAL_BAM \
     -tumor ${NORMAL_BAM2%.*} \
     --max-mnp-distance 0 \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     -L annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list \
     -O $RESDIR/$(basename $NORMAL_BAM .bam)\.vcf
fi
