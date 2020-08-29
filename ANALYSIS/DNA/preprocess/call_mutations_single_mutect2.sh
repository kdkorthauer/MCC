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
GATK=""
#module load jdk/1.8.0_45-fasrc01
module load gatk
module load samtools
module load bcftools
module load perl
  # gatk Mutect2 --help

ulimit -c unlimited

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
if [ \( -f $RESDIR/DFCI-5473-N-01.vcf \) ] ; then
if [ \( -f $RESDIR/DFCI-5474-N-01.vcf \) ] ; then
if [ ! \( -f $RESDIR/pon.vcf.gz \) ] ; then
  
  gatk GenomicsDBImport \
   -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
   -V $RESDIR/DFCI-5367-N-01.vcf \
   -V $RESDIR/DFCI-5368-N-01.vcf \
   -V $RESDIR/DFCI-5369-N-01.vcf \
   -V $RESDIR/DFCI-5350-N-01.vcf \
   -V $RESDIR/DFCI-5351-N-01.vcf \
   -V $RESDIR/DFCI-5473-N-01.vcf \
   -V $RESDIR/DFCI-5474-N-01.vcf \
   --genomicsdb-workspace-path pon_db_pdx \
   --merge-input-intervals \
   -L annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list

 gatk CreateSomaticPanelOfNormals \
   -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
   --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
   -V gendb://pon_db \
   -O $RESDIR/pon.vcf.gz
fi
fi
fi
fi
fi
fi
fi
fi 

# index bam
if [ ! \( -f /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$TUMOR_BAM\.bai \) ] ; then
  samtools index /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$TUMOR_BAM
fi

# run MuTect if output file does not already exist 
# no af filter needed: https://gatkforums.broadinstitute.org/gatk/discussion/comment/57078
# need -genotype-germline-sites flag to run with PureCN (https://www.bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/PureCN.pdf)
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \) ] ; then
if [ \( -f $RESDIR/pon.vcf.gz \) ] ; then

if [ ! -z "$NORMAL_BAM" ]; then
  # matched normal
  gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$TUMOR_BAM \
     -I /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$NORMAL_BAM \
     -tumor ${TUMOR_BAM%.*} \
     -normal ${NORMAL_BAM2%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     -L annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
else
  # no matched normal
  gatk Mutect2 \
     -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
     -I /rafalab/.keegan/MCC/PREPROCESS/DNA/pdx_filter/cleanbam/$TUMOR_BAM \
     -tumor ${TUMOR_BAM%.*} \
     --germline-resource annotation/GATK_bundle_b37/af-only-gnomad.raw.sites.b37.vcf \
     --panel-of-normals $RESDIR/pon.vcf.gz \
     -L annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list \
     -O $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf
fi
fi
fi

# filter variants
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf \) ] ; then
  gatk FilterMutectCalls \
    -R annotation/GATK_bundle_b37/human_g1k_v37.fasta \
    -V $RESDIR/$(basename $TUMOR_BAM .bam)\.vcf \
    -O $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf
fi

### REVISIT step below -- not sure if necessary
# after upgrading to version 4.1.2.0

# extract passing variants
if [ ! \( -f $RESDIR/$(basename $TUMOR_BAM .bam)_pass\.vcf \) ] ; then
   grep "#\|PASS" $RESDIR/$(basename $TUMOR_BAM .bam)_filt\.vcf > $RESDIR/$(basename $TUMOR_BAM .bam)_pass\.vcf
fi


##################### run vcf2maf
# VEP version 95 (VER=95 in gist https://gist.github.com/ckandoth/5390e3ae4ecf182fa92f6318cfa9fa97)
# first install VEP (https://github.com/Ensembl/ensembl-vep) - install perl module DBI
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf

#PLEASE REMEMBER TO 
#1. add /homes/keegan/vep to your PERL5LIB environment variable
#2. add /homes/keegan/vep/htslib to your PATH environment variable

#export VEP_PATH=$HOME/vep
#export VEP_DATA=/rafalab/keegan/ReferenceGenomes/.vep
#export VER=95

#curl -LO https://github.com/Ensembl/ensembl-vep/archive/release/95.3.tar.gz
#tar -zxf 95.3.tar.gz; rm -f 95.3.tar.gz; mv ensembl-vep-release-95.3 $VEP_PATH
#cd $VEP_PATH

#export PERL5LIB=$VEP_PATH:$PERL5LIB
#export PATH=$VEP_PATH/htslib:$PATH

# convert vcf to maf
if [ ! -f $RESDIR/$(basename $TUMOR_BAM .bam)\.maf ]; then
  #module load tabix
  perl mskcc-vcf2maf-5*/vcf2maf.pl \
    --input-vcf $RESDIR/$(basename $TUMOR_BAM .bam)\_pass.vcf \
    --output-maf $RESDIR/$(basename $TUMOR_BAM .bam)\.maf \
    --tumor-id $(basename $TUMOR_BAM .bam) \
    --normal-id $(basename $NORMAL_BAM .bam) \
    --vcf-normal-id $(basename $NORMAL_BAM2 .bam) \
    --ref-fasta /rafalab/.keegan/ReferenceGenomes/.vep/homo_sapiens/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    --vep-path ~/vep \
    --vep-data /rafalab/.keegan/ReferenceGenomes/.vep \
    --filter-vcf /rafalab/.keegan/ReferenceGenomes/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --vep-forks 1
fi

# extract HLA muts
#if [ ! -f mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf ]; then
#  head -2 mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf > mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#  grep HLA mutect-results/$(basename $TUMOR_BAM .bam)\_muts.maf >> mutect-results/$(basename $TUMOR_BAM .bam)\_HLA_muts.maf
#fi
