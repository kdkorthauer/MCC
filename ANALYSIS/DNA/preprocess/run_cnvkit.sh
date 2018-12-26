#!/bin/bash
#SBATCH -J MCCdna
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 20G
#SBATCH -t 0-99:00

cd ../../../PREPROCESS/DNA/

# activate python environment with CNVkit installed
source activate cnvkit

# build accessible regions
if [ ! -f annotation/access-excludes.b37.bed ]; then
  cnvkit.py access annotation/GATK_bundle_b37/human_g1k_v37.fasta -o annotation/access-excludes.b37.bed
fi


# Build a reference from normal samples and infer tumor copy ratios
if [ ! -f cnvkit-results/my_reference.cnn ]; then
  cnvkit.py batch ../../DATA/DNA/*T-01.bam --normal ../../DATA/DNA/*N-01.bam \
    --targets annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list \
    --fasta annotation/GATK_bundle_b37/human_g1k_v37.fasta \
    --access annotation/access-excludes.b37.bed -p 4 \
    --output-reference cnvkit-results/my_reference.cnn --output-dir cnvkit-results/
fi

# reuse reference.cnn for Cell line samples
if [! -f cnvkit-results/DFCI-5369-CL-01.targetcoverage.cnn ]; then
  cnvkit.py batch ../../DATA/DNA/*CL-01.bam -r cnvkit-results/my_reference.cnn -p 4 --scatter --diagram -d cnvkit-results/
fi
