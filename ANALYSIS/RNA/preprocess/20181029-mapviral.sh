#!/bin/bash
#SBATCH -J viralDNA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 20G
#SBATCH -t 0-5:00
#SBATCH -o ../slurm/expl-%j.out
#SBATCH -e ../slurm/expl-%j.err

## Script to align unmapped RNAseq reads to viral sequence 

# download Merkel Cell Polyomavirus reference sequence:
cd ../../../PREPROCESS/RNA/viral
mkdir ref
curl "https://www.ebi.ac.uk/ena/data/view/EU375804&display=fasta" > ref/mcc.fa


# combine and index viral genomes with bwa
module load bwa
cd ref
bwa index mcc.fa
cd ..

# pull out unmapped reads from original BAMs
module load samtools
module load bedtools2

for i in /n/irizarryfs01_backed_up/kkorthauer/MCC/DATA/RNA/*1.bam; do
    [ -f "$i" ] || break
    samp=$(basename $i .bam)
    pref=$(dirname $i)
    if [[ ! -e $pref/$samp\_unmapped.bam ]]; then
      samtools view -f 0x04 -h -b $i -o $pref/$samp\_unmapped.bam 
    fi
done


# realign unmapped reads to viral genomes
for i in /n/irizarryfs01_backed_up/kkorthauer/MCC/DATA/RNA/*unmapped.bam; do
    [ -f "$i" ] || break
    samp=$(basename $i .bam)
    pref=$(dirname $i)
    if [[ ! -e $pref/$samp\_virusmap.bam ]]; then
      bwa aln -b ref/mcc.fa $i > $pref/$samp\_virus.sai; 
      bwa samse ref/mcc.fa $pref/$samp\_virus.sai $i | tee >(samtools view - -Sb -f 0x04 -o $pref/$samp\_virusunmap.bam) | samtools view - -Sb -F 0x04 -o $pref/$samp\_virusmap.bam;
    fi
done

# bam to bedGraph
# 1. sort bam
# 2. index bam 
# 3. convert to bedGraph

for i in /n/irizarryfs01_backed_up/kkorthauer/MCC/DATA/RNA/*virusmap.bam; do
    samp=$(basename $i)
    pref=$(dirname $i)
    if [[ ! -e $pref/$samp\.sorted.bam ]]; then
      samtools sort $pref/$samp > $pref/$samp\.sorted.bam
    fi 

    if [[ ! -e $pref/$samp\.sorted.bam.bai ]]; then
      samtools index $pref/$samp\.sorted.bam 
    fi 

    if [[ ! -e $pref/$samp\.sorted.bam.bedGraph ]]; then
      bedtools genomecov -ibam $pref/$samp\.sorted.bam -bg > $pref/$samp\.sorted.bam.bedGraph
    fi 
done


# query which samples had reads aligned to the mcc polyomav
for i in /n/irizarryfs01_backed_up/kkorthauer/MCC/DATA/RNA/*virusmap.bam; do
    [ -f "$i" ] || break
    samtools view $i | cut -f3 | sort | uniq -c | awk -v i=$i '{print i"\t"$1"\t"$2}' >> virus.read.counts.txt
done