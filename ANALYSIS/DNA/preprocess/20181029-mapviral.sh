#!/bin/bash
#SBATCH -J viralDNA
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p shared,commons,serial_requeue
#SBATCH --mem 20G
#SBATCH -t 0-5:00
#SBATCH -o ../../../slurm/expl-%j.out
#SBATCH -e ../../../slurm/expl-%j.err


## Script to align unmapped DNA reads to viral sequence 

# download Merkel Cell Polyomavirus reference sequence:
cd ../../../PREPROCESS/DNA/viral/ref/
mkdir ref
curl "https://www.ebi.ac.uk/ena/data/view/EU375804&display=fasta" > mcc.fa


# combine and index viral genomes with bwa
module load bwa
cd ref
bwa index mcc.fa
cd ..

# pull out unmapped reads from original BAMs
module load samtools

for i in /n/irizarryfs01/kkorthauer/MCC/DATA/DNA/*1.bam; do
    [ -f "$i" ] || break
    samp=$(basename $i .bam)
    pref=$(dirname $i)
    if [[ ! -e $pref/$samp\_unmapped.bam ]]; then
      samtools view -f 0x04 -h -b $i -o $pref/$samp\_unmapped.bam 
    fi
done


# realign unmapped reads to viral genomes
for i in /n/irizarryfs01/kkorthauer/MCC/DATA/DNA/*unmapped.bam; do
    [ -f "$i" ] || break
    samp=$(basename $i .bam)
    pref=$(dirname $i)
    if [[ ! -e $pref/$samp\_virusmap.bam ]]; then
      bwa aln -b ref/mcc.fa $i > $pref/$samp\_virus.sai; 
      bwa samse ref/mcc.fa $pref/$samp\_virus.sai $i | tee >(samtools view - -Sb -f 0x04 -o $pref/$samp\_virusunmap.bam) | samtools view - -Sb -F 0x04 -o $pref/$samp\_virusmap.bam;
    fi
done

# query which samples had reads aligned to the mcc polyomav
for i in /n/irizarryfs01/kkorthauer/MCC/DATA/DNA/*virusmap.bam; do
    [ -f "$i" ] || break
    samtools view $i | cut -f3 | sort | uniq -c | awk -v i=$i '{print i"\t"$1"\t"$2}' >> virus.read.counts.txt
done