#!/bin/bash
#SBATCH -J ATAC
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem 20G
#SBATCH -t 0-5:00
#SBATCH -o slurm/MCCrna-%j.out
#SBATCH -e slurm/MCCrna-%j.err


cd /rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/

mkdir -p /rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph

for files in ./*.bigwig
do
  out=${files/.bigwig/.bedGraph}
  out=${out/.\//.\/bedgraph\/}
  if [ ! -f "$out" ]; then
    echo $out
    bigWigToBedGraph $files $out
  fi
done

R CMD BATCH --quiet --no-restore --no-save /rafalab/keegan/MCC/ANALYSIS/NIH/20191129_atac.R /rafalab/keegan/MCC/ANALYSIS/NIH/20191129_atac.Rout


