#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-10:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p commons,shared,serial_requeue   # Partition to submit to
#SBATCH --mem=20G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../_slurm/downloaddropbox_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../_slurm/downloaddropbox_%j.err  # File to which STDERR will be written, %j inserts jobid

module load rclone
cd /n/irizarryfs01/kkorthauer/MCC/DATA/mCG/wgbs

rclone copy "dropbox:MCC/Bismark Bisulfite Data - Ben Singer - wgbs" .
# wgbs in cov dir
