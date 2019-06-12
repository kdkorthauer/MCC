#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 4-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p commons,shared,serial_requeue   # Partition to submit to
#SBATCH --mem=20G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o scratch/slurm/mv_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e scratch/slurm/mv_%j.err  # File to which STDERR will be written, %j inserts jobid

/n/irizarryfs01_backed_up/kkorthauer/softwareTools/globusconnectpersonal-2.3.6/globusconnectpersonal -start &

globus transfer --recursive f6657116-4a1f-11e8-8f88-0a6d4e044368:for_keegan/ 766d20f2-87b4-11e9-8e6a-029d279f7e24:/n/irizarryfs01/kkorthauer/MCC/DATA/NIH/
