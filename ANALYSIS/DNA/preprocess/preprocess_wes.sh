# this script calls the "call_mutations_single.sh" script to call somatic 
# mutations using mutect. Then it annotates mutations using VEP and the
# vcf2maf tool. Next, it calls "run_cnvkit.sh" to run CNVkit.


# convert vcf to maf
# first install VEP with https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf:
module load samtools

# only run once to install vcf2maf. Assumes VEP is already installed.
if [ ! -f ../../../PREPROCESS/DNA/mskcc-vcf2maf-*/vcf2maf.pl  ]; then
  export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
  curl -L -o ../../PREPROCESS/DNA/mskcc-vcf2maf.tar.gz $VCF2MAF_URL
  tar -xvf ../../PREPROCESS/DNA/mskcc-vcf2maf.tar.gz -C ../../PREPROCESS/DNA/
  chmod +x ../../PREPROCESS/DNA/mskcc-vcf2maf-*/*.pl
fi

# call mutations - loop over all tumor / normal pairs 

for id in 5350 5351 5367 5368 5369; do

echo "Normal for ${id}"
export NORMAL_BAM=DFCI-${id}-N-01.bam
sbatch -o ../slurm/${id}N-std.out -e ../slurm/${id}N-std.err call_mutations_normal.sh
sleep 1 

done

### wait until above is complete before running next steps...

for id in 5350 5351 5367 5368 5369; do
echo "Tumor vs Normal for ${id}"
export NORMAL_BAM=DFCI-${id}-N-01.bam
export TUMOR_BAM=DFCI-${id}-T-01.bam

sbatch -o ../slurm/${id}T-std.out -e ../slurm/${id}T-std.err call_mutations_single_mutect2.sh
sleep 1 

echo "Cell line vs Normal for ${id}"
export NORMAL_BAM=DFCI-${id}-N-01.bam
export TUMOR_BAM=DFCI-${id}-CL-01.bam

if [ $id = "5350" ] || [ $id = "5351" ]; then
  export TUMOR_BAM=DFCI-${id}-C-01.bam
fi

sbatch -o ../slurm/${id}C-std.out -e ../slurm/${id}C-std.err call_mutations_single_mutect2.sh
sleep 1 
done


# call mutations - no matched normal

for id in 5473 5474; do
echo "Tumor for ${id}"
export NORMAL_BAM=
export TUMOR_BAM=DFCI-${id}-T-01.bam

sbatch -o ../slurm/${id}T-std.out -e ../slurm/${id}T-std.err call_mutations_single_mutect2.sh
sleep 1 

echo "Cell line for ${id}"
export NORMAL_BAM=
export TUMOR_BAM=DFCI-${id}-C-01.bam

sbatch -o ../slurm/${id}C-std.out -e ../slurm/${id}C-std.err call_mutations_single_mutect2.sh
sleep 1 
done



