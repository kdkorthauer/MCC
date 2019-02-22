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

for id in 5367 5368 5369; do
echo "Tumor vs Normal for ${id}"
export NORMAL_BAM=DFCI-${id}-N-01.bam
export TUMOR_BAM=DFCI-${id}-T-01.bam

sbatch -o ../slurm/${id}T.out -e ../slurm/${id}T.err call_mutations_single.sh
sleep 1 

echo "Cell line vs Normal for ${id}"
export NORMAL_BAM=DFCI-${id}-N-01.bam
export TUMOR_BAM=DFCI-${id}-CL-01.bam

sbatch -o ../slurm/${id}C.out -e ../slurm/${id}C.err call_mutations_single.sh
sleep 1 
done


# call mutations - no matched normal
# use 'contemporary normal' as recommended here: https://gatkforums.broadinstitute.org/gatk/discussion/2249/recommended-syntax-for-unmatched-tumor#latest
# using a normal sample from a different individual

for id in 5473 5474; do
echo "Tumor for ${id}"
export NORMAL_BAM=DFCI-5368-N-01.bam
export TUMOR_BAM=DFCI-${id}-T-01.bam

sbatch -o ../slurm/${id}T.out -e ../slurm/${id}T.err call_mutations_single.sh
sleep 1 

echo "Cell line for ${id}"
export NORMAL_BAM=DFCI-5368-N-01.bam
export TUMOR_BAM=DFCI-${id}-C-01.bam

sbatch -o ../slurm/${id}C.out -e ../slurm/${id}C.err call_mutations_single.sh
sleep 1 
done


##################### run CNVkit ############################
sbatch -o ../slurm/${id}cnvkit.out -e ../slurm/${id}cnvkit.err run_cnvkit.sh



