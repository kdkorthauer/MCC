# this script calls the "call_mutations_single.sh" script to call somatic 
# mutations using mutect. Next, it calls "run_cnvkit.sh" to run CNVkit.

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


# convert vcf to maf
# first install VEP with https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697
# ^ run into problems install.pl line still missing modules
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf


##################### run CNVkit ############################
sbatch -o ../slurm/${id}cnvkit.out -e ../slurm/${id}cnvkit.err run_cnvkit.sh

##################### run vcf2maf
# first install VEP (https://github.com/Ensembl/ensembl-vep) - install perl module DBI
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf

## VEP installed - need to install vcf2maf now https://github.com/mskcc/vcf2maf

# and add the following to bashrc
export PERL_MM_OPT="INSTALL_BASE=$LOCALPERL" 
export PERL_MB_OPT="--install_base $LOCALPERL" 
export PATH="$LOCALPERL/bin:$PATH" cd $HOME/apps/src


In your .bashrc include the following lines:
export LOCALPERL=$HOME/apps/perl
export PERL5LIB=$LOCALPERL:$LOCALPERL/lib/perl5:$PERL5LIB
export PERL_MM_OPT="INSTALL_BASE=$LOCALPERL"
export PERL_MB_OPT="--install_base $LOCALPERL"
export PATH="$LOCALPERL/bin:$PATH"