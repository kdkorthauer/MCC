
###############################
# This script is deprecated - see preprocess_wes.sh instead
###############################

# loop over all tumor / normal pairs 

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
# first first install perl5 with https://gist.github.com/ckandoth/c16b9a423b54dc0a7a37


export PERL_MM_OPT="INSTALL_BASE=$LOCALPERL" 
export PERL_MB_OPT="--install_base $LOCALPERL" 
export PATH="$LOCALPERL/bin:$PATH" cd $HOME/apps/src


In your .bashrc include the following lines:
export LOCALPERL=$HOME/apps/perl
export PERL5LIB=$LOCALPERL:$LOCALPERL/lib/perl5:$PERL5LIB
export PERL_MM_OPT="INSTALL_BASE=$LOCALPERL"
export PERL_MB_OPT="--install_base $LOCALPERL"
export PATH="$LOCALPERL/bin:$PATH"
... still debugging install.pl line - missing modules

# first install VEP with https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697
# next install vcf2maf tool with https://github.com/mskcc/vcf2maf



# CNVkit --> THetA2

# gene annotation downloaded from table browser:
# http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=698008551_5rAuSWReZpqbMZ4sRzebf3aGaii7&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=ensGene&hgta_table=0&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&hgta_outFileName=
# to file annotation/ensGenes_GRCh37.bed.gz - > unzipped
wget -P annotation http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip annotation/refFlat.txt.gz

# download bait bed file
# manually from google drive links to MCC/PREPROCESS/DNA/annotation/
# https://drive.google.com/file/d/1K2sx6rW8iELVWzeJE3reLOh8iyOx0vdU/view
# https://drive.google.com/file/d/1-79TQ7nKvBJtph-OKVznaSNpq15srizM/view

# activate python environment with CNVkit installed
source activate cnvkit

# build accessible regions
cnvkit.py access annotation/GATK_bundle_b37/human_g1k_v37.fasta -o annotation/access-excludes.b37.bed



# Build a reference from normal samples and infer tumor copy ratios
# first sample
cnvkit.py batch ../../DATA/DNA/$TUMOR_BAM --normal ../../DATA/DNA/$NORMAL_BAM \
    --targets annotation/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list \
    --fasta annotation/GATK_bundle_b37/human_g1k_v37.fasta \
    --access annotation/access-excludes.b37.bed \
    --output-reference my_reference.cnn --output-dir cnvkit-results/

#  --annotate annotation/refFlat.txt \ - don't have analogous file for b37 -- chromosome 
# names do not match 

# reuse reference.cnn for subsequent tumor samples
cnvkit.py batch ../../DATA/DNA/$TUMOR_BAM -r cnvkit-results/my_reference.cnn -p 4 --scatter --diagram -d cnvkit-results/