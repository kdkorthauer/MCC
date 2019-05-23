
# where data will be saved
outfile=/n/irizarryfs01/kkorthauer/MCC/DATA/DNA/

# assumes that gsutil is installed and account is authenticated to download the 
# data from the Broad's fire cloud -- requires specific python environment with 
# compiled version of crcmod (installed in env named gcloud2).
# steps to reproduce: 
#    1. download and untar standalone gustil 
#    2. load python 2.7 (module load command)
#    3. create new environment with conda (here called gcloud2)
#    4. install crcmod with pip
# subsequently, need to load appropriate version of python and this environment:
module load Anaconda/5.0.1-fasrc02
source activate gcloud2

fnames="$(cut -f8 /n/irizarryfs01/kkorthauer/MCC/DATA/DNA/sample.tsv)"

for f in $fnames 
do
  fout=${f##*/}
  if [ "$f" != "cram_or_bam_path" ] && [ ! -e "$outfile/$fout" ] 
  then
    gsutil cp $f $outfile
  fi
done

fnames="$(cut -f7 /n/irizarryfs01/kkorthauer/MCC/DATA/DNA/sample.tsv)"

# also get bais
for f in $fnames 
do
  fout=${f##*/}
  if [ "$f" != "crai_or_bai_path" ] && [ ! -e "$outfile/$fout" ] 
  then
    gsutil cp $f $outfile
  fi
done


# now add "chr" to chromosome names (to be compatible with GATK ref fasta)
module load samtools

for file in *.bam
do
  filename=`echo $file | cut -d "." -f 1`
  if [ ! -e "$outfile/${filename}_chr.bam" ]
  samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > $outfile/${filename}_chr.bam 
done
