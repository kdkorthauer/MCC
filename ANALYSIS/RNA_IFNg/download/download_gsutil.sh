
# where data will be saved
outfile=/rafalab/keegan/MCC/DATA/RNA_IFNg/

# assumes that gsutil is installed and account is authenticated to download the 
# data from the Broad's fire cloud -- requires specific python environment with 
# compiled version of crcmod (installed in env named gcloud2).
# steps to reproduce: 
#    1. download and untar standalone gustil 
#    2. load python 2.7 (module load command)
#    3. create new environment with conda (here called gcloud2)
#    4. install crcmod with pip
# subsequently, need to load appropriate version of python and this environment:
#module load Anaconda
#source activate gcloud2

#gcloud init

fnames="$(cut -f8 /rafalab/keegan/MCC/DATA/RNA_IFNg/sample_filt.tsv)"

rdate="$(cut -f142 /rafalab/keegan/MCC/DATA/RNA_IFNg/sample_filt.tsv)"


for f in ${fnames[@]} 
do
  fout=${f##*/}
  if [ "$f" != "cram_or_bam_path" ] && [ ! -e "$outfile/$fout" ] 
  then
    gsutil cp $f $outfile
  fi
done


fnames="$(cut -f7 /rafalab/keegan/MCC/DATA/RNA_IFNg/sample_filt.tsv)"

for f in $fnames 
do
  fout=${f##*/}
  if [ "$f" != "crai_or_bai_path" ] && [ ! -e "$outfile/$fout" ] 
  then
    gsutil cp $f $outfile
  fi
done
