# where data will be saved
outfile=/rafalab/keegan/MCC/DATA/DNA/

gsutil cp -r gs://fc-7303af08-c226-4422-99df-ea1673266748 $outfile


for dir in */; do
  file=$(find $dir -type f -name "*mutect.filtered.nonsyn.annotated")
  size=$(ls -l $file)
  echo $(du -h $file)
done
# no nonsyn mutations

for dir in */; do
  file=$(find $dir -type f -name "*mutect.filtered.syn.annotated")
  size=$(ls -l $file)
  echo $(du -h $file)
done
# no syn mutations

for dir in */; do
  file=$(find $dir -type f -name "*mutect.ambiguous.annotated")
  size=$(ls -l $file)
  echo $(du -h $file)
done
# no ambiguous mutations

touch $outfile/polysolverHLAType.txt

full_names=$(tail -n +2 $outfile/pair_polysolver.tsv | cut -f2)

# print out HLA types

polysolverWorkflow/d13fcf67-a537-4149-9b04-83c0cb6fccda/call-polysolverAnnot/h

for f in $full_names; do
  dir=$(echo $f | awk -F 'polysolverWorkflow/' '{print $2}')
  dir=$(echo $dir | awk -F '/call-polysolverAnnot' '{print $1}')

  logfile=$(find $dir -type f -name "*polysolverType.log")
  hlafile=$(find $dir -type f -name "*winners.hla.txt")

  name=$(echo $f | awk -F 'hla_annot_out/' '{print $2}')
  name=$(echo $name | awk -F '.mutect.ambiguous.annotated' '{print $1}')

  echo $name >> $outfile/polysolverHLAType.txt
  cat $hlafile >> $outfile/polysolverHLAType.txt
done

