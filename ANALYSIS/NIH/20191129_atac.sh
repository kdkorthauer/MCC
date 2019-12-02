
cd /rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/

mkdir /rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph

for files in ./*.bigwig
do
  out=${files/.bigwig/.bedGraph}
  out=${out/.\//.\/bedgraph\/}
  echo $out
  bigWigToBedGraph $files $out
done





