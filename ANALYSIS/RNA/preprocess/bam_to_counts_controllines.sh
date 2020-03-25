#!/bin/bash

######################################################
### parameters 
######################################################
# root directory of RNAseq data
export DATDIR='/scratch/st-kdkortha-1/MCC/PREPROCESS/RNA/control_lines/'
# root directory of RNAseq counts
export COUNTDIR='/scratch/st-kdkortha-1/MCC/PREPROCESS/RNA/control_lines/counts'
# reference genome assembly 
export REFDIR='/scratch/st-kdkortha-1/MCC/PREPROCESS/RNA/reference/'
# number of cores
export CORES=8
######################################################
### end of parameters             
######################################################

SCRIPTDIR=$PWD
mkdir -p $COUNTDIR

## activate conda environ with necessary tools installed
eval "$(conda shell.bash hook)"
conda activate cfmedip

# loop through bam files found and count features with featureCounts
for f in $DATDIR/*/*duplicates_marked.bam
do    
  fout=${f%%????}
  fout=${fout##*/}

  group="$(basename "$(dirname "$f")")"

  outfile=$COUNTDIR/$group\_$fout\.counts.txt
  outfile="$(echo $outfile | sed 's/_Aligned.out_duplicates_marked//g')"
  outfile="$(echo $outfile | sed 's/MKL1-RNAseq-//g')"
  outfile="$(echo $outfile | sed 's/Waga-RNAseq-//g')"

  if [[ ! -e $outfile ]];
  then
	echo "counting features of - $fout"
	featureCounts -T $CORES -g gene_id -a $REFDIR/gencode.v19.annotation_broad.gtf -o $outfile $f
  else
  	echo "count file $outfile already exists"
  fi
done


# R analysis
cat <<EOT > combineCounts_controllines.R
dat.dir <- Sys.getenv("COUNTDIR")
setwd(paste0(dat.dir));
# list all files that end in ".counts.txt";
# but do not contain ".summary" ;
count.files <- list.files(path = ".");
count.files <- count.files[count.files %in% list.files(path = ".", pattern = ".counts.txt")];
count.files <- count.files[!(count.files %in% list.files(path = ".", pattern = "summary"))];
count.files <- count.files[!(count.files %in% list.files(path = ".", pattern = "mRNA-Seq"))];
# build count table (combine all samples into one table);
counts <- NULL;
it <- 1;
for (f in count.files){;
  tab <- read.table(file=f, header = TRUE, stringsAsFactors=FALSE);
  colnames(tab)[ncol(tab)] <- gsub(".counts.txt", "", f)
    if (it > 1){;
      if (all.equal(counts\$Geneid, tab\$Geneid)){;
        counts <- cbind(counts, tab[,ncol(tab)]);
        colnames(counts)[ncol(counts)] <- colnames(tab)[ncol(tab)];
      }else{;
        stop("error: features not identical between count matrices");
      };
     }else{;
       counts <- tab;
     };
    it <- it + 1;
    rm(tab);
};
write.table(counts, file="mRNA-Seq.counts.txt", quote=FALSE, 
   row.names=FALSE, sep="\t")
EOT

R CMD BATCH combineCounts_controllines.R