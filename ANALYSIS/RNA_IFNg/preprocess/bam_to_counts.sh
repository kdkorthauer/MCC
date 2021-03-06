#!/bin/bash

######################################################
### parameters 
######################################################
# root directory of RNAseq data
export DATDIR='/rafalab/keegan/MCC/DATA/RNA_IFNg/'
# root directory of RNAseq counts
export COUNTDIR='/rafalab/keegan/MCC/PREPROCESS/RNA_IFNg'
# reference genome assembly 
export REFDIR='/rafalab/keegan/ReferenceGenomes/human'
# number of cores
export CORES=8
######################################################
### end of parameters             
######################################################

SCRIPTDIR=$PWD



# loop through bam files found and count features with featureCounts
for f in $DATDIR/*.bam
do    
  fout=${f%%????}
  fout=${fout##*/}
  outfile=$COUNTDIR/$fout\.counts.txt
  if [[ ! -e $outfile ]];
  then
	echo "counting features of - $fout"
	~/bin/subread-2.0.0-Linux-x86_64/bin/featureCounts -T $CORES -g gene_id -a $REFDIR/gencode.v19.annotation.gtf -p -s 2 -o $outfile $f
  else
  	echo "count file $outfile already exists"
  fi
done


# R analysis
cat <<EOT > combineCounts.R
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

R CMD BATCH combineCounts.R