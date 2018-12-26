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
      if (all.equal(counts$Geneid, tab$Geneid)){;
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
