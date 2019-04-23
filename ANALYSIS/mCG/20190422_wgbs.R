# 20190422_wgbs.R

library(readr)
library(dplyr)
library(tidyr)
library(bsseq)
library(dmrseq)

dat.dir <- file.path("../../DATA/mCG/wgbs/cov")
dat.files <- list.files(dat.dir, full.names = TRUE)
count.tab <- file.path(dat.dir, "counts.txt") 
bsseq.obj <- file.path(dat.dir, "bsseq.rds")


if (!file.exists(count.tab) | !file.exists(bsseq.obj)){
  # get counts table
  if (!file.exists(count.tab)){
	for (f in seq_along(dat.files)){
	  tab <- read_tsv(dat.files[f], col_names = FALSE,
	  	 col_types = list(col_character(), col_double(),
	  	 	              col_double(), col_double(), 
	  	 	              col_double(), col_double())) %>%
	     select(-X3,-X4) %>%
	     rename(chr = X1,
	     	    pos = X2) %>%
	     arrange(chr, pos)

	  colnames(tab)[3] <- paste0("M_", substr(
	    	gsub("../../DATA/mCG/wgbs/cov/1623D-52-", "", dat.files[f]), 1, 2))  
	  colnames(tab)[4] <- paste0("Cov_", substr(
	    	gsub("../../DATA/mCG/wgbs/cov/1623D-52-", "", dat.files[f]), 1, 2))

	  if (f > 1){
	    tab <- full_join(tab_prev, tab, by = c("chr", "pos")) 
	  } 

	  tab_prev <- tab
	}

	rm(tab_prev)
  	saveRDS(tab, file = count.tab)
  }else{
  	tab <- readRDS(count.tab)
  }

  # build bsseq obj
  if(!file.exists(bsseq.obj)){
    # replace NAs with zeros
    tab[is.na(tab)] <- 0

    chr <- tab$chr
    pos <- tab$pos
    meth.mat <- tab[,c(grepl("M", colnames(tab)))]
    cov.mat <- tab[,c(grepl("Cov", colnames(tab)))]
    colnames(meth.mat) <- colnames(cov.mat) <- gsub("M_", "", colnames(meth.mat))
    rm(tab)

    # filter out loci with zero cov
    ix <- which(rowSums(cov.mat) == 0)
    meth.mat <- meth.mat[-ix,]
    cov.mat <- cov.mat[-ix,]
    cov.mat <- cov.mat + meth.mat
    chr <- chr[-ix]
    pos <- pos[-ix]

    bs <- BSseq(M = as.matrix(meth.mat),
    	        Cov = as.matrix(cov.mat),
    	        pos = pos,
    	        chr = chr)

    # collapse strands
    library(BSgenome.Hsapiens.UCSC.hg38)  
	chrom <- names(Hsapiens)[1:25]
	cgs <- lapply(chrom, function(x) start(matchPattern("CG", Hsapiens[[x]])))

	cpgr <- do.call(c, lapply(1:25, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
	seqlevels(cpgr) <- gsub("chr", "", seqlevels(cpgr))
	seqlevels(cpgr) <- gsub("M", "MT", seqlevels(cpgr))
	end(cpgr) <- start(cpgr) # represent only one strand

	pos <- findOverlaps(bs, cpgr)
	strand(bs) <- "-"
	strand(bs[pos@from]) <- "+"
    bs <- strandCollapse(bs)

    pData(bs)$sample <- colnames(meth.mat)
    pData(bs)$virus <- "Positive"
    pData(bs)$virus[grepl("2|3|5|7", pData(bs)$sample)] <- "Negative"

    saveRDS(bs, file = bsseq.obj)
  }else{
  	bs <- readRDS(bsseq.obj)
  }
}



set.seed(486)
idx <- sample(1:nrow(bs), 1e6)
plotEmpiricalDistribution(bs[idx,], 
                          testCovariate = "virus",
                          type = "Cov")
plotEmpiricalDistribution(bs[idx,], 
                          bySample = TRUE,
                          testCovariate = "virus",
                          adj = 3)
plotEmpiricalDistribution(bs[idx,], 
                          testCovariate = "virus",
                          adj = 3)



