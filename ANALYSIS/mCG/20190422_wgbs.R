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
    chr <- chr[-ix]
    pos <- pos[-ix]

    bs <- BSseq(M = as.matrix(meth.mat),
    	        Cov = as.matrix(cov.mat),
    	        pos = pos,
    	        chr = chr)

    saveRDS(bs, file = bsseq.obj)
  }else{
  	bs <- readRDS(bsseq.obj)
  }
}



M <- matrix(0:8, 3, 3)
Cov <- matrix(1:9, 3, 3)
BS1 <- BSseq(chr = c("chr1", "chr2", "chr1"), pos = c(1,2,3),M = M, Cov = Cov, sampleNames = c("A","B", "C"))