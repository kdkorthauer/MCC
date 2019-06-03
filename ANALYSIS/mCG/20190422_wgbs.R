# 20190422_wgbs.R

library(readr)
library(dplyr)
library(tidyr)
library(bsseq)
library(dmrseq)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ComplexHeatmap)
library(RColorBrewer)
library(org.Hs.eg.db)

dat.dir <- file.path("/n/irizarryfs01/kkorthauer/MCC/DATA/mCG/wgbs/cov")
dat.files <- list.files(dat.dir, full.names = TRUE)
count.tab <- file.path(dat.dir, "counts.txt") 
bsseq.obj <- file.path(dat.dir, "bsseq.rds")
out.dir <- file.path("../../ANALYSIS/mCG/plots")
dir.create(out.dir, showWarnings=FALSE)

if (!file.exists(count.tab) | !file.exists(bsseq.obj)){
  # get counts table
  if (!file.exists(count.tab)){
	for (f in seq_along(dat.files)){
	  tab <- read_tsv(dat.files[f], col_names = FALSE,
	  	 col_types = list(col_character(), col_double(),
	  	 	              col_double(), col_double(), 
	  	 	              col_double(), col_double())) %>%
	     dplyr::select(-X3,-X4) %>%
	     dplyr::rename(chr = X1,
	     	      pos = X2) %>%
	     arrange(chr, pos)

	  colnames(tab)[3] <- paste0("M_", substr(
	    	gsub("/n/irizarryfs01/kkorthauer/MCC/DATA/mCG/wgbs/cov/1623D-52-", "", dat.files[f]), 1, 2))  
	  colnames(tab)[4] <- paste0("Cov_", substr(
	    	gsub("/n/irizarryfs01/kkorthauer/MCC/DATA/mCG/wgbs/cov/1623D-52-", "", dat.files[f]), 1, 2))

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
    pData(bs)$mcc.id <- c("275", "282", "290", "301", "320",
                          "336", "350", "367", "2314")
    pData(bs)$ifng <- "++"
    pData(bs)$ifng[grepl("336|350", pData(bs)$mcc.id)] <- "-"
    pData(bs)$ifng[grepl("2314|320", pData(bs)$mcc.id)] <- "+"
  }
}

# filter (at least 4 samples have coverage)
ix <- which(rowSums(getCoverage(bs, type = "Cov") > 0) >= 4)
bs <- bs[ix,]

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


plotEmpiricalDistribution(bs[idx,], 
                          testCovariate = "ifng",
                          type = "Cov")
plotEmpiricalDistribution(bs[idx,], 
                          bySample = TRUE,
                          testCovariate = "ifng",
                          adj = 3)
plotEmpiricalDistribution(bs[idx,], 
                          testCovariate = "ifng",
                          adj = 3)


goi <- read_tsv("../genes_of_interest.txt", col_names = FALSE)
colnames(goi) <- c("symbol", "description")

# get coords of genes of interest 
db <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx = transcriptsBy(db)
prom = promoters(tx)


x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
names(prom) <- unlist(xx[names(prom)])
seqlevelsStyle(prom) <- "NCBI"
prom = unlist(prom[names(prom) %in% goi$symbol])
prom$symbol <- names(prom)

# overlap with meth
bs.prom <- subsetByOverlaps(bs, prom)
prom.bs <- subsetByOverlaps(prom, bs)

overlaps <- findOverlaps(bs, prom)
signal <- bs[overlaps@from]
averagedSignal <- aggregate(getMeth(bs, type="raw")[overlaps@from,], 
                            list(overlaps@to), mean, na.rm = TRUE)

averagedSignal <- aggregate(averagedSignal[,-1],
                            list(prom.bs$symbol), mean, na.rm = TRUE)
rownames(averagedSignal) <- averagedSignal$Group.1

# heatmap of goi with viral status and ifng labels
ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) 
virus <- c("lightgrey", "black")
ifng <- c("darkred", "darkblue", "lightblue")
names(virus) <- c("Positive", "Negative")
names(ifng) <- c("-", "++", "+")

# all goi
goi$description[grepl("DC", goi$description)] <- "DC/Mono/Macro"
goi$description[grepl("antigen", goi$description)] <- "MCPyV antigen"
goi$description[grepl("Cytokine", goi$description)] <- "Cytokine"
goi$description[grepl("cell", goi$description)] <- "B/T cell"

ha_column = HeatmapAnnotation(df = data.frame(Virus = pData(bs)$virus,
                                              IFNg = pData(bs)$ifng),
                              col = list(Virus = virus,
                                         IFNg = ifng))

x <- match(goi$symbol, averagedSignal$Group.1)
x <- x[!is.na(x)]
rowcol <- colorRampPalette( rev(brewer.pal(length(unique(goi$description[x])),
             "Paired")) )(length(unique(goi$description[x])))
names(rowcol) <- sort(unique(goi$description[x]))
ha_row = rowAnnotation(df = data.frame(Category = goi$description[x]),
                       col = list(Category = rowcol))


ht = Heatmap(averagedSignal[,-1], name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = FALSE, show_column_names = FALSE)

heatmap.file <- file.path(out.dir, "heatmap_all.pdf")

pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# just class I
y <- goi$description[x] == "Class I"
ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE)
ha_row = rowAnnotation(df = data.frame(Category = goi$description[x][y]),
                       col = list(Category = rowcol))

heatmap.file <- file.path(out.dir, "heatmap_classI.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in virus pos vs neg
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Positive")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Negative"))
  )
byGroupDiff <- abs(byGroupAvg[,1] - byGroupAvg[,2])
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$virus == "Positive"),
                                which(pData(bs)$virus == "Negative")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_classI_orderedV.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in IFNg category
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "++")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "+")),
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "-"))
  )
byGroupDiff <- rowSums(cbind(abs(byGroupAvg[,1] - byGroupAvg[,2]),
                     abs(byGroupAvg[,2] - byGroupAvg[,3]),
                     abs(byGroupAvg[,1] - byGroupAvg[,3])))
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$ifng == "++"),
                                which(pData(bs)$ifng == "+"),
                                which(pData(bs)$ifng == "-")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_classI_orderedI.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# just class II
y <- goi$description[x] == "Class II"
ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE)
ha_row = rowAnnotation(df = data.frame(Category = goi$description[x][y]),
                       col = list(Category = rowcol))

heatmap.file <- file.path(out.dir, "heatmap_classII.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in virus pos vs neg
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Positive")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Negative"))
  )
byGroupDiff <- abs(byGroupAvg[,1] - byGroupAvg[,2])
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$virus == "Positive"),
                                which(pData(bs)$virus == "Negative")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_classII_orderedV.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in IFNg category
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "++")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "+")),
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "-"))
  )
byGroupDiff <- rowSums(cbind(abs(byGroupAvg[,1] - byGroupAvg[,2]),
                     abs(byGroupAvg[,2] - byGroupAvg[,3]),
                     abs(byGroupAvg[,1] - byGroupAvg[,3])))
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$ifng == "++"),
                                which(pData(bs)$ifng == "+"),
                                which(pData(bs)$ifng == "-")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_classII_orderedI.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()


# other
y <- !(goi$description[x] %in% c("Class I", "Class II"))
ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE)
ha_row = rowAnnotation(df = data.frame(Category = goi$description[x][y]),
                       col = list(Category = rowcol))

heatmap.file <- file.path(out.dir, "heatmap_other.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in virus pos vs neg
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Positive")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$virus == "Negative"))
  )
byGroupDiff <- abs(byGroupAvg[,1] - byGroupAvg[,2])
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$virus == "Positive"),
                                which(pData(bs)$virus == "Negative")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_other_orderedV.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# ordered by dmr signal in IFNg category
byGroupAvg <- cbind(rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "++")), 
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "+")),
  rowMeans(subset(averagedSignal[y,-1], select = pData(bs)$ifng == "-"))
  )
byGroupDiff <- rowSums(cbind(abs(byGroupAvg[,1] - byGroupAvg[,2]),
                     abs(byGroupAvg[,2] - byGroupAvg[,3]),
                     abs(byGroupAvg[,1] - byGroupAvg[,3])))
diffRank <- order(byGroupDiff, decreasing = TRUE)

ht = Heatmap(averagedSignal[y,-1], 
               name = "Promoter\nmCG", 
               top_annotation = ha_column, col = ecolors,
               show_row_names = TRUE, show_column_names = FALSE,
               row_order = diffRank,
               column_order = c(which(pData(bs)$ifng == "++"),
                                which(pData(bs)$ifng == "+"),
                                which(pData(bs)$ifng == "-")),
               cluster_rows = FALSE,
               cluster_columns = FALSE)
heatmap.file <- file.path(out.dir, "heatmap_other_orderedI.pdf")
pdf(heatmap.file)
   draw(ht + ha_row)
dev.off()

# traceplot meth colored by viral status

 
# traceplot meth colored by ifng status


