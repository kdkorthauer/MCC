library(Gviz)
library(data.table)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rlist)
library(AnnotationDbi)
library(genomation)

gen <- "hg38"
gtrack <- GenomeAxisTrack()

outdir <- "/rafalab/jill/MCC/ANALYSIS/NIH/plots/"
outdir <- "/arc/project/st-kdkortha-1/MCC/ANALYSIS/NIH/plots/"
dir.create(outdir, showWarnings = FALSE)

# fetch gene annot
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# get gene list
datdir <- "/rafalab/jill/MCC/ANALYSIS/"
datdir <- "/arc/project/st-kdkortha-1/MCC/ANALYSIS/"
genelist <- read.table(file.path(datdir, "ATAC-seq Gene Set for Class I APM - Sheet1.tsv"), sep="\t", stringsAsFactors=FALSE)$V1

tss <- c("ENST00000376809.10", "ENST00000412585.7", "ENST00000376228.9", "ENST00000376630.5")
names(tss) <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E")


bgdir <- "/rafalab/jill/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph"
bgdir <- "/arc/project/st-kdkortha-1/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph"

files.pval <- list.files(path = bgdir, pattern = "pval", 
	full.names = TRUE)

files.fc <- list.files(path = bgdir, pattern = "fc", 
	full.names = TRUE)


grl <- vector("list", length(files.fc))
for(f in seq_len(length(files.fc))){
	grl[[f]] <- readGeneric(files.fc[f], header = TRUE, 
		keep.all.metadata = TRUE)
}

idmap <- data.frame(
    nih = paste0("0", 1:9),
    mcc = c("275", "282", "290", "301", "320",
    	    "336", "350", "367", "2314"),
    stringsAsFactors=FALSE
	)


dtrk <- vector("list", length(tss))

names(dtrk) <- names(tss)
i <- 0
axisTrack <- GenomeAxisTrack()
tssTrack <- GRanges(seqnames = "chr6", 
              ranges = IRanges(start = 3000, width = 1))
tssTrack <- AnnotationTrack(tssTrack, name="TSS")

for (gene in names(tss)){

	i <- i + 1

	fetchgene <- select(org.Hs.eg.db, keys=gene, 
		keytype="SYMBOL", columns="ENTREZID")


    if(length(fetchgene) > 0){

    # pull out promoter for tx id'd by Patrick
    txtab <- transcripts(txdb, filter=list(gene_id = fetchgene$ENTREZID))
    txtab <- txtab[txtab$tx_name == tss[gene]]
    prom <- promoters(txtab, upstream = 3000, downstream=3000)
    grl_gene <- vector("list", length(grl))
    for (g in 1:length(grl_gene)){
      grl_gene[[g]] <- subsetByOverlaps(grl[[g]], prom)
      start(grl_gene[[g]][1,]) <- start(prom)
      end(grl_gene[[g]][length(grl_gene[[g]])]) <- end(prom)
    }

    # collapse grl_gene - get bp level data & combine
    grl_gene <- GRangesList(grl_gene)
    tiles <- unlist(tile(range(unlist(grl_gene)),
     		n=width(range(unlist(grl_gene)))))

    for (f in 1:length(grl_gene)){
        ol <- findOverlaps(grl_gene[[f]], tiles)
        mcols(tiles) <- cbind(mcols(tiles), fc = rep(NA, length(tiles)))
        names(mcols(tiles))[f] <- paste0(names(mcols(tiles))[f], "_", f)
        mcols(tiles)[ol@to,f] <- grl_gene[[f]]$X[ol@from]
    }

    ymax <- max(rowMeans(data.frame(mcols(tiles))))
    st <- start(tiles)
    start(tiles) <- start(tiles)-min(st)
    end(tiles) <- end(tiles)-min(st)
    dtrk[gene] <- DataTrack(tiles, name=gene)

    if (i < length(tss)){
      plotTracks(dtrk[gene], 
      	groups=rep(1,length(grl)), 
      	type=c("a"), ylim=c(0,9))
    }else{
      plotTracks(list.append(dtrk[gene], tssTrack, axisTrack), 
      	groups=rep(1,length(grl)), 
      	  type=c("a"), ylim=c(0,9))
    }
  }

}

pdf(file.path(outdir, paste0("ATAC_fc_average_tss.pdf")), 
		height=5, width=6)
plotTracks(list.append(dtrk, tssTrack, axisTrack), 
      	groups=rep(1,length(grl)), 
      	  type=c("a"), ylim=c(0,9))
dev.off()

#######################################################################
#######################################################################
#######################################################################
#######################################################################

# plot full length transcripts

for (gene in genelist){

	fetchgene <- select(org.Hs.eg.db, keys=gene, 
		keytype="SYMBOL", columns="ENTREZID")


    if(length(fetchgene) > 0){

    # pull out promoter for tx id'd by Patrick
    txtab <- transcripts(txdb, filter=list(gene_id = fetchgene$ENTREZID))
    txtab <- txtab[txtab$tx_name == tss[gene]]
  
    # pull out gene annotation from txdb
	genetab <- exons(txdb, filter=list(gene_id = fetchgene$ENTREZID))
	genetab$gene = fetchgene$ENTREZID
    genetab$symbol = gene
	chr <- as.character(unique(seqnames(genetab)))
	if (max(nchar(chr))>5){
      chr <- chr[nchar(chr) <= 6]
      genetab <- genetab[seqnames(genetab) %in% chr,]
    }

    # remove tx that span more than 10 times the next longest transcript
    test <- subsetByOverlaps(exonsBy(txdb), genetab)
    spans <- unlist(width(range(test)))
    if (max(spans) > 10*max(spans[-which.max(spans)])){
    	rmv <- test[[which.max(spans)]]$exon_id
    	if (sum(genetab$exon_id %in% rmv)>0){
    	  genetab <- genetab[-(genetab$exon_id %in% rmv),]
    	}
    }

    beg <- min(start(genetab), end(genetab))
    end <- max(start(genetab), end(genetab))
    w <- end-beg
    beg <- beg - 0.3*w
    end <- end + 0.3*w

    thistx <- transcripts(txdb, filter=list(gene_id = fetchgene$ENTREZID))$tx_name
	
    grtrack <- GeneRegionTrack(txdb, chromosome = chr,
    	feature=thistx,
    	transcriptAnnotation = "symbol", name = gene)
    grtrack <- subset(grtrack, from=beg, to=end)
    isthistx <- transcript(grtrack) %in% thistx 
    grtrack_target <- grtrack[isthistx,]
    grtrack_others <- grtrack[!isthistx,]
    names(grtrack_others) <- "Others"

    # get chr drawing 
	itrack <- IdeogramTrack(genome=gen, chromosome=chr)

    # get ATAC pval signal for plotting
	dtrackList.pval <- tab_all <- vector("list", length(files.pval))
	ymax <- 0

	for (f in seq_along(files.pval)){
		tab <- fread(files.pval[f])
		colnames(tab) <- c("chr", "start", "end", "pval")

		samp <- gsub("/rafalab/jill/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph/1623XD-53-", "", files.pval[f])

		samp <- gsub("/arc/project/st-kdkortha-1/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph/1623XD-53-", "", files.pval[f])
		samp <- substr(samp, 1, 2)
		samp <- idmap$mcc[match(samp, idmap$nih)]
	    keep <- (((tab$start >= beg) & (tab$end <= end)) |
            ((tab$start <= end) & (tab$end >= end)) |
            ((tab$start <= beg) & (tab$end >= beg)) |
            ((tab$start <= beg) & (tab$end >= end))) & 
            (tab$chr %in% chr)
		tab <- tab[keep,] 
        dtrackList.pval[[f]] <- DataTrack(data=tab$pval, 
			start=tab$start, 
			end=tab$end, 
			chromosome=chr, 
			genome=gen, 
			name=samp,
			type="h")
       

		ymax <- max(ymax, tab$pval)
	}
	
	pdf(file.path(outdir, paste0("ATAC_pval_", gene, ".pdf")), 
		height=6, width=6)
	plotTracks(list.append(itrack, grtrack_others, grtrack_target, dtrackList.pval),
		from=beg, to=end, col=NULL, ylim=c(0,ymax+1), type="h")
    dev.off()

    
     # get ATAC fc signal for plotting
	dtrackList.fc <-  tab_all <- vector("list", length(files.fc))
    ymax <- 0 

	for (f in seq_along(files.fc)){
		tab <- fread(files.fc[f])
		colnames(tab) <- c("chr", "start", "end", "fc")

		samp <- gsub("/rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph/1623XD-53-", "", files.pval[f])
		samp <- substr(samp, 1, 2)
		samp <- idmap$mcc[match(samp, idmap$nih)]
		keep <- (((tab$start >= beg) & (tab$end <= end)) |
            ((tab$start <= end) & (tab$end >= end)) |
            ((tab$start <= beg) & (tab$end >= beg)) |
            ((tab$start <= beg) & (tab$end >= end))) & 
            (tab$chr %in% chr)
		tab <- tab[keep,]

		dtrackList.fc[[f]] <- DataTrack(data=tab$fc, 
			start=tab$start, 
			end=tab$end, 
			chromosome=chr, 
			genome=gen, 
			name=samp,
			type="h")

		ymax <- max(ymax, tab$fc)

	}
	
	pdf(file.path(outdir, paste0("ATAC_foldchange_", gene, ".pdf")), 
		height=6, width=6)
	plotTracks(list.append(itrack, grtrack_others, grtrack_target, dtrackList.fc),
		from=beg, to=end, col=NULL, ylim=c(0,ymax + 1), type="h",
		chromosome=chr)
    dev.off()


  }
}








