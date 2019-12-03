library(Gviz)
library(data.table)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rlist)
library(AnnotationDbi)

gen <- "hg38"
gtrack <- GenomeAxisTrack()

outdir <- "/rafalab/keegan/MCC/ANALYSIS/NIH/plots/"
dir.create(outdir, showWarnings = FALSE)

# fetch gene annot
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# get gene list
datdir <- "/rafalab/keegan/MCC/ANALYSIS/"
genelist <- read.table(file.path(datdir, "ATAC-seq Gene Set for Class I APM - Sheet1.tsv"), sep="\t", stringsAsFactors=FALSE)$V1

bgdir <- "/rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph"

files.pval <- list.files(path = bgdir, pattern = "pval", 
	full.names = TRUE)

files.fc <- list.files(path = bgdir, pattern = "fc", 
	full.names = TRUE)


idmap <- data.frame(
    nih = paste0("0", 1:9),
    mcc = c("275", "282", "290", "301", "320",
    	    "336", "350", "367", "2314"),
    stringsAsFactors=FALSE
	)

for (gene in genelist){

	fetchgene <- select(org.Hs.eg.db, keys=gene, 
		keytype="SYMBOL", columns="ENTREZID")

    if(length(fetchgene) > 0){

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

    grtrack <- GeneRegionTrack(txdb, chromosome = chr,
    	transcriptAnnotation = "symbol", name = gene)
  
    # get chr drawing 
	itrack <- IdeogramTrack(genome=gen, chromosome=chr)

    # get ATAC pval signal for plotting
	dtrackList.pval <- vector("list", length(files.pval))
	ymax <- 0

	for (f in seq_along(files.pval)){
		tab <- fread(files.pval[f])
		colnames(tab) <- c("chr", "start", "end", "pval")

		samp <- gsub("/rafalab/keegan/MCC/DATA/NIH/atac_out/signal/macs2/rep1/bedgraph/1623XD-53-", "", files.pval[f])
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
	plotTracks(list.append(itrack, gtrack, grtrack, dtrackList.pval),
		from=beg, to=end, col=NULL, ylim=c(0,ymax+1), type="h")
    dev.off()

     # get ATAC fc signal for plotting
	dtrackList.fc <- vector("list", length(files.fc))
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
	plotTracks(list.append(itrack, gtrack, grtrack, dtrackList.fc),
		from=beg, to=end, col=NULL, ylim=c(0,ymax + 1), type="h",
		chromosome=chr)
    dev.off()

  }
}


