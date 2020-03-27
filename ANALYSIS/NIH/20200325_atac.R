library(Gviz)
library(data.table)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rlist)
library(AnnotationDbi)
library(genomation)
library(Rsamtools)
library(readxl)

gen <- "hg38"
gtrack <- GenomeAxisTrack()

outdir <- "/rafalab/jill/MCC/ANALYSIS/NIH/plots/"
outdir <- "/arc/project/st-kdkortha-1/MCC/ANALYSIS/NIH/plots/"
dir.create(outdir, showWarnings = FALSE)

# fetch gene annot
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# get gene list
datdir <- "/arc/project/st-kdkortha-1/MCC/ANALYSIS/"
genelist_tab <- read_excel(file.path(datdir, 
    "ATAC-seq Gene Set for Class I APM.xlsx"))
genelist <- genelist_tab$`Gene Symbol`
gene_tx <- genelist_tab$`Transcript Name`

tss <- genelist_tab$`ENSEMBL Transcript ID (to use for TSS)`
names(tss) <- genelist

bamdir <- "/arc/project/st-kdkortha-1/MCC/DATA/NIH/atac_out/align/rep1"
files.bam <- list.files(path = bamdir, pattern = "nodup.bam", 
    full.names = TRUE)
files.bam <- files.bam[!grepl(".bai", files.bam)]


idmap <- data.frame(
    nih = paste0("0", 1:9),
    mcc = c("275", "282", "290", "301", "320",
    	    "336", "350", "367", "2314"),
    stringsAsFactors=FALSE
	)


#######################################################################
#######################################################################
#######################################################################
#######################################################################

# get library sizes
libsize <- rep(NA, length(files.bam))
for(f in seq_len(length(files.bam))){
      libsize[f] <- countBam(files.bam[f])$records
}

  # "HLA-A" "HLA-B" "HLA-C" "HLA-E" "HLA-F" 
maxy <- c(16,20,20,24,18,
  # "HLA-G" "B2M" "TAP1"  "TAP2" "TAPBP"
          10,80,55,40,36,
  # "PSMB8" "PSMB9" "NLRC5-201" "NLRC5-203" "NLRC5-217"    
          50,60,45,10,10,
  # "GAPDH" "ACTB" 
          100,60)



# plot full length transcripts
inc <- 0

for (gene in genelist){
    inc <- inc + 1 

    fetchgene <- select(org.Hs.eg.db, keys=gene, 
        keytype="SYMBOL", columns="ENTREZID")


    if(length(fetchgene) > 0){

    # pull out promoter for tx id'd by Patrick
    txtab <- transcripts(txdb, filter=list(gene_id = fetchgene$ENTREZID))
    txtab <- txtab[txtab$tx_name == tss[inc]]
  
    # pull out gene annotation from txdb
    genetab <- exons(txdb, filter=list(tx_name = tss[inc]))
    genetab$gene = fetchgene$ENTREZID
    genetab$symbol = gene
    chr <- as.character(unique(seqnames(genetab)))
    if (max(nchar(chr))>5){
      chr <- chr[nchar(chr) <= 6]
      genetab <- genetab[seqnames(genetab) %in% chr,]
    }


    spn <- promoters(txtab, upstream = 3000, downstream=3000)
    spn$gene = fetchgene$ENTREZID
    spn$symbol = gene
    beg <- min(start(spn), end(spn))
    end <- max(start(spn), end(spn))

    grtrack <- GeneRegionTrack(genetab, chromosome = chr, 
        name = gene_tx[inc], genome=gen,
        transcript=genetab$gene,
        rotation.title=0,
        cex.title=0.33)

    # get chr drawing 
    itrackfile <- file.path(paste0("/scratch/st-kdkortha-1/MCC/tmp/itrack_", 
        chr, ".rds"))
    if (!file.exists(itrackfile)){
      itrack <- IdeogramTrack(genome=gen, chromosome=chr, 
        name=gsub("chr", "", chr))
      saveRDS(itrack, itrackfile)
    }else{
      itrack <- readRDS(itrackfile)
    }

    # get ATAC pval signal for plotting
    bTrack <- vector("list", length(files.bam))

    for(f in seq_len(length(files.bam))){
      samp <- gsub("/arc/project/st-kdkortha-1/MCC/DATA/NIH/atac_out/align/rep1/1623XD-53-", "", files.bam[f])
      samp <- substr(samp, 1, 2)
      samp <- idmap$mcc[match(samp, idmap$nih)]

      bTrack[[f]] <- AlignmentsTrack(files.bam[f], isPaired = FALSE,
                     genome = "hg38",
                     chromosome  = chr,
                     transformation=function(x) {x/libsize[f]*median(libsize)},
                     name=samp)
    }

    
    pdf(file.path(outdir, paste0("ATAC_coverage_", gene_tx[inc], ".pdf")), 
        height=5, width=2)
    plotTracks(list.append(itrack, grtrack, bTrack),
        from=beg, to=end, ylim=c(0,maxy[inc]), 
        type="coverage",
        sizes=c(0.6, 0.3, rep(1, length(bTrack))))
    dev.off()
  }
}




