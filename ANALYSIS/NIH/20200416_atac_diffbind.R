library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(GenomicRanges)
library(biomaRt)

# construct samplesheet with file paths
datdir <- "/arc/project/st-kdkortha-1/MCC/DATA/NIH/atac_out"

id <- list.files(path=file.path(datdir,"peak/macs2/overlap"),
  pattern="*filt.narrowPeak.xls")
id <- gsub("1623XD-53-", "", id)
id <- substr(id, 1, 2)
resp <- rep("Repsonsive", 9)
names(resp) <- c("275", "282", "290", "301", "320",
          "336", "350", "367", "2314")
resp[names(resp) %in% c("350", "336")] <- "NonResponsive"

resp_rank <- c(NA, 6, 7, 5, 4, 2, 1, 8, 3)
resp_val <- c(NA, 7.263374486, 9, 6.605166052, 3.066914498, 
              1.100671141, 1.092592593, 12.38267148, 2.48496994)

viral <- rep("ViralPos", 9)
viral[names(resp) %in% c("350", "320", "282", "290")] <- "ViralNeg"


peakfiles <- list.files(path=file.path(datdir, "peak/macs2/overlap"),
  pattern="*filt.narrowPeak.xls", full.names = TRUE)

bamfiles <- list.files(path=file.path(datdir, "align/rep1"),
  pattern="*.nodup.bam$", full.names = TRUE)

sampsheet <- data.frame(SampleID=names(resp),
  Condition = factor(ifelse(resp_val<=4, "NonResponsive", "Responsive")),
  Treatment = factor(viral),
  bamReads = bamfiles,
  Peaks = peakfiles)

# exclude 275/277
sampsheet <- sampsheet[!sampsheet$SampleID == "275",]
sampsheet2 <- sampsheet[,-4]
sampsheet2$Peaks <- gsub(".xls", "", basename(as.character(sampsheet2$Peaks)))
write.table(sampsheet2, file = "plots/atac_metadata.txt", 
  sep="\t", quote=FALSE, row.names=FALSE)

if(!file.exists("plots/ATAC_DiffPeaks_Viral.rds") | 
   !file.exists("plots/ATAC_DiffPeaks_IFNgResp.rds") ){

  # read in peaksets
  ps <- dba(sampleSheet = sampsheet, peakCaller = "macs")

  # count reads in peak sets from bams
  ps <- dba.count(ps)

  # differential binding affinity analysis
  ps1 <- dba.contrast(ps, categories="Treatment", minMembers=2,
                      block = DBA_CONDITION)
  ps1 <- dba.analyze(ps1, bReduceObjects=FALSE)

  ps2 <- dba.contrast(ps, categories="Condition", minMembers=2,
                         block = DBA_TREATMENT)
  ps2 <- dba.analyze(ps2, bReduceObjects=FALSE)

  # plotting and reporting
  pdf("plots/ATAC_DiffPeaks_Viral_FDR_0.10.pdf")
  plot(ps1, contrast=1, ColAttributes="Treatment", th=0.1)
  dev.off()


  pdf("plots/ATAC_DiffPeaks_IFNgResp_FDR_0.10.pdf")
  plot(ps2, contrast=1, ColAttributes="Condition", th=0.1)
  dev.off()

  # save rds
  saveRDS(ps1, file = "plots/ATAC_DiffPeaks_Viral.rds")
  saveRDS(ps2, file = "plots/ATAC_DiffPeaks_IFNgResp.rds")
}

add_annot <- function(gr){
  
  # Give ranges numeric names in order
  names(gr) <- c(1:length(gr))
  
  # Create GRanges object with annotations from TxDb database
  annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature="gene")
  
  # Annotate granges with the nearest TSS
  annot <- annotatePeakInBatch(gr, 
                               AnnotationData=annoData, 
                               featureType = "TSS",
                               output="nearestLocation",
                               PeakLocForDistance = "start")
  
  # Load mart
  ensembl <- useMart(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl")
  
  # Add gene information
  annot <- addGeneIDs(annot, mart = ensembl, feature_id_type = "entrezgene_id",
                      IDs2Add = c("hgnc_symbol"))
  annot$ENTREZ_ID <- unlist(sapply(sapply(names(annot), function(x) strsplit(x, "\\.")),
                     function(x) x[2]))

  return(annot)
}

ps1.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_Viral.rds"), th=0.1))
write.table(ps1.DB, quote=FALSE, sep="\t", row.names=FALSE,
  file="plots/ATAC_DiffPeaks_Viral_FDR_0.10.txt")

ps1.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_Viral.rds"), th=0.05))
write.table(ps1.DB, quote=FALSE, sep="\t", row.names=FALSE,
  file="plots/ATAC_DiffPeaks_Viral_FDR_0.05.txt")

ps1.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_Viral.rds"), th=1))
write.table(ps1.DB, quote=FALSE, sep="\t", row.names=FALSE,
  file="plots/ATAC_DiffPeaks_Viral_All.txt")



ps2.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_IFNgResp.rds"), th=0.1))
write.table(ps2.DB, quote=FALSE, sep="\t", 
  file="plots/ATAC_DiffPeaks_IFNgResp_FDR_0.10.txt")

ps2.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_IFNgResp.rds"), th=0.05))
write.table(ps2.DB, quote=FALSE, sep="\t", 
  file="plots/ATAC_DiffPeaks_IFNgResp_FDR_0.05.txt")

ps2.DB <- add_annot(dba.report(readRDS("plots/ATAC_DiffPeaks_IFNgResp.rds"), th=1))
write.table(ps2.DB, quote=FALSE, sep="\t", 
  file="plots/ATAC_DiffPeaks_IFNgResp_All.txt")