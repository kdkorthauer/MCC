library(PureCN) # use version 1.13.19 or higher
library(genomation)
library(BiocParallel)

annot.dir <- "../../../PREPROCESS/DNA/annotation"
out.dir <- "../../../PREPROCESS/DNA/purecn-mt2"
plot.dir <- file.path(out.dir, "plots")
mutect.dir <- "../../../PREPROCESS/DNA/mutect2-gatk4-results"
dir.create(out.dir, showWarnings = FALSE)
dir.create(plot.dir, showWarnings = FALSE)
bam.dir <- "../../../DATA/DNA"

register(MulticoreParam(6))

## ----preprocess:gc-------------------------------------------------------------

reference.file <- file.path(annot.dir, "GATK_bundle_b37",
	                        "human_g1k_v37.fasta")

bed.file <- file.path(annot.dir,
	                  "whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list_noheader")

# get mappability file: 
mappability.file <- file.path(annot.dir,
	                          "wgEncodeCrgMapabilityAlign100mer.bigWig")

if(!file.exists(mappability.file)){
  download.file("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig",
	              destfile = mappability.file)
}

# get replication timing file:
reptime.file <- file.path(annot.dir,
                          "wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig")

if(!file.exists(reptime.file)){
  download.file("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
                destfile = reptime.file)
}

# convert to GRanges
intervals <- readBed(bed.file)

# remove extraneous columns and make sure to keep strand info
strand(intervals) <- mcols(intervals)$name
intervals <- intervals[,-c(1,2)]

# remove MT baits (strongly recommended by PureCN author)
intervals <- intervals[seqnames(intervals) != "MT",]
intervals <- dropSeqlevels(intervals, "MT")

# file to save intervals
interval.file <- file.path(out.dir, "intervals.txt")

# preprocess intervals
if (!file.exists(interval.file)){
  message("Creating file ", interval.file, "...")

  mappability <- import(mappability.file)
  seqlevels(mappability) <- mapSeqlevels(seqlevels(mappability), "NCBI")

  #remove chrMT from mappability:
  mappability <- mappability[seqnames(mappability) != "MT",]
  mappability <- dropSeqlevels(mappability, "MT")

  reptime <- import(reptime.file)
  seqlevels(reptime) <- mapSeqlevels(seqlevels(reptime), "NCBI")

  preprocessIntervals(intervals, reference.file, 
      mappability = mappability, 
      reptiming = reptime,
      output.file = interval.file,
      off.target = TRUE,
      min.target.width = 100)

  rm(mappability)
  message("File ", interval.file, " finished.")
}else{
  message("File ", interval.file, " already created.")
}

## ----coverage-------------------------------------------------------

bam.files <- list.files(path = bam.dir, pattern = "01.bam")
t.bam.files <- bam.files[grepl("CL", bam.files)]
c.bam.files <- bam.files[grepl("T", bam.files)]
n.bam.files <- bam.files[grepl("N", bam.files)]

for (f in bam.files){
  this.out <- file.path(out.dir, 
 	                    paste0(gsub(".bam", "", f), "_coverage.txt"))	
  if (!file.exists(this.out)){	
    message("Creating file ", this.out, "...")

    calculateBamCoverageByInterval(bam.file = file.path(bam.dir, f), 
     interval.file = interval.file, 
     output.file = this.out)

    message("File ", this.out, " finished.")
  }else{
  	message("File ", this.out, " already created.")
  }
}

## ----library-specific gc-correction------------------------------------------

coverage.files <- list.files(path = out.dir, pattern = "_coverage.txt")

for (f in coverage.files){
  this.out <- file.path(out.dir, 
                        gsub("_coverage.txt", "_coverage_loessnorm.txt", f)) 
  if (!file.exists(this.out)){  
    message("Creating file ", this.out, "...")
    
    pdf(file.path(plot.dir,
                 gsub("_coverage.txt", "_gcbias.pdf", f)),
        width = 10, height = 9.5)
    correctCoverageBias(file.path(out.dir,f), interval.file, 
      output.file = this.out, plot.bias = TRUE)
    dev.off()

    message("File ", this.out, " finished.")
  }else{
    message("File ", this.out, " already created.")
  }
}

## ----normaldb--------------------------------------------------------------

#### uncomment following after loessnorm generated
coverage.files <- list.files(path = out.dir, pattern = "_loessnorm.txt")
n.coverage.files <- coverage.files[grepl("N", coverage.files)]
t.coverage.files <- coverage.files[grepl("T", coverage.files)]
cl.coverage.files <- coverage.files[grepl("CL", coverage.files)]

ndb.file <- file.path(out.dir, "normalDB.rds")
if (!file.exists(ndb.file)){
  message("Creating ", ndb.file, "...")
  normalDB <- createNormalDatabase(file.path(out.dir, n.coverage.files))
 
  # serialize, so that we need to do this only once for each assay
  saveRDS(normalDB, file = ndb.file)
  message("Finished and saved ", ndb.file)
}else{
  message(ndb.file, " already created. Loading.")
  normalDB <- readRDS(ndb.file)
}

## ----normaldbpca-----------------------------------------------------------

for (f in c(t.coverage.files, cl.coverage.files)){
  this.out <- file.path(out.dir, 
                        gsub("coverage_loessnorm.txt", 
                             "coverage_poolnorm.txt", f)) 
  if (!file.exists(this.out)){  
    message("Creating file ", this.out, "...")

    pool <- calculateTangentNormal(file.path(out.dir, f), normalDB,
                                   ignore.sex = TRUE)
    PureCN:::.writeCoverage(pool, output.file = this.out) # save to file

    message("File ", this.out, " finished.")
  }else{
    message("File ", this.out, " already created.")
  }
}


## ----calculatemappingbias--------------------------------------------------
# Skip - can't run mutect with sample of ~40 PON (panel of normals)


## ----expected variance - target weight file-----------------------------------
interval.weight.file <- file.path(out.dir, "interval_weights.txt")
if (!file.exists(interval.weight.file)){
  pdf(file.path(plot.dir, "interval_weight_plot.pdf"))
    calculateIntervalWeights(normalDB$normal.coverage.files, 
                         interval.weight.file,
                         plot = TRUE)
  dev.off()
}

## ----ucsc_segmental--------------------------------------------------------

# find low quality regions
snp.blacklist.file <- file.path(annot.dir,
                                "hg19_simpleRepeats.bed")

if (!file.exists(snp.blacklist.file)) {
    library(rtracklayer)
    mySession <- browserSession("UCSC")
    genome(mySession) <- "hg19"
    simpleRepeats <- track(ucscTableQuery(mySession, 
        track="Simple Repeats", table="simpleRepeat"))
    export(simpleRepeats, snp.blacklist.file)
}

vcf.files <- list.files(path = mutect.dir, 
                        pattern = "01.vcf")

for (f in vcf.files){
  this.out <- file.path(out.dir, 
                        gsub(".vcf", 
                             "_absolute.rds", f)) 

  mutect.stats.file <- file.path(mutect.dir,                         
                                 gsub(".vcf", 
                                      "_call_stats.txt", f))

  ncov <- file.path(out.dir, 
                    gsub(".vcf", "_coverage_poolnorm.txt", f)) 

  tcov <- file.path(out.dir, 
                    gsub(".vcf", "_coverage_loessnorm.txt", f)) 

  plot.file <- file.path(plot.dir, 
                    gsub(".vcf", "_plots.pdf", f))
  
  if (!file.exists(tcov)){
    message("Tumor coverage file missing for ", f, ". Skipping.")
    next;
  }

  if (!file.exists(ncov)){
    message("Normal coverage file missing for ", f, ". Using process-matched")
    ncov <- file.path(out.dir, 
                      gsub(".vcf", "_coverage_poolnorm.txt", 
                           "DFCI-5368-CL-01.vcf")) 
  }

 
  if (grepl("-C-|-CL-", this.out)){
    tp=seq(0.9, 0.99,by=0.01)
    message("Cell line sample; Restricting purity to [",
            min(tp), ",", max(tp), "]")
  }else{
    tp=seq(0.3, 0.98, by = 0.01)
    message("Tumor sample; Restricting purity to [",
            min(tp), ",", max(tp), "]")
  }

  if (!file.exists(this.out)){  
    message("Creating file ", this.out, "...")

    pdf(plot.file, width = 10, height = 7)
    ret <- runAbsoluteCN(normal.coverage.file = ncov, 
                         tumor.coverage.file = tcov, 
                         vcf.file = file.path(mutect.dir, f), 
                         genome = "hg19", 
                         sampleid = gsub(".vcf", "", f), 
                         interval.file = interval.file, 
                         normalDB = normalDB,
                         args.filterVcf = list(snp.blacklist = snp.blacklist.file, 
                                              stats.file = mutect.stats.file),
                         post.optimize = FALSE,
                         test.purity = tp) 
    dev.off()
    saveRDS(ret, file = this.out)

    pdf(file.path(plot.dir, gsub(".vcf", "_plots_overview.pdf", f)), 
      width = 6, height = 6)
    plotAbs(ret, type="overview")
    dev.off()

    pdf(file.path(plot.dir, gsub(".vcf", "_plots_hist.pdf", f)), 
      width = 6, height = 6)
    plotAbs(ret, 1, type="hist")
    dev.off()

    pdf(file.path(plot.dir, gsub(".vcf", "_plots_baf.pdf", f)), 
      width = 8, height = 8)
    plotAbs(ret, 1, type="BAF")
    dev.off()

    message("File ", this.out, " and plots finished.")
  }else{
    message("File ", this.out, " already created.")
  }
}

