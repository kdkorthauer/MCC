library(PureCN)
library(genomation)

annot.dir <- "../../../PREPROCESS/DNA/annotation"
out.dir <- "../../../PREPROCESS/DNA/purecn"
dir.create(out.dir, showWarnings = FALSE)
bam.dir <- "../../../DATA/DNA"
bam.files <- list.files(path = bam.dir, pattern = "01.bam")
t.bam.files <- bam.files[grepl("CL", bam.files)]
c.bam.files <- bam.files[grepl("T", bam.files)]
n.bam.files <- bam.files[grepl("N", bam.files)]

## ----preprocess:gc-------------------------------------------------------------

reference.file <- file.path(annot.dir, "GATK_bundle_b37",
	                        "human_g1k_v37.fasta")

bed.file <- file.path(annot.dir,
	                  "whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.baits.interval_list_noheader")

# get mappability file: 
mappability.file <- file.path(annot.dir,
	                          "wgEncodeCrgMapabilityAlign36mer.bigWig")

if(!file.exists(mappability.file)){
  download.file("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig",
	            destfile = mappability.file)
}

# convert to GRanges
intervals <- readBed(bed.file)

# remove extraneous columns and make sure to keep strand info
strand(intervals) <- mcols(intervals)$name
intervals <- intervals[,-c(1,2)]

mappability <- import(mappability.file)

# file to save intervals
interval.file <- file.path(out.dir, "intervals.txt")

# preprocess intervals
if (!file.exists(interval.file)){
  message("Creating file ", interval.file, "...")

  preprocessIntervals(intervals, reference.file, 
      mappability=mappability, 
      output.file = interval.file)

  message("File ", interval.file, " finished.")
}else{
  message("File ", interval.file, " already created.")
}

## ----coverage-------------------------------------------------------

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






