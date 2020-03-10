
# script to download fastq files for control lines to compare MCC line RNAseq

library(SRAdb)
library(DBI)

srafile <- "/scratch/st-kdkortha-1/SRAmetadb.sqlite"
if (!file.exists(srafile)){
  srafile = getSRAdbFile(destdir = "/scratch/st-kdkortha-1")
}
con = dbConnect(RSQLite::SQLite(), srafile)


# fibroblasts
fib <- listSRAfile('SRP126422', con)
fib <- fib[fib$experiment %in% c("SRX3459554", "SRX3459555", "SRX3459556", "SRX3459557"),]


# keratinocytes
ker <- listSRAfile('SRP131347', con)
ker <- ker[ker$experiment %in%
    c("SRX3145870", "SRX3145871", "SRX3145872", "SRX3145873", "SRX3145874", "SRX3145875"),]

datdir <- "../../../DATA/RNA/control_lines"
dir.create(datdir)

# download fib
dir.create(file.path(datdir, "fibroblasts"))
getSRAfile(fib$run, destDir = file.path(datdir, "fibroblasts"), con, fileType='fastq')

# download keratin
dir.create(file.path(datdir, "keratinocytes"))
getSRAfile(sort(ker$run), destDir = file.path(datdir, "keratinocytes"), con, fileType='fastq')



# download files from google drive (patrick)
prefix <- "https://drive.google.com/uc?export=download&id="
prefix_esc <- "https://drive.google.com/uc\\?export=download&id="
mkl1 <- paste0(prefix, "MKL1-RNAseq-replicate", 1:3, ".fastq.gz")
waga <- paste0(prefix, "Waga-RNAseq-replicate", 1:3, ".fastq.gz")
dir.create(file.path(datdir, "mkl1"))
dir.create(file.path(datdir, "waga"))


download.file("https://www.dropbox.com/s/f60ewzyzvelo34t/MKL1-RNAseq-replicate1.fastq.gz?dl=1", 
	destfile = file.path(datdir, "mkl1", gsub(prefix_esc, "", mkl1[1])))
download.file("https://www.dropbox.com/s/8910gwpe4dhtyb2/MKL1-RNAseq-replicate2.fastq.gz?dl=1", 
	destfile = file.path(datdir, "mkl1", gsub(prefix_esc, "", mkl1[2])))
download.file("https://www.dropbox.com/s/6ww9u4pbinbipyt/MKL1-RNAseq-replicate3.fastq.gz?dl=1", 
	destfile = file.path(datdir, "mkl1", gsub(prefix_esc, "", mkl1[3])))

download.file("https://www.dropbox.com/s/k12ruwcfhlugfsl/Waga-RNAseq-replicate1.fastq.gz?dl=1", 
	destfile = file.path(datdir, "waga", gsub(prefix_esc, "", waga[1])))
download.file("https://www.dropbox.com/s/qeb1hfb6mcoavnq/Waga-RNAseq-replicate2.fastq.gz?dl=1", 
	destfile = file.path(datdir, "waga", gsub(prefix_esc, "", waga[2])))
download.file("https://www.dropbox.com/s/emhs8yh3is2shmh/Waga-RNAseq-replicate3.fastq.gz?dl=1", 
	destfile = file.path(datdir, "waga", gsub(prefix_esc, "", waga[3])))



