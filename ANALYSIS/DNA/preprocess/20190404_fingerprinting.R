library(readxl)
library(tidyverse)
library(colordistance)

datdir <- file.path("/n/irizarryfs01_backed_up/kkorthauer",
	"MCC/DATA/DNA")

fp1 <- read_xls(file.path(datdir, "20190403_fingerprinting_candace_1.xls"))
fp2 <- read_xls(file.path(datdir, "20190403_fingerprinting_candace_2.xls"))

# can we merge?
identical(colnames(fp1), colnames(fp2))

fp <- rbind(fp1, fp2)

# which is 5370 
ix <- which(grepl("5370", fp$'Collaborator Sample ID'))
fp[ix,]

# fingerprinting cols 
fix <- which(grepl("rs", colnames(fp)))
id <- which(grepl("Collaborator Sample ID", colnames(fp)))

fp <- fp[,c(id,fix)] %>%
  rename( sampleID = 'Collaborator Sample ID') %>%
  filter(grepl("DFCI", sampleID)) %>%
  unique() %>%
  column_to_rownames(var = "sampleID") %>%
  #replace(., .=="--", NA) %>%
  unite(genotype, sep="")

sdist <- adist(fp$genotype, fp$genotype)
rownames(sdist) <- colnames(sdist) <- rownames(fp)

ident <- function(x){
	as.dist(x)
}

pdf(file.path("/n/irizarryfs01_backed_up/kkorthauer",
	          "MCC/ANALYSIS/DNA/plots", "heatmap_fingerprinting_mcc.pdf"))
heatmap.2(sdist, margins = c(9, 9), symm = TRUE, trace="none", 
	      distfun = ident)
dev.off()

# repeat for all samps
fp <- rbind(fp1, fp2)

# which is 5370 
ix <- which(grepl("5370", fp$'Collaborator Sample ID'))
fp[ix,]

# fingerprinting cols 
fix <- which(grepl("rs", colnames(fp)))
id <- which(grepl("Collaborator Sample ID", colnames(fp)))

fp <- fp[,c(id,fix)] %>%
  rename( sampleID = 'Collaborator Sample ID') %>%
 #filter(grepl("DFCI", sampleID)) %>%
  unique() %>%
  column_to_rownames(var = "sampleID") %>%
  #replace(., .=="--", NA) %>%
  unite(genotype, sep="")

sdist <- adist(fp$genotype, fp$genotype)
rownames(sdist) <- colnames(sdist) <- rownames(fp)

ident <- function(x){
	as.dist(x)
}

pdf(file.path("/n/irizarryfs01_backed_up/kkorthauer",
	          "MCC/ANALYSIS/DNA/plots", "heatmap_fingerprinting_all.pdf"),
   width = 9, height = 9)
heatmap.2(sdist, margins = c(10, 8), symm = TRUE, trace="none", 
	      distfun = ident)
dev.off()



#  gather("snp", "genotype", -sampleID)
