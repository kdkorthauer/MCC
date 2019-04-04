library(readxl)
library(tidyverse)

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
  column_to_rownames(var = "sampleID")

  gather("snp", "genotype", -sampleID)

rownames(mat) <- fp$'Collaborator Sample ID'

