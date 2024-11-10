# SETUP ####

## Packages
library(readxl)
library(tidyverse)

## Functions ####
source("./R/functions.R")

## Metadata ####
dat <- read_xlsx("./data/MIP_Data.xlsx",sheet = "Colonization")
meta <- read_xlsx("./data/Sequencing_Palouse2021_metadata.xlsx")
meta$`Sample ID`

# remove missing samples (no seq file associated)



## Variables ####
dl_path_18S <- "./data/seq/18S"

accessions_18S <- dat$SRA_Accession_18S
filenames_f8S <- file.path(str_remove(dat$Fwd_18S,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")

# Download from SRA ####

## 18S ####
for(i in seq_along(accessions_18S)){
  system2("fasterq-dump",args = accessions_18S[i])
  system2("gzip", args = c(dl_filenames_1[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./data/ITS/",fung$filename[i]))
}


