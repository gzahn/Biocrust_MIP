# SETUP ####

## Packages
library(readxl)
library(tidyverse)

## Functions ####
source("./R/functions.R")

## Metadata ####
dat <- read_xlsx("./data/full_metadata.xlsx")

## Variables ####
dl_path_18S <- "./data/seq/18S"
accessions_18S <- dat$sra_accession
fwd_filenames_18S <- file.path(dl_path_18S,dat$fwd_name)
rev_filenames_18S <- file.path(dl_path_18S,dat$rev_name)
dl_filenames_fwd <- paste0(accessions_18S,"_1.fastq")
dl_filenames_rev <- paste0(accessions_18S,"_2.fastq")

# Download from SRA ####
system("which ls")
system2("which",args = "ls")
## 18S ####
for(i in seq_along(accessions_18S)){
  system2(fasterq_dump,args = accessions_18S[i])
  system2("gzip", args = c(dl_filenames_fwd[i]))
  system2("gzip", args = c(dl_filenames_rev[i]))
  file.rename(paste0(dl_filenames_fwd[i],".gz"),fwd_filenames_18S[i])
  file.rename(paste0(dl_filenames_rev[i],".gz"),rev_filenames_18S[i])
}


