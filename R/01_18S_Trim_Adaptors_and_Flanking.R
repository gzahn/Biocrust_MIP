# SETUP ####

## Packages ####
library(tidyverse); packageVersion('tidyverse')
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel); packageVersion("parallel")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")

## Functions ####
source("./R/functions.R")

## Data ####

# sample metadata and filepaths
metadata <- read_xlsx("./data/Palouse2021_18S_sequencing_metadata.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::filter(sample != "Rhizosphere") %>%  # leave out rhizosphere data from this analysis
  mutate(fwd_18s_filepath = file.path("./data/seq/18S",fwd_18s_filepath), # add paths to file names
         rev_18s_filepath = file.path("./data/seq/18S",rev_18s_filepath)) %>% 
  dplyr::rename("crust" = "sample", # change column names to be more intuitive
                "invasion" = "treatment")

# add more metadata
metadata$amplicon <- "SSU"
metadata$run <- 1

## primer sequences ####

# SSU
WANDAf <- "CAGCCGCGGTAATTCCAGCT"  
AML2r <- "GAACCCAAACACTTTGGTTTCC"  

# RUN CUTADAPT ####

## On SSU Samples ####
remove_primers(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
               amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
               amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
               sampleid.colname = "sample_id", # column name in metadata containing unique sample identifier
               fwd.fp.colname = "fwd_18s_filepath", # name of column in metadata indicating fwd filepath to raw data
               rev.fp.colname = "rev_18s_filepath",
               fwd_pattern="_R1_",
               rev_pattern="_R2_",
               fwd_primer=WANDAf,
               rev_primer=AML2r,
               multithread=parallel::detectCores()-1)

