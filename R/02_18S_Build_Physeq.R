# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(dada2)
library(decontam)
library(readxl)

# functions
source("./R/functions.R")

# metadata
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

# add cutadapt file paths
cutadapt_dir <- list.dirs(full.names = TRUE)[grep("cutadapt$",list.dirs(full.names = TRUE))]
filtN_dir <- list.dirs(full.names = TRUE)[grep("filtN$",list.dirs(full.names = TRUE))]
metadata$cutadapt_fwd_paths <- paste0(cutadapt_dir,"/",metadata$sample_id,"_cutadapt_fwd.fastq.gz")
metadata$cutadapt_rev_paths <- paste0(cutadapt_dir,"/",metadata$sample_id,"_cutadapt_rev.fastq.gz")

# subset metadata to samples clearly present in cutadapt
metadata <- metadata[file.exists(metadata$cutadapt_fwd_paths),]

# add control column
metadata$control <- FALSE
metadata$control[grep("BLANK|NEG",metadata$sample_id)] <- TRUE

# RUN DADA2 ####
build_asv_table(metadata=metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                run.id.colname = "run", # name of column in metadata indicating which sequencing run a sample comes from
                run.id = "1", # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
                amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                sampleid.colname = "sample_id", # column name in metadata containing unique sample identifier
                fwd.fp.colname = "cutadapt_fwd_paths", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
                rev.fp.colname = "cutadapt_rev_paths", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
                fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
                rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
                maxEE = c(2,2), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
                truncQ = 2, # special value denoting "end of good quality sequence"
                trim.right = 0, # trim extra off right end for really crappy seq runs?
                rm.phix = TRUE, # remove phiX sequences?
                compress = TRUE, # gzip compression of output?
                multithread = (parallel::detectCores() -1), # how many cores to use? Set to FALSE on windows
                single.end = FALSE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
                filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
                asv.table.dir = "./data/ASV_Tables", # path to directory where final ASV table will be saved
                random.seed = 666
)

# reload point
seqtab.nochim <- readRDS("./data/ASV_Tables/SSU_ASV_Table.RDS")


# ASSIGN TAXONOMY ####

# load path to database 
SSU_DB <- "./taxonomy/Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta.gz"

# set output file path
outfile <- "./data/ASV_Tables/SSU_Taxonomy_Table.RDS"

# assign taxonomy
tax <- assign_taxonomy_to_asv_table(asv.table=seqtab.nochim,
                                    tax.database=SSU_DB,
                                    multithread=(parallel::detectCores()-2),
                                    random.seed=666,
                                    try.rc = TRUE,
                                    min.boot=50)
# export as RDS
saveRDS(tax,outfile)

# BUILD PHYSEQ ####

# metadata
meta <- sample_data(metadata)
sample_names(meta) <- meta$sample_id

# remove neg controls
meta <- meta[grep("^NEG|^BLANK",sample_names(meta),invert = TRUE),]

# tax table
taxa <- tax_table(tax)

# asv table
asv <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
sample_names(asv)

if(identical(sample_names(asv),sample_names(meta))){
  ps <- phyloseq(meta,taxa,asv)
} else (cat("Check sample names match"))

ps
saveRDS(ps,"./data/physeq_18S_not_cleaned.RDS")
