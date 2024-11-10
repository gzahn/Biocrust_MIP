# SETUP ####

# install mycobank package, if not already installed
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}
  
if(!"mycobank" %in% installed.packages()[,"Package"]){
  devtools::install_github("gzahn/mycobank")
}


if(!"ShortRead" %in% installed.packages()[,"Package"]){
  BiocManager::install("ShortRead")
}

# packages
library(tidyverse)
library(mycobank)
library(ShortRead)

# paths
maarjam_ssu <- "./maarjam_database_SSU_2021.fasta"
maarjam_ssu_type <- "./maarjam_database_SSU_TYPE_2021.fasta"
maarjam_its <- "./maarjam_database_full_ITS_2021.fasta"

# load fasta files
ssu <- readFasta(maarjam_ssu)
ssu_type <- readFasta(maarjam_ssu_type)
its <- readFasta(maarjam_its)

# get fasta headers
ssu_headers <- ssu@id %>% as.character()
ssu_type_headers <- ssu_type@id %>% as.character()
its_headers <- its@id %>% as.character()

# load mycobank database
db <- mycobank::get_mycobank_db()

# ADD TAXONOMY ####

## Pull genera names ####
ssu_genera <- ssu_headers %>% str_split(" ") %>% map_chr(3)
ssu_type_genera <- ssu_type_headers %>% str_split(" ") %>% map_chr(3)
its_genera <- its_headers %>% str_split(" ") %>% map_chr(3)

## Pull VTX names ####
ssu_vtx <- ssu_headers %>% str_extract("VTX.*")
ssu_type_vtx <- ssu_type_headers %>% str_extract("VTX.*")
  its.pattern <- paste0("(",paste0(its_genera %>% unique %>% paste(collapse = "|"),")",".*"))
its_vtx <- its_headers %>% str_extract(its.pattern)

## Extract taxonomic ranks ####

# unique list of genera present in maarjAM
unique_genera <- c(ssu_genera,ssu_type_genera,its_genera) %>% unique

# mycobank taxonomy retrieval (for all higher ranks)
genera_taxonomy_list <- mycobank::get_mycobank_taxonomy(taxa = unique_genera,db)

## Build lookup table ####
genera_taxonomy_list <- data.frame(genus=names(genera_taxonomy_list),taxonomy=genera_taxonomy_list)

# split string into separate columns
lookup_df <- 
genera_taxonomy_list %>% 
  separate(taxonomy, into = c("Kingdom","Extra1","Phylum","Extra2","Class","Order","Family")) %>% 
  dplyr::select(-starts_with("Extra"))

# create Eukaryome-style format column
# >k__cf.Fungi;p__unclassified;c__unclassified;o__unclassified;f__unclassified;g__unclassified;s__unclassified
lookup_df <- 
  lookup_df %>% 
  mutate(taxonomy = paste0("k__",Kingdom,";",
                           "p__",Phylum,";",
                           "c__",Class,";",
                           "o__",Order,";",
                           "f__",Family,";"))

## Build new headers ####
ssu_new_headers <- 
data.frame(genus=ssu_genera) %>% 
  left_join(lookup_df) %>% 
  mutate(new_taxonomy=paste0(taxonomy,"g__",genus,";",
                             "s__",ssu_vtx)) %>% 
  pluck('new_taxonomy')

ssu_type_new_headers <- 
  data.frame(genus=ssu_type_genera) %>% 
  left_join(lookup_df) %>% 
  mutate(new_taxonomy=paste0(taxonomy,"g__",genus,";",
                             "s__",ssu_type_vtx)) %>% 
  pluck('new_taxonomy')

its_new_headers <- 
  data.frame(genus=its_genera) %>% 
  left_join(lookup_df) %>% 
  mutate(new_taxonomy=paste0(taxonomy,"g__",genus,";",
                             "s__",its_vtx)) %>% 
  pluck('new_taxonomy')

# ADD NEW HEADERS ####
ssu@id <- BStringSet(ssu_new_headers)
ssu_type@id <- BStringSet(ssu_type_new_headers)
its@id <- BStringSet(its_new_headers)

# EXPORT UPDATED FASTAS ####
writeFasta(ssu,"./maarjam_database_SSU_2021_newheaders.fasta",width=10000)
writeFasta(ssu_type,"./maarjam_database_SSU_TYPE_2021_newheaders.fasta",width=10000)
writeFasta(its,"./maarjam_database_full_ITS_2021_newheaders.fasta",width=10000)

