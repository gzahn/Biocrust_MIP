# SETUP ####

# packages
library(tidyverse)
library(phyloseq)

# functions
source("./R/functions.R")

# data
ps <- readRDS("./data/physeq_18S_not_cleaned.RDS")

# REMOVE NON-AMF
ps <- 
  ps %>% 
  subset_taxa(Phylum == "p__Glomeromycota")

# CHECK POS CONTROLS ####

# add pos control designations based on 'note' column 
ps@sam_data$control <- !is.na(ps@sam_data$note)

pos <- 
  ps %>% 
  transform_sample_counts(ra) %>% 
  subset_samples(control)
keepers <- taxa_sums(pos) > 0
names(keepers) <- NULL
pos %>% 
  subset_taxa(keepers) %>% 
  plot_bar2(fill = "Genus")
ggsave("./output/figs/18S_positive_control_barchart.png")


# remove positive controls
ps <- 
  ps %>% 
  subset_samples(control == FALSE)
# remove empty taxa 
ps <- 
  ps %>% 
  subset_taxa(taxa_sums(ps) > 0)
# check sample sums
all(ps %>% sample_sums() > 0)

# CLEAN ####

# clean taxa ranks
ps <- clean_ps_taxonomy(ps)

# remove useless columns
ps@sam_data$note <- NULL
ps@sam_data$run <- NULL
ps@sam_data$control <- NULL
ps@sam_data <- ps@sam_data[,grep("path",ps@sam_data %>% names,invert = TRUE)]

# convert factors
ps@sam_data$invasion <- factor(ps@sam_data$invasion,levels=c("Native","Transition","Invaded"))


# build better species-level names
ps <- 
  ps %>% 
  make_sane_taxa_names()

# add plant data
V <- ps@sam_data$sample_id %>% grepl(pattern="V$")
P <- ps@sam_data$sample_id %>% grepl(pattern="P$")
ps@sam_data <- 
  ps@sam_data %>% 
  as('data.frame') %>% 
  mutate(plant=case_when(V ~ "V",
                         P ~ "P")) %>% 
  sample_data()


# SAVE CLEAN PHYSEQ ####
saveRDS(ps,"./data/physeq_18S_clean.RDS")
