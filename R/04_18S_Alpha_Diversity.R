# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(lmerTest)
library(corncob)

# functions
source("./R/functions.R")

# theme
theme_set(theme_minimal() +
            theme(strip.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=14)
            )
)

# data
ps <- readRDS("./data/physeq_18S_clean.RDS")


# RARECURVE ####
mat <- otu_table(ps) %>% as('matrix')
png("./output/figs/18S_rarefaction_curve.png")
rarecurve(mat,step = 100,label = FALSE)
dev.off()

# SIMPLE ALPHA DIV ####
ps@sam_data$asv_richness <- 
  ps %>% 
  estimate_richness(measures = c("Observed")) %>% 
  pluck("Observed")
plot_richness(ps,color='crust',measures = "Observed") +
  facet_wrap(~site,scales = 'free')

ps@tax_table[,7] %>% unname


# TAX GLOM ####
ps_spp <- 
  ps %>% 
  tax_glom("Species",NArm = FALSE,
           bad_empty = c("unclassified","", " ", "\t"))

# species richness
ps_spp@sam_data$species_richness <- 
  ps_spp %>% 
  estimate_richness(measures = "Observed") %>% 
  pluck("Observed")
# melt for plotting and modeling
dat <- 
  ps_spp %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  dplyr::rename("relative_abundance" = "Abundance")

# PLOT & MODEL "SPECIES" RICHNESS ####

# site, invasion, crust

dat %>% 
  ggplot(aes(x=crust,y=asv_richness,fill=invasion)) +
  geom_boxplot() +
  facet_wrap(~site,scales = 'free') +
  scale_fill_manual(values=pal$pal.earthtones) +
  labs(fill="Invasion\nstatus",y="Species richness") 
ggsave("./output/figs/18S_ASV_Richness_Boxplot.png")
dat %>% 
  ggplot(aes(x=crust,y=species_richness,fill=invasion)) +
  geom_boxplot() +
  facet_wrap(~site,scales = 'free') +
  scale_fill_manual(values=pal$pal.earthtones) +
  labs(fill="Invasion\nstatus",y="Species richness") 
ggsave("./output/figs/18S_Species_Richness_Boxplot.png")

dat %>% 
  lmer(data=.,formula=species_richness ~ invasion * crust + (1|site)) %>% 
  summary


# BARPLOTS ####
ps_spp %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Genus") + 
  facet_wrap(~site*invasion,scales = 'free') +
  labs(y="Relative abundance") +
  theme(legend.text = element_text(face='bold.italic',size=12)) +
  scale_fill_viridis_d(option = 'turbo')
ggsave("./output/figs/18S_Barplot_by_site-invasion.png",height = 8,width = 12)




