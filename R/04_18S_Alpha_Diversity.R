# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(lmerTest)
library(corncob)
library(patchwork)

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

ps %>% 
  plot_richness(color='invasion',measures = "Observed",x = 'invasion') +
  geom_boxplot(alpha=.5) +
  facet_wrap(~site,scales = 'free_x') +
  labs(y="ASV richness",x="",color="Invasion\nstatus") +
  scale_color_manual(values = pal$pal.earthtones) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size=10),
        axis.title.y = element_text(face='bold',size=14))
ggsave("./output/figs/18S_ASV_Richness_Boxplot_by_site-invasion.png")



# look at just transition plots (and plant types in them)
ps %>% 
  subset_samples(invasion == "Transition") %>% 
  plot_richness(color='plant',measures = "Observed",x = 'invasion') +
  geom_boxplot(alpha=.5) +
  facet_wrap(~site,scales = 'free_x') +
  labs(y="ASV richness",x="",color="Plant") +
  scale_color_manual(values = pal$pal.earthtones) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face='bold',size=10),
        axis.title.y = element_text(face='bold',size=14))
ggsave("./output/figs/18S_ASV_Richness_Boxplot_Transition_by_plant.png")

# TAX GLOM ####
ps_spp <- 
  ps %>% 
  tax_glom("Species",NArm = FALSE,
           bad_empty = c("unclassified","", " ", "\t","sp.",NA))
# change NA to 'unclassified' at genus level
ps_spp@tax_table[,6][is.na(ps_spp@tax_table[,6])] <- "unclassified"

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
  labs(fill="Invasion\nstatus",y="ASV richness") 
ggsave("./output/figs/18S_ASV_Richness_Boxplot.png")
dat %>% 
  ggplot(aes(x=crust,y=species_richness,fill=invasion)) +
  geom_boxplot() +
  facet_wrap(~site,scales = 'free') +
  scale_fill_manual(values=pal$pal.earthtones) +
  labs(fill="Invasion\nstatus",y="Species richness") 
ggsave("./output/figs/18S_Species_Richness_Boxplot.png")

# boxplot to match Rae's figure 2

dat %>% 
  ggplot(aes(x=invasion,y=species_richness,color=crust,fill=crust)) +
  geom_boxplot(alpha=.5) +
  # facet_wrap(~invasion,scales = 'free_x') +
  labs(y="Species richness",x="",color="Biocrust",fill="Biocrust") +
  scale_color_manual(values = pal$pal.okabe) +
  scale_fill_manual(values = pal$pal.okabe) +
  theme(axis.text.x = element_text(face='bold',angle=0,size=14,hjust=.5),
        axis.text.y = element_text(face='bold',size=10),
        axis.title.y = element_text(face='bold',size=14),
        plot.title = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(face='bold',size=12),
        legend.title = element_text(face='bold',size=14))
ggsave("./output/figs/18S_Species_Richness_Boxplot_by_Crust_and_Invasion.png",height = 8,width = 8)


m <- dat %>% 
  lmer(data=.,formula=species_richness ~ invasion * crust + (1|site))
equatiomatic::extract_eq(m,wrap = TRUE)

# BARPLOTS ####
ps_spp %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Genus") + 
  facet_wrap(~site*invasion,scales = 'free') +
  labs(y="Relative abundance") +
  theme(legend.text = element_text(face='bold.italic',size=12)) +
  scale_fill_viridis_d(option = 'turbo')
ggsave("./output/figs/18S_Barplot_by_site-invasion.png",height = 8,width = 12)

spp_plot <- 
  ps_spp %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Species") + 
  facet_wrap(~site*invasion,scales = 'free') +
  labs(y="Relative abundance") +
  theme(legend.text = element_text(face='bold.italic',size=12)) +
  scale_fill_viridis_d(option = 'turbo')

spp_legend <- ggpubr::get_legend(spp_plot,position = 'right') %>% ggpubr::as_ggplot()

p1 <- spp_plot + theme(legend.position = 'none')

p1 / spp_legend
ggsave("./output/figs/18S_Barplot_by_site-invasion_spp-level.png",height = 14, width = 22)
p1
ggsave("./output/figs/18S_Barplot_by_site-invasion_spp-level_part1.png",height = 10, width = 18)
spp_legend
ggsave("./output/figs/18S_Barplot_by_site-invasion_spp-level_part2.png",height = 8, width = 22)
