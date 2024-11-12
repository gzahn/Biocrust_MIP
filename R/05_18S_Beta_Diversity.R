# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(lmerTest)
library(corncob)
library(patchwork)
library(broom)

# functions
source("./R/functions.R")

# theme
theme_set(theme_minimal() +
            theme(strip.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=14),
                  axis.title = element_text(face='bold',size=14)
            )
          )

# data
ps <- readRDS("./data/physeq_18S_clean.RDS")
# "OTU Table matrix"
mat <- 
  ps %>% 
  transform_sample_counts(ra) %>% 
  otu_table() %>% 
  as('matrix')

# ORDINATIONS ####  
ord <- 
  ps %>% 
  transform_sample_counts(ra) %>% 
  ordinate(method = "NMDS",distance = 'bray')
plot_ordination(ps,ord,color='site') +
  stat_ellipse()
plot_ordination(ps,ord,color='invasion') +
  stat_ellipse()
plot_ordination(ps,ord,color='crust') +
  stat_ellipse()
plot_ordination(ps,ord,color='invasion') +
  stat_ellipse() +
  facet_wrap(~site,scales = 'free') +
  theme()
data.frame(site=ps@sam_data$site,
           invasion=ps@sam_data$invasion,
           crust=ps@sam_data$crust,
           MDS1=ord$points[,1],
           MDS2=ord$points[,2]) %>% 
  ggplot(aes(x=MDS1,y=MDS2,color=invasion)) +
  geom_point(size=3) +
  stat_ellipse(linetype=2) + 
  facet_wrap(~site,scales = 'free') +
  scale_color_manual(values = pal$pal.earthtones) +
  labs(color="Invasion\nstatus")
ggsave("./output/figs/18S_NMDS_Plot_site_by_invasion.png", height = 8,width = 10)

data.frame(site=ps@sam_data$site,
           invasion=ps@sam_data$invasion,
           crust=ps@sam_data$crust,
           plant=ps@sam_data$plant,
           MDS1=ord$points[,1],
           MDS2=ord$points[,2]) %>%
  dplyr::filter(invasion == "Transition") %>% 
  ggplot(aes(x=MDS1,y=MDS2,color=plant)) +
  geom_point(size=3) +
  stat_ellipse(linetype=2) + 
  facet_wrap(~site,scales = 'free') +
  scale_color_manual(values = pal$pal.okabe) +
  labs(color="Plant",title = "Transition plots only")

# PermANOVA ####
perm.mod <- 
  adonis2(mat ~ ps@sam_data$invasion * ps@sam_data$crust, strata = ps@sam_data$site) %>% 
  broom::tidy() %>% 
  mutate(term=term %>% str_remove_all("ps@sam_data\\$"))
perm.mod
perm.mod %>% 
  write_csv("./output/18S_PermANOVA_Table.csv")


# DIFFABUND ####

ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Species"))
# use raw count data for corncob
da_analysis_ASV <- differentialTest(formula = ~ invasion + invasion:site, #abundance
                                             phi.formula = ~ 1, #dispersion
                                             formula_null = ~ 1, #mean
                                             phi.formula_null = ~ 1,
                                             test = "Wald", boot = FALSE,
                                             data = ps,
                                             fdr_cutoff = 0.05,
                                             full_output = TRUE)
dat <- corncob:::plot.differentialTest(da_analysis_ASV,level = "Species",data_only = TRUE)
dat$taxa <- dat$taxa %>% str_remove(" \\(.*")
dat$variable <- dat$variable %>% str_remove_all("invasion|site") %>% 
  str_remove("\nDifferential Abundance") 
dat %>% 
  group_by(taxa,variable) %>% 
  summarize(x=mean(x,na.rm=TRUE),
            xmin=mean(xmin),
            xmax=mean(xmax)) %>% 
  mutate(effect=if_else(x>0,"Positive","Negative")) %>%
  mutate(effect= case_when(x>0 & xmin <= 0 ~ "None",
                           x>0 & xmin > 0 ~ "Positive",
                           x<0 & xmax < 0 ~ "Negative",
                           x<0 & xmax >= 0 ~ "None",)) %>% 
  mutate(variable=case_when(variable == "Invaded" ~ "Invaded:Kamiak",
                            variable == "Transition" ~ "Transition:Kamiak",
                            TRUE ~ variable)) %>% 
  dplyr::filter(grepl(":",variable)) %>% 
  separate(variable, into = c("invasion", "site")) %>% 
  mutate(invasion=invasion %>% factor(levels=c("Transition","Invaded"))) %>% 
  dplyr::filter(invasion != "Native") %>% 
  ggplot(aes(x=x,xmin=xmin,xmax=xmax,y=taxa,color=effect)) +
  geom_vline(xintercept = 0,linetype=21) +
  geom_errorbar(width=.25,linewidth=1) +
  facet_wrap(~site*invasion,ncol = 2,scales = 'free') +
  labs(x="Model estimate",y="Taxa",color="Effect of\n invasion",
       caption = "corncob model results showing differential abundance parameter estimates\nfor taxa with significant differences (P<0.05) from Native plots.") +
  scale_color_manual(values = c("darkred","gray","green4")) +
  theme(axis.text.y = element_text(face='bold.italic',size=7))
ggsave("./output/figs/18S_Sig-Diff_Taxa__by_site_Plot.png",height = 16, width = 12)




# pull model info out for reporting significant ASVs
mods <- da_analysis_ASV$significant_models
bbdml_mods <- map(mods,capture_mods)

# find the significant taxa
sig_taxa <- da_analysis_ASV$significant_taxa %>% otu_to_taxonomy(data=ps)

names(bbdml_mods) <- sig_taxa
joined_mods <- bbdml_mods %>% 
  purrr::map(as.data.frame) %>% 
  purrr::reduce(full_join)
joined_mods$taxon <- names(bbdml_mods)

# just look at sig. taxa in barplot
mergevar <- paste(ps_sig@sam_data$site,ps_sig@sam_data$invasion,sep="_")
ps@sam_data$mergevar <- mergevar
ps_sig <- 
  ps %>% 
  merge_samples("mergevar") %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(ps) %in% names(sig_taxa))
ps_sig@sam_data$site <- sample_names(ps_sig) %>% str_split("_") %>% map_chr(1)
ps_sig@sam_data$invasion <- sample_names(ps_sig) %>% str_split("_") %>% map_chr(2)

ps_sig <- 
  ps_sig %>% 
  subset_samples(sample_sums(ps_sig) > 0)
ps_sig %>% 
  plot_bar2(fill="Species") +
  facet_wrap(~invasion,scales = 'free_x')


sig_melt <- 
  ps_sig %>% 
  psmelt()
new.sample.order <- 
  sig_melt %>% 
  arrange(invasion,site,crust) %>% 
  pluck("Sample") %>% unique
sig_melt$Sample <- factor(sig_melt$Sample,levels = new.sample.order)

sig_melt$invasion <- sig_melt$invasion %>% factor(levels = c("Native","Transition","Invaded"))
sig_melt %>% 
  ggplot(aes(x=Abundance,y=Species)) +
  geom_boxplot(aes(color=invasion))
