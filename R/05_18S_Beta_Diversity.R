# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(lmerTest)
library(corncob)
library(patchwork)
library(broom)
library(gdm)

# functions
source("./R/functions.R")

# variables
set.seed(666)

# theme
theme_set(theme_bw() +
            theme(strip.text = element_text(face='bold',size=14),
                  legend.title = element_text(face='bold',size=14),
                  axis.title = element_text(face='bold',size=14),
                  axis.text = element_text(face='bold',size=10)
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
  scale_color_manual(values=c("#6E016B","#9EBCDA","gray")) +
  labs(color="Invasion\nstatus") +
  theme(strip.background = element_blank(), legend.text = element_text(face='bold',size=12))

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
ggsave("./output/figs/18S_NMDS_Plot_transitionplots_by_plant.png",height = 8, width = 10)
# PermANOVA ####
perm.mod <- 
  adonis2(mat ~ ps@sam_data$invasion * ps@sam_data$crust, strata = ps@sam_data$site) %>% 
  broom::tidy() %>% 
  mutate(term=term %>% str_remove_all("ps@sam_data\\$"))
perm.mod
perm.mod %>% 
  write_csv("./output/18S_PermANOVA_Table.csv")


ps2 <- ps %>% subset_samples(invasion == "Transition") %>% transform_sample_counts(ra)
mat2 <- ps2 %>% otu_table() %>% as("matrix")
perm.mod2 <- 
  adonis2(mat2 ~ ps2@sam_data$plant, strata = ps2@sam_data$site) %>% 
  broom::tidy() %>% 
  mutate(term=term %>% str_remove_all("ps2@sam_data\\$"))
perm.mod2
perm.mod2 %>% 
  write_csv("./output/18S_PermANOVA_Table_plant.csv")

# GDM MODEL ####
# extract species by site info
ps_melt <- 
  ps %>% 
  tax_glom("Species",NArm = FALSE,
           bad_empty = c("unclassified","", " ", "\t","sp.",NA)) %>% 
  transform_sample_counts(ra) %>% 
  psmelt()
names(ps_melt)

# biological data
# get columns with xy, site ID, and species data
sppTab <- ps_melt %>% dplyr::select(OTU,Sample,longitude,latitude,site,invasion,crust,Abundance)
sppTab_native <- ps_melt %>% dplyr::select(OTU,Sample,longitude,latitude,site,invasion,crust,Abundance) %>% 
  dplyr::filter(invasion == "Native")
sppTab_transition <- ps_melt %>% dplyr::select(OTU,Sample,longitude,latitude,site,invasion,crust,Abundance) %>% 
  dplyr::filter(invasion == "Transition")
sppTab_invaded <- ps_melt %>% dplyr::select(OTU,Sample,longitude,latitude,site,invasion,crust,Abundance) %>% 
  dplyr::filter(invasion == "Invaded")


# get columns with site ID, env. data, and xy-coordinates
envTab <- ps_melt %>% dplyr::select(Sample,longitude,latitude,31:45)
envTab_native <- ps_melt %>% dplyr::filter(invasion == "Native") %>% dplyr::select(Sample,longitude,latitude,31:45)
envTab_transition <- ps_melt %>% dplyr::filter(invasion == "Transition") %>% dplyr::select(Sample,longitude,latitude,31:45)
envTab_invaded <- ps_melt %>% dplyr::filter(invasion == "Invaded") %>% dplyr::select(Sample,longitude,latitude,31:45)



# format for gdm
gdmTab <- formatsitepair(bioData=sppTab, 
                              bioFormat=2, #x-y spp list
                              XColumn="longitude", 
                              YColumn="latitude",
                              sppColumn="OTU", 
                              siteColumn="Sample", 
                              predData=envTab,
                              abundance = TRUE,
                              abundColumn = "Abundance")
gdmTab_native <- formatsitepair(bioData=sppTab_native, 
                         bioFormat=2, #x-y spp list
                         XColumn="longitude", 
                         YColumn="latitude",
                         sppColumn="OTU", 
                         siteColumn="Sample", 
                         predData=envTab_native,
                         abundance = TRUE,
                         abundColumn = "Abundance")
gdmTab_transition <- formatsitepair(bioData=sppTab_transition, 
                         bioFormat=2, #x-y spp list
                         XColumn="longitude", 
                         YColumn="latitude",
                         sppColumn="OTU", 
                         siteColumn="Sample", 
                         predData=envTab_transition,
                         abundance = TRUE,
                         abundColumn = "Abundance")
gdmTab_invaded <- formatsitepair(bioData=sppTab_invaded, 
                         bioFormat=2, #x-y spp list
                         XColumn="longitude", 
                         YColumn="latitude",
                         sppColumn="OTU", 
                         siteColumn="Sample", 
                         predData=envTab_invaded,
                         abundance = TRUE,
                         abundColumn = "Abundance")

# fit GDM models
gdm <- gdm(data = gdmTab,geo = TRUE)
gdm_native <- gdm(data = gdmTab_native,geo = TRUE)
gdm_transition <- gdm(data = gdmTab_transition,geo = TRUE)
gdm_invaded <- gdm(data = gdmTab_invaded,geo = TRUE)

# quick look at model fits
summary(gdm)

# predictions from model (using same distances)
gdm_pred <- predict(object=gdm, data=gdmTab)
gdm_pred_native <- predict(object=gdm_native, data=gdmTab_native)
gdm_pred_transition <- predict(object=gdm_transition, data=gdmTab_transition)
gdm_pred_invaded <- predict(object=gdm_invaded, data=gdmTab_invaded)


preds <- data.frame(observed = gdmTab$distance,
                         predicted = gdm_pred,
                         dist = gdm$ecological,gdmTab)
preds_native <- data.frame(observed = gdmTab_native$distance,
                    predicted = gdm_pred_native,
                    dist = gdm_native$ecological,gdmTab_native)
preds_transition <- data.frame(observed = gdmTab_transition$distance,
                    predicted = gdm_pred_transition,
                    dist = gdm_transition$ecological,gdmTab_transition)
preds_invaded <- data.frame(observed = gdmTab_invaded$distance,
                    predicted = gdm_pred_invaded,
                    dist = gdm_invaded$ecological,gdmTab_invaded)

p1 <- preds %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() + ggtitle("Overall")
p2 <- preds_native %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() + ggtitle("Native")
p3 <- preds_transition %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() + ggtitle("Transition")
p4 <- preds_invaded %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() + ggtitle("Invaded")
(p1 + p2) / (p3 + p4)
ggsave("./output/figs/18S_GDM_Prediction_Fit.png")


## Extract splines ####
isplines <- gdm::isplineExtract(gdm) %>% as.data.frame()
isplines_native <- gdm::isplineExtract(gdm_native) %>% as.data.frame()
isplines_transition <- gdm::isplineExtract(gdm_transition) %>% as.data.frame()
isplines_invaded <- gdm::isplineExtract(gdm_invaded) %>% as.data.frame()

names(isplines)
newnames <- c("geographic_actual","NO3_actual","NH4_actual","K_actual","P_actual","OM_actual",
              "pH_actual","Ca_actual","Mg_actual","Na_actual","SO4_actual","B_actual",
              "Zn_actual","Mn_actual","Cu_actual","Fe_actual",
              "geographic_partial","NO3_partial","NH4_partial","K_partial","P_partial","OM_partial",
              "pH_partial","Ca_partial","Mg_partial","Na_partial","SO4_partial","B_partial",
              "Zn_partial","Mn_partial","Cu_partial","Fe_partial")
names(isplines) <- newnames
names(isplines_native) <- newnames
names(isplines_transition) <- newnames
names(isplines_invaded) <- newnames

isplines <- 
isplines %>% 
  mutate(group="Overall") %>% 
  full_join(mutate(isplines_native,group="Native")) %>% 
  full_join(mutate(isplines_transition,group="Transition")) %>% 
  full_join(mutate(isplines_invaded,group="Invaded")) %>% 
  mutate(group=factor(group,levels=c('Overall','Native','Transition',"Invaded")))



## Plot splines ####

gdm %>% summary

top_predictors <- # contributing at least 5% explanatory power overall
  data.frame(
    coefs=gdm$coefficients,
    predictor=rep(gdm$predictors,each=3)
  ) %>% 
  group_by(predictor) %>% 
  summarize(sum_coef=sum(coefs)) %>% 
  arrange(desc(sum_coef)) %>% 
  dplyr::filter(sum_coef > 0.05 | predictor == "Geographic") # but also include geog dist
top_predictors$predictor

p1 <- isplines %>% 
  ggplot(aes(x=geographic_actual*10000,y=geographic_partial,color=group)) +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  labs(y="f(Geographic distance)",x="Geographic distance (m)",color="Taxa") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p2 <- isplines %>% 
  ggplot(aes(x=Ca_actual,y=Ca_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(Ca)",x="Ca",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p3 <- isplines %>% 
  ggplot(aes(x=NO3_actual,y=NO3_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(NO3)",x="NO3",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p4 <- isplines %>% 
  ggplot(aes(x=P_actual,y=P_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(P)",x="P",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p5 <- isplines %>% 
  ggplot(aes(x=pH_actual,y=pH_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(pH)",x="pH",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p6 <- isplines %>% 
  ggplot(aes(x=Mg_actual,y=Mg_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(Mg)",x="Mg",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")


p7 <- isplines %>% 
  ggplot(aes(x=Fe_actual,y=Fe_partial,color=group)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75,aes(linetype=group)) +
  scale_color_manual(values=c("gray20",pal$pal.earthtones)) +
  # facet_wrap(~group) +
  labs(y="f(Fe)",x="Fe",color="Treatment") +
  scale_linetype_manual(values = c(21,1,1,1)) + guides(linetype = "none")



(p2 | p3) / (p4 | p5) / (p6 | p7) + plot_layout(guides='collect')
ggsave("./output/figs/GDM_splines_plots_soilvars.png",dpi=400,height = 8,width = 12)



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
  scale_color_manual(values = c("black","gray","#F0027F")) +
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
mergevar <- paste(ps@sam_data$site,ps@sam_data$invasion,sep="_")
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
  facet_wrap(~invasion,scales = 'free_x') +
  scale_fill_viridis_d(option='turbo')
ggsave("./output/figs/significantly_different_taxa_relabund_barplot.png", height = 8, width = 10)

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

sig_melt
joined_mods <- 
  joined_mods %>% 
  mutate(term = term %>% str_remove_all("mu.") %>% 
           str_remove_all("invasion") %>% 
           str_remove_all("site")) %>% 
  mutate(across(where(is.numeric),function(x){round(x,4)})) 
joined_mods %>% 
  write_csv("./output/significant_taxa_stats_table.csv")
joined_mods %>% 
  mutate(taxon = taxon %>% 
           sub(".*_", "", .))
