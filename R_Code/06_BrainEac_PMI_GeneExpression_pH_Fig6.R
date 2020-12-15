library(here)
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

##Section: PMI influences pH-eQTL strength

#three-way interaction between SNP, gene expression and PMI in BrainEac data
interaction_PMIGeneSNP <- function(gene) {
  lm_all <- lm(as.formula(paste("ph ~ rs12456492*PMI*",  gene,paste("+ RIN + age +sex"))), filter(BE_allRegions,region=="substantia nigra" & !phOutlier))  
  sum_mod <- summary(lm_all)
  confid_inter <- confint(lm_all, level=0.95)
  print(gene)
  print(sum_mod)
  print(confid_inter)
}

#run interaction_PMIGeneSNP function for each gene
interaction_PMIGeneSNP(gene = "rit2") 
interaction_PMIGeneSNP(gene = "syt4") 


#function for SNP*gene interaction, stratify based on PMI group (group 1 PMI =< 49h, group 2 = PMI > 49h)
PMI_interactionGeneSNP <- function(gene) {
  short_pmi <- lm(as.formula(paste("ph ~ rs12456492*",  gene,paste("+ sex + RIN + age"))), filter(BE_allRegions,region=="substantia nigra", !phOutlier,PMI_Group==1))
  sum_mod_short_pmi <- summary(short_pmi)
  confid_inter_short_pmi <- confint(short_pmi, level=0.95)
  print(paste0("Short_PMI: ", gene ))
  print(sum_mod_short_pmi)
  print(confid_inter_short_pmi)
  long_pmi <- lm(as.formula(paste("ph ~ rs12456492*",  gene,paste("+ sex + RIN + age"))), filter(BE_allRegions,region=="substantia nigra", !phOutlier,PMI_Group==2))
  sum_mod_long_pmi <- summary(long_pmi)
  confid_inter_long_pmi <- confint(long_pmi, level=0.95)
  print(paste0("Long_PMI: ", gene ))
  print(sum_mod_long_pmi)
  print(confid_inter_long_pmi)
}

#run interaction_PMIGeneSNP function for each gene
PMI_interactionGeneSNP(gene = "rit2") 
PMI_interactionGeneSNP(gene = "syt4") 


#Figure 6 Scatter plot of pH and gene expression grouped by risk SNP genotype and stratified PMI within BrainEac sample RIT2 and SYT4

RIT2_PMI <- ggplot(filter(BE_allRegions,region=="substantia nigra",!phOutlier), aes(x = rit2, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point(size=0.8) + 
  facet_wrap( ~ PMI_Group, labeller= labeller(PMI_Group=c("1" = "PMI =< 49h", "2" = "PMI > 49h"))) +
  geom_smooth(method= "lm", size = 0.6) + 
  theme_bw() + 
  xlab("RIT2 Expression") + ylab("pH") + scale_colour_discrete(name = "rs12456492") 

SYT4_PMI <- ggplot(filter(BE_allRegions,region=="substantia nigra",!phOutlier), aes(x = syt4, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point(size=0.8) + 
  facet_wrap( ~ PMI_Group, labeller= labeller(PMI_Group=c("1" = "PMI =< 49h", "2" = "PMI > 49h"))) +
  geom_smooth(method= "lm", size = 0.6) + 
  theme_bw() + 
  xlab("SYT4 Expression") + ylab("pH") + scale_colour_discrete(name = "rs12456492") 

#add RIT2 and SYT4 plot side by side 
RIT2_SYT4_pH_PMIplot<- plot_grid(RIT2_PMI + theme(legend.position = "none") , SYT4_PMI +theme(legend.position = "none"), labels = c('A', 'B'), label_size = 12, nrow=2)

#add legend 
Figure6 <- plot_grid(RIT2_SYT4_pH_PMIplot, get_legend(SYT4_PMI + theme(legend.position="bottom")), ncol=1, rel_heights = c(1, 0.1))

ggsave(Figure6, filename = here('Figures', 'Figure 6 PMI RIT2 and SYT4 gene expression vs pH in SN.pdf'), dpi=300)
