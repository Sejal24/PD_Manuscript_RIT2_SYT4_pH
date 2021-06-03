library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

##Section: Sex-specific signals

#Function for three way interaction between gene expression, genotype, and gender
interaction_SexGeneSNP <- function(gene) {
  lm_all <- lm(as.formula(paste("ph ~ rs12456492*sex*",  gene,paste("+ PMI + RIN + age +sex"))), filter(BE_allRegions,region=="substantia nigra" & !phOutlier))  
  sum_mod <- summary(lm_all)
  confid_inter <- confint(lm_all, level=0.95)
  print(gene )
  print(sum_mod)
  print(confid_inter)
}

#run interaction_SexGeneSNP function for each gene
interaction_SexGeneSNP(gene = "rit2") 
interaction_SexGeneSNP(gene = "syt4") 

#function for SNP gene expression interaction, stratify based on gender
gender_GeneSNP_interaction <- function(gene) {
  male <- lm(as.formula(paste("ph ~ rs12456492*",  gene,paste("+ PMI + RIN + age"))), filter(BE_allRegions,region=="substantia nigra", !phOutlier,sex==0))
  sum_mod_male <- summary(male)
  confid_inter_male <- confint(male, level=0.95)
  print(paste0("Male: ", gene ))
  print(sum_mod_male)
  print(confid_inter_male)
  female <- lm(as.formula(paste("ph ~ rs12456492*",  gene,paste("+ PMI + RIN + age"))), filter(BE_allRegions,region=="substantia nigra", !phOutlier,sex==1))
  sum_mod_female <- summary(female)
  confid_inter_female <- confint(female, level=0.95)
  print(paste0("Female: ", gene ))
  print(sum_mod_female)
  print(confid_inter_female)
}

#run gender_GeneSNP_interaction function for each gene
gender_GeneSNP_interaction(gene = "rit2") 
gender_GeneSNP_interaction(gene = "syt4") 


#Figure 4: Scatter plot of pH and gene expression grouped by risk SNP genotype and stratified sex for RIT2 and SYT4 

RIT2_MF<- ggplot(filter(BE_allRegions,region=="substantia nigra", !phOutlier), aes(x = rit2, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point() + 
  facet_wrap( ~ sex, labeller= labeller(sex=c("0" = "Male", "1" = "Female"))) +
  geom_smooth(method= "lm") + 
  theme_bw() + 
  xlab("RIT2 Expression") + ylab("pH") + scale_colour_discrete(name = "rs12456492") 

SYT4_MF <- ggplot(filter(BE_allRegions,region=="substantia nigra", !phOutlier), aes(x = syt4, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point() + 
  facet_wrap( ~ sex, labeller= labeller(sex=c("0" = "Male", "1" = "Female"))) +
  geom_smooth(method= "lm") + 
  theme_bw() + 
  xlab("SYT4 Expression") + ylab("pH") + scale_colour_discrete(name = "rs12456492") 

#add RIT2 and SYT4 plot side by side 
RIT2_SYT4_pH_Genderplot <- plot_grid(RIT2_MF + theme(legend.position = "none") , SYT4_MF +theme(legend.position = "none"), labels = c('A', 'B'), label_size = 12, nrow=2)

#add legend 
Figure4 <- plot_grid(RIT2_SYT4_pH_Genderplot, get_legend(SYT4_MF + theme(legend.position="bottom")), ncol=1, rel_heights = c(0.5, 0.05))

ggsave(Figure4, filename = here('Figures','Figure 6 Gender RIT2 and SYT4 gene expression vs pH in SN.pdf'), dpi=300)


##NOTE: Above analysis can also be run with GTEx data