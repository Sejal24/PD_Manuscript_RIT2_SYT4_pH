library(here)
library(tidyr)
library(dplyr)

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

##Section: Rs12456492 influences the association between RIT2 and SYT4 expression and pH in the substantia nigra

#function for linear model between pH and averaged expression of RIT2, and STY4 with risk SNP for each brain region excluding outliers 
# store data in a table format: region tested, Beta value, P value, adjusted p value (fdr)
gene_SNP_region <- function(gene) {
  GeneSNP_LinearModel_table = data.frame()
  for(regionTested in unique(BE_allRegions$region)) {
    lm_Sum <- summary(lm(as.formula(paste("ph ~ rs12456492*",  gene)), filter(BE_allRegions,region==regionTested & !phOutlier)))  
    B_value<- lm_Sum$coefficients[(paste0("rs12456492:",gene)), "Estimate"]
    p_value <-lm_Sum$coefficients[(paste0("rs12456492:",gene)), "Pr(>|t|)"]
    df_row <- data.frame(regionTested, B_value, p_value)
    GeneSNP_LinearModel_table <- rbind(GeneSNP_LinearModel_table,df_row)
  }
  GeneSNP_LinearModel_table$adjusted_p <- p.adjust(GeneSNP_LinearModel_table$p_value, method = "fdr")
  return(GeneSNP_LinearModel_table)
}

#run gene_SNP_region function for each gene
gene_SNP_region(gene = "rit2") 
gene_SNP_region(gene = "syt4") 

#Figure 2 and 3: Scatter plot of pH and gene expression based on rs12456492 genotype for 10 brain regions

RIT2_region <- ggplot(filter(BE_allRegions,!phOutlier), aes(x = rit2, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point(size=0.3) + facet_wrap( ~ region) + 
  geom_smooth(method= "lm", size= 0.5) + theme_bw() + xlab("RIT2 Expression") + ylab("pH") + 
  scale_colour_discrete(name = "rs12456492") + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))

ggsave(RIT2_region, filename = here('Figures','Figure 4 pH with RIT2 gene expression across all brain regions.pdf'), height = 6, width = 8.25, dpi=300)

SYT4_region <-ggplot(filter(BE_allRegions,!phOutlier), aes(x = syt4, y = ph, group= rs12456492Named,color=rs12456492Named)) +
  geom_point(size=0.3) + facet_wrap( ~ region) + 
  geom_smooth(method= "lm", size= 0.5) + theme_bw() + xlab("SYT4 Expression") + ylab("pH") + 
  scale_colour_discrete(name = "rs12456492") + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))

ggsave(SYT4_region, filename = here('Figures', 'Figure 5 pH with SYT4 gene expression across all brain regions.pdf'),  height = 6, width = 8.25, dpi=300)
