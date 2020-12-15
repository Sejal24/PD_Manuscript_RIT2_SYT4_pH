library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
#includes demographics, gene expression, and genotype data    
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

## Section: Brainwide RIT2 and SYT4 gene expression is correlated with pH 

#we need to reformat the data to obtain mean averages of each gene expression across all 10 brain regions for each donor
#store the data into a new dataframe "BE_donorStats"
BE_donorStats <- group_by(BE_allRegions, donor) %>% dplyr::summarise(
  count = n(),
  ph = mean(ph, na.rm = TRUE),
  rit2 = mean(rit2, na.rm = TRUE),
  syt4 = mean(syt4, na.rm = TRUE),
  ca10 = mean(ca10, na.rm = TRUE),
  rs12456492 = first(rs12456492),
  rs12456492Named = first(rs12456492Named),
  sex = first(sex),
  PMI = first(PMI),
  RIN = mean(RIN),
  age = first(age),
  phOutlier = first(phOutlier),
  brainBank = first(brainBank),
  hasSN = "substantia nigra" %in% region,
  G_allele = first(G_allele), 
  Allele_Named = first(Allele_Named)
)

#Can write data into a csv file, for future use:
write.csv(BE_donorStats, file = here("Data","BrainEac_DonorStats_MeanExp.csv"), row.names = F)

#Pearson Correlation between pH and averaged expression of RIT2, and STY4 across all 10 brain regions including outliers 
cor.test(BE_donorStats$rit2, BE_donorStats$ph,m='p')
#cor = 0.5851373, p-value < 0.0001
cor.test(BE_donorStats$syt4, BE_donorStats$ph,m='p')
#cor = 0.5819768, p-value < 0.0001

#Pearson Correlation between pH and averaged expression of RIT2, and STY4 across all 10 brain regions excluding outliers 
cor.test(filter(BE_donorStats, !phOutlier)$rit2, filter(BE_donorStats, !phOutlier)$ph,m='p')
#cor = 0.220528, p-value = 0.03175
cor.test(filter(BE_donorStats, !phOutlier)$syt4, filter(BE_donorStats, !phOutlier)$ph,m='p')
#cor = 0.1431437, p-value = 0.1664

#Pearson Correlation between pH and averaged expression of RIT2, and STY4 for each brain region excluding outliers 
# store data in a table format: region tested, correlation value, P value, adjusted p value (fdr) 
correlation_function_region <- function(gene) {
  correlation_table = data.frame()
  for(regionTested in unique(BE_allRegions$region)) {
    sumCor <- cor.test(as.formula(paste("~", gene, paste("+ ph"))), filter(BE_allRegions,region==regionTested & !phOutlier), m="p")  
    r_value <-sumCor$estimate
    p_value<- sumCor$p.value
    df_row <- data.frame(regionTested, r_value, p_value)
    correlation_table <- rbind(correlation_table,df_row)
  }
  correlation_table$adjusted_p <- p.adjust(correlation_table$p_value, method = "fdr")
  return(correlation_table)
}

#run function for each gene
correlation_function_region(gene = "rit2") 
correlation_function_region(gene = "syt4") 


#Figure 1 Scatter plots of the relationship between pH with RIT2 and SYT4 gene expression

RIT2_avg <- ggplot(filter(BE_donorStats), aes(x = rit2, y = ph, colour=phOutlier)) + geom_point() + xlab("Mean RIT2 expression across regions") + ylab("Brain pH") +  
  scale_color_manual(values=c("blue", "red"), name= NULL, breaks=c(TRUE), labels=c("pH outliers"))+
  geom_smooth(method='lm', data =filter(BE_donorStats), colour="red")+
  geom_smooth(method='lm', data =filter(BE_donorStats,!phOutlier), colour="blue")+
  theme_bw()

SYT4_avg <- ggplot(filter(BE_donorStats), aes(x = syt4, y = ph, colour=phOutlier)) + geom_point() + xlab("Mean SYT4 expression across regions") + ylab("Brain pH") +  
  scale_color_manual(values=c("blue", "red"))+
  geom_smooth(method='lm', data =filter(BE_donorStats), colour="red")+
  geom_smooth(method='lm', data =filter(BE_donorStats,!phOutlier), colour="blue")+
  theme_bw()

#add RIT2 and SYT4 plot side by side 
RIT2_SYT4_pH_plot<- plot_grid(RIT2_avg + theme(legend.position = "none") , SYT4_avg +theme(legend.position = "none"), labels = c('A', 'B'), label_size = 12)

#add legend 
Figure1 <- plot_grid(RIT2_SYT4_pH_plot, get_legend(RIT2_avg), 
                     nrow =1, rel_widths  = c(1, 0.15))

ggsave(Figure1, filename = here('Figures','Figure 1 pH with average RIT2 and SYT4 gene expression.pdf'), height = 4, width = 9.5, dpi=300)

