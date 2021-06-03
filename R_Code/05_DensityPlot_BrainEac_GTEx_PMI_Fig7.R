library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

##Section: PMI influences pH-eQTL strength

##NOTE: need to obtain GTEx data for this script to run

#Figure 5  Density plot of PMI values in hours of GTEx and BrainEac samples

##code includes GTEx data that is not publicly available, when data is obtained 
##and structured the same way as the BrainEac data table, run code below for density PMI plot   

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)
#Load GTEx data
#GT <- read.csv(here("Data","GTEx_SN_data.csv"), header = T)


pmiBE <- BE_allRegions %>% filter(region=="substantia nigra",!phOutlier) %>% select(PMI) %>% mutate( type = "BE")
pmiGT <- GT %>% filter(!phOutlier) %>% select(PMI) %>% mutate( type = "GT")
plot_pmi_df <- rbind(pmiBE, pmiGT)

Figure7 <- ggplot(plot_pmi_df,aes(x=PMI, fill=type)) + 
  geom_histogram(aes(y=(..density..),  binwidth = 0.05, position = "dodge", kernel="r" )) +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size = 14)) +
  labs(y="Density")+
  scale_fill_discrete(labels = c("BrainEAC", "GTEx"))

ggsave(Figure7, filename = here('Figures','Figure 7 PMI distribution.pdf'), dpi=300)

