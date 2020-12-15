library(tidyr)
library(here)
library(ggplot2)
library(dplyr)

#load data, BrainEac csv data file "BrainEac_allRegions_data.csv" from Data folder:
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

##Section: Rs12456492 influences the association between RIT2 and SYT4 expression and pH in the substantia nigra

#correlation between pH and PD risk SNP in substantia nigra without covariates (age, sex, RIN, PMI) 
summary(lm(ph ~ rs12456492, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))

#correlation between pH and gene expression (RIT2 or SYT4) with PD risk SNP in substantia nigra without covariates (age, sex, RIN, PMI) 
summary(lm(ph ~ rs12456492*rit2, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))
summary(lm(ph ~ rs12456492*syt4, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))

#Table 1: linear models 1 to 6 accounitng for covariates of sex, age, PMI and RIN for both gene expression in substantia nigra 
#Each model gives summary statistics (summary()) and confidence intervals (confint())

Model_1 <- lm(ph ~ rs12456492 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_1)
confint(Model_1, level=0.95)

Model_2 <- lm(ph ~ rs12456492*rit2 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_2)
confint(Model_2, level=0.95)

Model_3 <- lm(ph ~ rs12456492*syt4 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_3)
confint(Model_3, level=0.95)

Model_4 <- lm(ph ~ rs12456492*rit2 + syt4 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_4)
confint(Model_4, level=0.95)

Model_5 <- lm(ph ~ rs12456492*syt4 + rit2 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_5)
confint(Model_5, level=0.95)

Model_6 <- lm(ph ~ rs12456492*rit2 + rs12456492*syt4 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier))
summary (Model_6)
confint(Model_6, level=0.95)


##NOTE: Above analysis can also be run with GTEx data