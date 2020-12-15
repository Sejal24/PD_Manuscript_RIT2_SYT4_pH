library(here)
library(tidyr)
library(dplyr)
library(data.table)

##Section: CA10 is co-expressed with SYT4 and RIT2
## Top ranked co-expression gene is CA10 and the top ranked database is GSE20146 study from the SEEK webtool
##GSE20146 study raw data was obtained form Gemma web tool

#read in the GSE20146 data from Data folder 
PD_data_exp <- read.csv(here('Data','1881_GSE20146_expmat.data.txt'), header = T, sep = "\t")

#extract probes for targeted genes 
gene_exp <- PD_data_exp %>% filter(GeneSymbol %in% c("RIT2", "SYT4", "CA10"))

# transpose the data table
gene_exp_T <- transpose(gene_exp)
colnames(gene_exp_T) <- rownames(gene_exp)
rownames(gene_exp_T) <- colnames(gene_exp)

#set index number
gene_exp_Tindex <- setDT(gene_exp_T, keep.rownames = TRUE)
gene_exp_Tindex_use <- gene_exp_Tindex %>% slice(-c(2,4:6,15))

#set column names
colnames(gene_exp_Tindex_use ) <- gene_exp_Tindex_use[1, ]
gene_exp_Tindex_use <- gene_exp_Tindex_use[-1, ] 
donor_pheno <- gene_exp_Tindex_use %>% rename( Donor = "Probe", RIT2 = "206984_s_at", SYT4 = "223529_at", CA10a = "220889_s_at", CA10b = "223550_s_at")
donor_pheno <- donor_pheno[-1, ] 

#set expression values as numeric
donor_pheno$CA10a <- as.numeric(donor_pheno$CA10a)
donor_pheno$CA10b <- as.numeric(donor_pheno$CA10b)
donor_pheno$RIT2 <- as.numeric(donor_pheno$RIT2)
donor_pheno$SYT4 <- as.numeric(donor_pheno$SYT4)

#average the two CA10 expression values into one value
donor_pheno$CA10=rowMeans(donor_pheno[,c("CA10a","CA10b")], na.rm=TRUE)

#extract data from text from Donor column
#extract ID, and Disease state 

reformat_donor_pheno <- donor_pheno %>% 
  mutate(ID = as.character(gsub("GSE20146_Biomat_[0-9]|[0-9]__|_BioAssayId.|__|[0-9]*Name.", "", Donor))) %>% 
  mutate(ID = as.character(gsub("GPi", " GPi ", ID))) %>% 
  mutate(ID = as.character(gsub("Cm|Cf|Pm|Pf", "", ID))) %>% 
  mutate(disease = as.character(gsub("GSE20146_Biomat_[0-9]|[0-9]__|_BioAssayId.|__|[0-9]*Name.[0-9]*GPi", "", Donor)))

#code the disease as 0 (C) or 1 (P) in separate column called "disease_code"
reformat_donor_pheno[reformat_donor_pheno$disease == "Cm","disease_code"] <- 0
reformat_donor_pheno[reformat_donor_pheno$disease == "Pm","disease_code"] <- 1
reformat_donor_pheno[reformat_donor_pheno$disease == "Cf","disease_code"] <- 0
reformat_donor_pheno[reformat_donor_pheno$disease == "Pf","disease_code"] <- 1

#code the disease as  control or Parkinson's Disease in separate column called "disease_state"
reformat_donor_pheno[reformat_donor_pheno$disease == "Cm","disease_state"] <- "control"
reformat_donor_pheno[reformat_donor_pheno$disease == "Pm","disease_state"] <- "Parkinson's Disease"
reformat_donor_pheno[reformat_donor_pheno$disease == "Cf","disease_state"] <- "control"
reformat_donor_pheno[reformat_donor_pheno$disease == "Pf","disease_state"] <- "Parkinson's Disease"

#code sex as  male (0) female(1) in separate column
reformat_donor_pheno[reformat_donor_pheno$disease == "Cm","sex"] <- 0
reformat_donor_pheno[reformat_donor_pheno$disease == "Pm","sex"] <- 0
reformat_donor_pheno[reformat_donor_pheno$disease == "Cf","sex"] <- 1
reformat_donor_pheno[reformat_donor_pheno$disease == "Pf","sex"] <- 1

#create data frame and select variables needed for correlation analysis 
PD_gene_exp <- reformat_donor_pheno %>% select(ID, sex, disease_state, disease_state, disease_code, RIT2, SYT4, CA10) 

#Pearson correlation with all the samples 
cor.test(PD_gene_exp$RIT2, PD_gene_exp$CA10,m='p')
cor.test(PD_gene_exp$RIT2, PD_gene_exp$SYT4,m='p')
cor.test(PD_gene_exp$SYT4, PD_gene_exp$CA10,m='p')

#Pearson correlation of RIT2, SYT4 and CA10 in controls
cor.test(filter(PD_gene_exp, disease_code == 0)$RIT2, filter(PD_gene_exp, disease_code == 0)$CA10,m='p')
cor.test(filter(PD_gene_exp, disease_code == 0)$RIT2, filter(PD_gene_exp, disease_code == 0)$SYT4,m='p')
cor.test(filter(PD_gene_exp, disease_code == 0)$SYT4, filter(PD_gene_exp, disease_code == 0)$CA10,m='p')

#Pearson correlation of RIT2, SYT4 and CA10 in PD
cor.test(filter(PD_gene_exp, disease_code == 1)$RIT2, filter(PD_gene_exp, disease_code == 1)$CA10,m='p')
cor.test(filter(PD_gene_exp, disease_code == 1)$RIT2, filter(PD_gene_exp, disease_code == 1)$SYT4,m='p')
cor.test(filter(PD_gene_exp, disease_code == 1)$SYT4, filter(PD_gene_exp, disease_code == 1)$CA10,m='p')

#test relatioship between coexpressed and interaction with disease state
test_ca10 <- function(gene1, gene2) {
    lm_all <- lm(as.formula(paste(gene1, "~ disease_code*", gene2)), filter(PD_gene_exp)) 
    print(summary(lm_all))
}

#use function test_ca10 with different gene combinations
test_ca10("RIT2", "SYT4")
test_ca10("RIT2", "CA10")
test_ca10("SYT4", "CA10")
test_ca10("SYT4", "RIT2")
test_ca10("CA10", "RIT2")
test_ca10("CA10", "SYT4")


#Analysis between CA10 and pH relationship in BrainEac data
#Note: Below analysis can be also applied to GTEx data
BE_allRegions <- read.csv(here("Data","BrainEac_allRegions_data.csv"), header = T)

#Pearson Correlation between pH and  CA10 gene expression 
cor.test(filter(BE_allRegions, region=="substantia nigra", !phOutlier)$ca10, filter(BE_allRegions,region=="substantia nigra", !phOutlier)$ph,m='p')
#r = 0.1545715, p-value = 0.2014

#Linear models: between pH and  CA10 gene expression when added to the models in the substantia nigra 
summary(lm(ph ~ ca10 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))
summary(lm(ph ~ rs12456492 * rit2 + ca10 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))
summary(lm(ph ~ rs12456492 * syt4 + ca10 + PMI + RIN + age + sex, filter(BE_allRegions,region=="substantia nigra" & !phOutlier)))

