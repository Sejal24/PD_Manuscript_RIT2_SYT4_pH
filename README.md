The following R scripts were used for the analysis performed in the manuscript (https://www.biorxiv.org/content/10.1101/2020.12.16.423140v1). The analysis consisted of using the BrainEac and GTEx datasets. 

Data obtained from BrainEac database is publicly available (http://www.braineac.org/). To run the analysis, BrianEac data is already provided in the data folder. Raw data can be obtained at (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46706). Pre-processed steps were taken to create the data file. The data file in the Data folder is pre-formatted with demographics, gene expression (RIT2, SYT4 and CA10) and genotype data in one data frame (BrainEac_allRegions_data.csv). 

Demographic data from the GTEx database is restricted and request of data is needed. Therefore the R scripts are written based on BrainEac data analysis. Exact same analysis can be applied to GTEx data. Based on the restrictions obtaining GTEx data, we are not able to post pre-formatted data table in the Data folder. Information from the data files were obtained and pre-formatted to obtain similar data structure as seen in the BrainEac data table from the Data folder. For more information on the steps taken to obtain GTEx data and format the raw data please contact us. 

GTEx data can be obtained at the GTEx Portal (https://gtexportal.org/home/):
- Demographics data: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
- Gene Expression data: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct
- Genotype data for PD risk SNP (rs12456492): GTEx_Analysis_2017-06-05_v8_WholeExomeSeq_979Indiv_VEP_annot.vcf

Script can be run in the following order :
- 01_BrainEac_CorrelationTesting_Fig1.R
- 02_BrainEac_SNPGeneExpression_pH_AllBrainRegions_Fig2and3.R
- 03_BrainEac_SNPGeneExpression_pH_SubstantiaNigra.R
- 04_BrainEac_SexSpecfic_pH_Fig4.R
- 05_DensityPlot_BrainEac_GTEx_PMI_Fig5.R
- 06_BrainEac_PMI_GeneExpression_pH_Fig6.R
- 07_CA10_CoExpression_Analysis.R

Scripts 03 to 07 can also be used on GTEx data set

In script "07_CA10_CoExpression_Analysis", data for the top co-expression data set for RIT2 was obtained from Gemma with ID GSE20146 (https://gemma.msl.ubc.ca/home.html). 

