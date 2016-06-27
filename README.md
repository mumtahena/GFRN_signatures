# GFRN_signatures

This repository presents all analyses done for the paper "Discrete growth and survival phenotypes link to apoptotic mechanism and drug response in breast cancer." `ASSIGN` folder describes pathway prediction generation in cell lines and in patient samples.`GFRN_characterization_in_breast_cancer_Final.Rmd`script has code for all the results described in the paper. 

# Datasets

You will need the following datasets for running the `GFRN_characterization_in_breast_cancer_Final.Rmd` script:

1. ICBP breast cancer cell line gene expression dataset: [icbp_Rsubread_tpmlog.txt](https://www.dropbox.com/sh/moyt4evz0cowbl1/AACZKh5uti9Lsc4xv2dYhRgNa?dl=0)
2. ICBP breast cancer cell line drug response dataset: [ICBP_drugs.txt](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0658-5/MediaObjects/13059_2015_658_MOESM2_ESM.xlsx)
3. TCGA cancer tumor gene expression dataset: [GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1536837&format=file&file=GSM1536837%5F06%5F01%5F15%5FTCGA%5F24%2Etumor%5FRsubread%5FTPM%2Etxt%2Egz)
4. TCGA PAM50 classification of tumor samples: [BRCA.547.PAM50.SigClust.Subtypes.txt](https://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt)
5. TCGA clinical dataset: [GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz](http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE62944&format=file&file=GSE62944%5F06%5F01%5F15%5FTCGA%5F24%5F548%5FClinical%5FVariables%5F9264%5FSamples%2Etxt%2Egz) 
6. ICBP cell line single pathway optimized predictionL:optimized_single_pathway_icbp.txt (in `Datasets` downloadable directory below)
7. TCGA cell line single pathway optimized prediction: optimized_single_pathway_tcga.txt (in `Datasets` downloadable directory below)
8. Drug response assay: Drug_response_assay.txt (in `Datasets` downloadable directory below) 
9. Gene Expression Data for GFRN Signatures: [GSE83083_GFP18_AKT_BAD_HER2_IGF1R_RAF_GFP30_KRAS-G12V_tpmlog.txt.gz] (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?token=avanmkmqrnmrdwn&acc=GSE83083)

To download all the required files click [Datasets](https://www.dropbox.com/sh/ltfubdiodti5yx0/AAAuVRh34mOOQYq7s7jF6IQJa?dl=0).

# Required R packages
1. gplots
2. RColorBrewer
3. data.table
4. mclust
5. ggplot2
6. gridExtra
7. ComplexHeatmap
8. Cluster

# Other Required Files:
`Key_ASSIGN_functions_balancedsig.R` - This file is available in this repository.
