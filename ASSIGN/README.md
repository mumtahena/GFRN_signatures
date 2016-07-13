ASSIGN Analysis Scripts
=======================

To recreate the ASSIGN analysis performed on the ICBP and TCGA test data, you
can use the scripts provided here. Please contact [David Jenkins](mailto:dfj@bu.edu)
and [Mumtahena Rahman](mailto:mumtahena@gmail.com) with questions or issues.

> Before running the scripts, you must modify the file input and output locations
> in all scripts

## Input files

The same ASSIGN scripts can be used to analyze the 55 sample ICBP dataset and
the 1119 sample TCGA BRCA dataset used in the paper. Modify the path for the
```testFile``` parameter to the location of the file on your system.

## R Packages

Before running the scripts, make sure the following R packages are installed:

__CRAN Packages__:

1. [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
2. [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

__Development Packages__:

Development versions of _SVA_ and _ASSIGN_ were used for the analysis in this
paper. Please install the development versions of these packages, available
on GitHub:

```
library(devtools)
install_github("wevanjohnson/sva-devel", ref="2e87233c3797c3a630d3cf662e2ee08dfe42c7b5")
install_github("wevanjohnson/ASSIGN", ref="fa6621a31629d9a122b9a2bcaa724f7d8f2c02ce")
```

## Running the ASSIGN Scripts

The script ```ASSIGN_merge_and_combat.R``` will load the data, run batch
correction, and save an rsession that can be loaded by
```ASSIGN_run_predictions_single.R```. Before running either script, modify the file
paths in the _Input Files_ and _Output Files_ section of the scripts to
locations on your system. Once the file paths have been set, you can run the
analysis with the following commands:

> Run the optimized pathway lengths with the commands below. For each pathway, modify
> the second parameter of ASSIGN_run_predictions_single.R  to change the number 
> of genes in the signature:

```
# Load, merge, and ComBat data:
Rscript ASSIGN_merge_and_combat.R

# Run predictions for AKT, 20 genes
Rscript ASSIGN_run_predictions_single.R 1 20

# Run predictions for BAD, 250 genes
Rscript ASSIGN_run_predictions_single.R 2 250

# Run predictions for EGFR, 50 genes
Rscript ASSIGN_run_predictions_single.R 3 50

# Run predictions for HER2, 10 genes
Rscript ASSIGN_run_predictions_single.R 4 10

# Run predictions for IGF1R, 100 genes
Rscript ASSIGN_run_predictions_single.R 5 100

# Run predictions for KRASGV, 175 genes
Rscript ASSIGN_run_predictions_single.R 6 175

# Run predictions for RAF, 350 genes
Rscript ASSIGN_run_predictions_single.R 9 350
```

The pathway prediction values for each ASSIGN run can be found in the
```pathway_activity_testset.csv``` file in each pathway's subdirectory.

## Protein and Drug Correlation Scripts

After running ASSIGN, we correlate the pathway activity predictions to drug response and RPPA data

> Run the following scripts in your ASSIGN working directory:

```
Rscript icbp_cor.R
Rscript icbp_prot_analysis.R
Rscript tcga_prot_analysis.R
```

## Other Required Files:

* ```Key_ASSIGN_functions_balancedsig.R``` - Internal R functions used by
```ASSIGN_merge_and_combat.R```,```ASSIGN_run_predictions_single.R``` and the correlation scripts. This file is available in this repository
* ICBP breast cancer cell line drug response dataset: ```ICBP_drugs.txt``` for ```icbp_cor.R```
* ICBP breast cancer cell line log2(TPM) expression data: ```icbp_Rsubread_tpmlog.txt```
* TCGA BRCA log2(TPM) expression data: ```PANCAN24_BRCA_1119_TPMlog2.txt```
* ICBP breast cancer cell line RPPA data: ```proteomics.txt``` for ```icbp_prot_analysis.R```
* TCGA breast cancer RPPA data: ```TCGA-BRCA-RBN.csv``` for ```tcga_prot_analysis.R```

These files are available [here](https://www.dropbox.com/sh/ltfubdiodti5yx0/AAAuVRh34mOOQYq7s7jF6IQJa?dl=0)
