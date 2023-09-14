# Methyl-TWAS

Methyl-TWAS is a powerful method for precise in silico transcriptome-wide association studies (TWAS) using long-range DNA methylation.

## Description

We propose Methyl-TWAS, the first in silico TWAS method utilizing long-range DNA methylation instead of genotype data by three steps. In step 1, Methyl-TWAS estimates weights of the long-range DNA methylation per gene g using the elastic-net penalized regression models using the training data set. Methyl-TWAS can estimate epigenetically regulated/associated expression (eGReX), which includes a part of the genetically regulated expression (GReX) and environmentally regulated expression,  trait-altered expression, tissue-specific expression since DNA methylation is affected by all of these factors. Note that the traditional TWAS methods only can estimate genetically regulated expression (GReX). In step 2, once we drive cis- and trans- (<10MB) weights of DNA methylation levels at CpG sites for each gene, Methyl-TWAS predicts eGReX using the estimated coefficients in the first step and DNA methylation data in a test data set. In step 3, using the imputed gene expression from step 2, Methyl-TWAS conducts TWAS to identify genes that are differentially expressed in a disease or a phenotype using the accurately predicted gene expressions.  

## Getting Started

### Dependencies

* Libraries used in the package:
```
library(utils)
library(parallel)
library(glmnet)
library(limma)
```

### Installing

* How/where to download the package:
```
install.packages("devtools")
library(devtools)
install_github("YidiQin/MethylTWAS", force = TRUE)
library(MethylTWAS)
```
* Any modifications needed to be made to files/folders

### Executing program

* How to run Methyl-TWAS: \
The public available Methyl-TWAS R-package offers 2 modes. In the first mode, users can train the model and conduct TWAS using any type of tissue. Users input the training (reference) data consisting of DNA methylation and gene expression data for the same subjects to estimate tissue-specific effects of DNA methylation on expression levels of each gene and then input testing data of DNA methylation and a phenotype to impute gene expression and then conduct in silico TWAS. In the second mode, users can skip step 1 and just conduct TWAS using a test data set (DNA methylation) of nasal epithelium (upper airway) tissue. Users input DNA methylation and a phenotype and then, Methyl-TWAS uses the estimated coefficients from our reference data, EVA-PR, to impute gene expression and conduct in silico TWAS without requiring further reference data.

* Example data: \
To help users gain a better understanding of input and output files, we provide example datasets generated from Yang et al. (GSE65205) which can be directly used after package installation. Example datasets include two subsetted methylation matrix (train.meth.1.rda, train.meth.2.rda), a subsetted gene expression matrix (train.exp.rda) as trainning data, a subsetted methylation matrix (test.meth.rda) as test data, and a phenotype data file (pheno.rda) to be used in TWAS analysis.

* Mode 1: \
By running function "MethylTWAS", users are able to input customized training data, test data, and phenotype data. To use example data in this mode, users can set "example = T" when running the function. An demo of running the function is showed below:
```
MethylTWAS(example = F,
           train.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.meth.rda",
           test.meth.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/test.meth.rda",
           train.exp.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/train.exp.rda",
           pheno.file = "/ix/ksoyeon/YQ/code/MethylTWAS/data/pheno.rda",
           output.file.path = "/ix/ksoyeon/YQ/results/test/",
           TWAS = T,
           phenotype = "cc_new",
           confounder = "gender, age")
```
There are two output files in the folder specified by output.file.path:
```
ls /ix/ksoyeon/YQ/results/test/
prediction.Rdata  TWAS.result.txt
```
prediction.Rdata is a matrix containing gene expression value predicted using train.meth.rda and train.exp.rda. TWAS.result.txt is a matrix containing TWAS results.

* Mode 2: \
By running function "MethylTWAS_BySummary", users are able to directly use pre-generated coefficients for prediction. To use example data in this mode, users can set "example = T" when running the function. An demo of running the function is showed below:
```
MethylTWAS_BySummary(example = F,
           test.meth.file = "/ix/ksoyeon/YQ/data/Yang.meth.rda",
           TWAS = T,
           pheno.file = "/ix/ksoyeon/YQ/data/Yang.pheno.rda",
           output.file.path = "/ix/ksoyeon/YQ/results/test_EVA_PR_on_Yang/",
           phenotype = "cc_new",
           confounder = "gender, age",
           core.num = 1)
```
There are two output files in the folder specified by output.file.path:
```
ls /ix/ksoyeon/YQ/results/test_EVA_PR_on_Yang/
prediction.Rdata  TWAS.result.txt
```
prediction.Rdata is a matrix containing gene expression value predicted using train.meth.rda and train.exp.rda. TWAS.result.txt is a matrix containing TWAS results.

## Authors

### Contributors names and contact info:

* Yidi Qin \
yiq22@pitt.edu


* Soyeon Kim \
soyeon.kim21@chp.edu

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT.A short and simple permissive license with conditions only requiring preservation of copyright and license notices. Licensed works, modifications, and larger works may be distributed under different terms and without source code.
