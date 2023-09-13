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

* Mode 1:
```
code blocks for commands
```

* Mode 2:
```
code blocks for commands
```

## Authors

* Contributors names and contact info:

Yidi Qin \
yiq22@pitt.edu


Soyeon Kim \
soyeon.kim21@chp.edu

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments
