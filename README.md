# MESiCA

This repository contains code and custom scripts for the manuscript titled "Mutational Signatures identification in clinical-setting assays using neural embedding-based representations". 

MESiCA is a tool for identifying the dominant mutational signature in clinical-setting gene panels, even with only a few mutations available. 
This repository contains [MESiCA development](#Model) and [Validation schemes and downstream analyses](#Downstream). 

<a name="Model"/>

# System requirements 

### Software requirements
1. **Python 3.7.3.** 
    + **Packages:** pytorch 1.3.1, argparse 1.4.0, json 0.1.1, numpy 1.16.4, tqdm 4.32.1, pandas 1.2.4.
2. **R 4.1.0.** 
    + **Packages:** BiocManager 1.30.16, ggplot2 3.3.5, caret 6.0.90, dplyr 1.0.8, tibble 3.1.6, rtracklayer 1.52.1, reshape2 1.4.4, Biostrings 2.60.1, plyr 1.8.6, Rsamtools 2.8.0, VariantAnnotation 1.38.0, SomaticSignatures 2.28.0, MutationalPatterns 3.2.0, survival 3.2.11, survminer 0.4.9, BSgenome 1.60.0, gridExtra 2.3, precrec 0.12.9, pheatmap 1.0.12. 
3. **Jupyter Notebook**

### OS requirements
+ The embedding model has been tested on Linux: Ubuntu 20.04
+ The downstream analyses have been tested on Windows 10 and 11

### Hardware requirements
A standard computer with enough RAM to support the in-memory operations.


# Model training
### Installation
Once the above Python and associated packages are installed, the code can be used
```
git clone https://github.com/Adarya/MESiCA
```

### Running
+ In order to create representations, the main function to run is **train.py**. 
+ This script takes as input mutational classification by penta-nucleotide (1,536 possibilities), as well as sample tag, cancer type, and dominant mutational signature. A processed demo data is available in the repository.  
+ The model creates representations for each entity, in such manner that related entities get closer representations (mathematically) than non-related entities. Please refer to the manuscript for a detailed description. 
+ The output of train.py is representations for each mutational signature, mutation class, cancer type and sample tag.

The overall scheme and code for creating representations was adopted and modified from [MutSpace](https://github.com/ma-compbio/MutSpace).

<a name="Downstream"/>

# Downstream
The Downstream directory contains codes and notebooks describing the workflow of:
  1. Independent validation of MESiCA's prediction capabilities in various Whole Exome and Whole Genome samples and settings.
  2. Predictions of mutational signatures in cancer gene panels cohorts.
  3. Survival analysis using our predictions in various indications. 

### Installation
Once the above R environment and associated packages are installed, and the repository has been cloned, no further installation is needed.

### Running
+ The outputs of the training are the inputs for the **model_output_to_csv.ipynb**, which produce .csv files
+ These files are then the inputs for predicting the dominant signature in new samples, as well as for all downstream tasks.

# Notes
- This repository will become publicly available once the manuscript is published
- Other supporting data are availabe at [doi.org/10.5281/zenodo.7047369](https://doi.org/10.5281/zenodo.7047369) with restricted access until publication.
