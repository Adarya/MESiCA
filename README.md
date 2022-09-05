# MESiCA

This repository contains code and custom scripts for the manuscript titled "Mutational Signatures identification in clinical-setting assays using neural embedding-based representations". 

MESiCA is a tool for identifying the dominant mutational signature in clinical-setting gene panels, even with only a few mutations available. 
This repository contains [MESiCA development](#Model) and [Validation schemes and downstream analyses](#Downstream). 

<a name="Model"/>

# Model training
In order to create representations, the main function to run is train.py. This script takes as input mutational classification by penta-nucleotide (1,536 possibilities), as well as sample tag, cancer type, and dominant mutational signature. 
The model creates representations for each entity, in such manner that related entities get closer representations (mathematically) than non-related entities. Please refer to the manuscript for a detailed description. 
The output of train.py is representations for each mutational signature, mutation class, cancer type and sample tag.

The overall scheme and code for creating representations was adopted and modified from [MutSpace](https://github.com/ma-compbio/MutSpace).

<a name="Downstream"/>

# Downstream
The Downstream directory contains codes and notebooks describing the workflow of:
  1. Independent validation of MESiCA's prediction capabilities in various Whole Exome and Whole Genome samples and settings.
  2. Predictions of mutational signatures in cancer gene panels cohorts.
  3. Survival analysis using our predictions in various indications. 

# Notes
- This repository will become publicly available once the manuscript is published
- Other supporting data are availabe at [doi.org/10.5281/zenodo.7047369](https://doi.org/10.5281/zenodo.7047369) with restricted access until publication.
