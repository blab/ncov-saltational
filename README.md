# Saltational events in SARS-CoV-2

## Abstract

## Repository organization
This repository is organized around a Snakemake workflow that enables to process and analyze the viridian data, download SRA metadata, run the Bayesian model we developed, identify adaptive branches, visualize the adaptive branches in the context of the global SARS-CoV-2 phylogeny and investigate mutation patterns in adaptive branches. 

This repository is organized in sub-folders as follows:
- ```figures/``` contains the figures and the scripts to generate the figures (both from the main text and the supplementary information) associated with the manuscript. Further information is available on the folder-level [README]() file. 
- ```manuscript/``` contains the manuscript.
- ```scripts/``` contains the custom scripts used in the workflow. Further information is available on the folder-level [README]() file. 
- ```workflow/``` contains the Snakemake files and environment files used in the main Snakefile.
