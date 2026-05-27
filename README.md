# SARS-CoV-2 saltational events are recurrent and consistent with evolution during prolonged human infections

## Abstract

SARS-CoV-2 evolution is characterized by the gradual accumulation of mutations but has also been punctuated by the emergence of highly mutated variants. Whether such saltational jumps are a broad feature of SARS-CoV-2 evolution or rare anomalies remains unclear. Here, we perform a systematic analysis of SARS-CoV-2 saltational evolution. We develop a scalable framework to detect saltational events and apply it to 4.4 million high-quality genomes. We find that saltational events occurred at a low but detectable rate throughout the pandemic and across geographies. They harbor a distinct mutational signature strongly resembling the one  observed in persistent infections, supporting the role of prolonged human infections in their emergence. While most saltational events don’t show evidence of onward transmission, those that do tend to carry mutations also found in successful clades. Our work suggests that the emergence of highly mutated SARS-CoV-2 variants reflects a persistent evolutionary process, with implications for epidemic preparedness. 


## Repository organization
This repository is organized around a Snakemake workflow that enables to process and analyze the viridian data, download SRA metadata, run the Bayesian model we developed, identify adaptive branches, visualize the adaptive branches in the context of the global SARS-CoV-2 phylogeny and investigate mutation patterns in adaptive branches. 

This repository is organized in sub-folders as follows:
- ```data/``` stores input data.  Further information is available on the folder-level [README](data/README.md) file.
- ```figures/``` contains the figures and the scripts to generate the figures (both from the main text and the supplementary information) associated with the manuscript. Further information is available on the folder-level [README](figures/README.md) file.
- ```manuscript/``` contains the manuscript.
- ```scripts/``` contains the custom scripts used in the workflow. Further information is available on the folder-level [README](scripts/README.md) file.
- ```workflow/``` contains the Snakemake files and environment files used in the main Snakefile.

## Running the workflow

### 0. Requirements
- [Snakemake](https://snakemake.readthedocs.io/) with conda support
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (for environment management)

### 1. Downloading required input data

The workflow requires the following (large) input files for the **SARS-CoV-2 viridian phylogeny**. The files correspond to a large mutation-annotated tree from 4.4 million high-quality SARS-CoV-2 samples (described in [Hunt et al.](https://doi.org/10.1038/s41592-025-02947-1)). Hunt et al. provide a nice detailed description of available files from their processing pipeline. This large tree is also publicly available on Taxonium at [this link](https://www.taxonium.org/sars-cov-2/viridian).

To be able to rerun the workflow, these files should be downloaded from https://doi.org/10.6084/m9.figshare.27194547 for both the `.pb` and `.jsonl` format. For the workflow to run properly, these two files should be put at the following paths: `data/viridian/tree.all_viridian.202409.pb.gz` and `data/viridian/tree.all_viridian.202409.jsonl.gz`.

### 2. Running the worfklow
The workflow can be run using:
```bash
snakemake --use-conda --conda-frontend conda --cores <number_of_cores>
```

## Reproducing the figures

The workflow generates files in the `results/` folder which are subsequently used to perform descriptive analysis and figures (code in the `scripts/` folder). Our version of the `results` folder can be downloaded at **xxx (Add link)**.

