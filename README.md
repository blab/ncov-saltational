# Code and data for: SARS-CoV-2 saltational events are recurrent and trace to persistent human infections

## Abstract

SARS-CoV-2 evolution is characterized by gradual mutation accumulation but has been punctuated by rare yet impactful highly mutated variants. Whether such saltational jumps are a broad feature of SARS-CoV-2 evolution or rare anomalies remains unclear. We systematically investigate SARS-CoV-2 saltational evolution by developing a scalable framework to detect saltational events from 4.4 million high-quality viral genomes. Saltational events occurred at low but detectable rates during the pandemic and post-pandemic periods and across geographies. Their mutational signature closely matches that seen in persistent human infections but is inconsistent with the signatures of mink or deer infections. This points to persistent infection, rather than reverse zoonosis, as their primary source. While most saltational events lack evidence of onward transmission, those that do tend to carry mutations found in successful clades. Our work demonstrates that the emergence of highly mutated SARS-CoV-2 variants reflects a recurrent evolutionary process, with implications for preparedness. 


## Repository organization
This repository is organized around a Snakemake workflow that enables to process and analyze the viridian data, download SRA metadata, run the Bayesian model we developed, identify adaptive branches, visualize the adaptive branches in the context of the global SARS-CoV-2 phylogeny and investigate mutation patterns in adaptive branches. 

This repository is organized in sub-folders as follows:
- ```data/``` stores input data.  Further information is available on the folder-level [README](data/README.md) file.
- ```figures/``` contains the figures and the scripts to generate the figures (both from the main text and the supplementary information) associated with the manuscript. Further information is available on the folder-level [README](figures/README.md) file.
- ```manuscript/``` contains the manuscript.
- ```scripts/``` contains the custom scripts used in the workflow. Further information is available on the folder-level [README](scripts/README.md) file.
- ```workflow/``` contains the Snakemake files and environment files used in the main Snakefile.

## Running the workflow

### 0. System requirements

**Operating system**: Tested on macOS 26.5.2.

**Required software** (versions used during development):
- [Snakemake](https://snakemake.readthedocs.io/) with conda support (tested with v8.11.1)
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (for environment management)
- R (tested with v4.3.1)

All other dependencies (Python, bioinformatics tools, and R packages) are managed automatically via conda environments defined in `workflow/envs/`.

**Typical installation time**: ~5 minutes for conda environment setup.

### 1. Downloading required input data

The workflow requires the following (large) input files for the **SARS-CoV-2 viridian phylogeny**. The files correspond to a large mutation-annotated tree from 4.4 million high-quality SARS-CoV-2 samples (described in [Hunt et al.](https://doi.org/10.1038/s41592-025-02947-1)). Hunt et al. provide a nice detailed description of available files from their processing pipeline. This large tree is also publicly available on Taxonium at [this link](https://www.taxonium.org/sars-cov-2/viridian).

To be able to rerun the workflow, these files should be downloaded from https://doi.org/10.6084/m9.figshare.27194547 for both the `.pb` and `.jsonl` format. For the workflow to run properly, these two files should be put at the following paths: `data/viridian/tree.all_viridian.202409.pb.gz` and `data/viridian/tree.all_viridian.202409.jsonl.gz`.

This can be done using the command line with:
```bash
curl -L -o data/viridian/tree.all_viridian.202409.pb.gz "https://ndownloader.figshare.com/files/49691037"
curl -L -o data/viridian/tree.all_viridian.202409.jsonl.gz "https://ndownloader.figshare.com/files/49691040"
```

Downloading metadata from the SRA required some manual curation to collate both metadata directly available on the SRA but also adding metadata that we couldn't download but that were abailable in the Viridian metadata. This results in a complicated and somewhat messy workflow that we don't recommend rerunning. The resulting metadata files are available with the other workflow output on Figshare [here](https://doi.org/10.6084/m9.figshare.32736657).

### 2. Running the worfklow
The workflow can be run using:
```bash
snakemake --use-conda --conda-frontend conda --cores <number_of_cores>
```

## Reproducing the figures

The workflow generates files in the `results/` folder which are subsequently used to perform descriptive analysis and figures (code in the `figures/scripts/` folder). Our version of the `results` folder can be downloaded [here](https://doi.org/10.6084/m9.figshare.32736657).

To reproduce all manuscript figures:

1. Download the pre-computed results from [Figshare](https://doi.org/10.6084/m9.figshare.32736657) and place in `results/`
2. Run the figure scripts in `figures/scripts/` using R v4.3.1:
   ```bash
   Rscript figures/scripts/<script_name>.R
   ```

