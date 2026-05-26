# SARS-CoV-2 saltational events are recurrent and consistent with evolution during prolonged human infections

## Abstract

SARS-CoV-2 evolution is characterized by the gradual accumulation of mutations but has also been punctuated by the emergence of highly mutated variants. Whether such saltational jumps are a broad feature of SARS-CoV-2 evolution or rare anomalies remains unclear. Here, we perform a systematic analysis of SARS-CoV-2 saltational evolution. We develop a scalable framework to detect saltational events and apply it to 4.4 million high-quality genomes. We find that saltational events occurred at a low but detectable rate throughout the pandemic and across geographies. They harbor a distinct mutational signature strongly resembling the one  observed in persistent infections, supporting the role of prolonged human infections in their emergence. While most saltational events don’t show evidence of onward transmission, those that do tend to carry mutations also found in successful clades. Our work suggests that the emergence of highly mutated SARS-CoV-2 variants reflects a persistent evolutionary process, with implications for epidemic preparedness. 


## Repository organization
This repository is organized around a Snakemake workflow that enables to process and analyze the viridian data, download SRA metadata, run the Bayesian model we developed, identify adaptive branches, visualize the adaptive branches in the context of the global SARS-CoV-2 phylogeny and investigate mutation patterns in adaptive branches. 

This repository is organized in sub-folders as follows:
- ```figures/``` contains the figures and the scripts to generate the figures (both from the main text and the supplementary information) associated with the manuscript. Further information is available on the folder-level [README](figures/README.md) file.
- ```manuscript/``` contains the manuscript.
- ```scripts/``` contains the custom scripts used in the workflow. Further information is available on the folder-level [README](scripts/README.md) file.
- ```workflow/``` contains the Snakemake files and environment files used in the main Snakefile.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) ≥7 with conda support
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (for environment management)
- Python 3.11 and R ≥4.2 (installed automatically via conda environments)

## Getting started

**1. Obtain input data**

The workflow requires two large input files that are not tracked in this repository:

- **Viridian phylogeny** (`data/viridian/tree.all_viridian.202409.pb.gz` and `data/viridian/tree.all_viridian.202409.jsonl.gz`): available from the [Viridian World project](https://viridian.genomics.ox.ac.uk/).
- **Nextstrain global tree** (`data/nextstrain/ncov_open_global_all-time.json`): downloadable from [nextstrain.org](https://nextstrain.org/ncov/open/global/all-time) via the JSON export.

Place these files in the `data/viridian/` and `data/nextstrain/` directories respectively.

**2. Run the workflow**

```bash
snakemake --use-conda -j <number_of_cores>
```

**3. Generate figures**

Figure scripts are located in `figures/scripts/` and must be run from that directory:

```bash
cd figures/scripts
Rscript full_figure_mutational_profile.R
```

See [`figures/README.md`](figures/README.md) for the full list of figure scripts.
