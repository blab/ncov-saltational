# 

## Setup & installation
We use conda to define an environment containing all the required depencencies to run the Snakemake pipeline.
We can then use conda to install the required dependencies.
The environment can be created and activated using the following commands:

```bash
# Install
conda env create -f envs/bte.yaml
# Activate (not necessary here as the environment will be automatically activated within the Snakemake workflow)
source activate bte
```

If needed, the environment can be deactivated / deleted using the following command:
```bash
# Deactivate 
conda deactivate
# Remove
conda env remove --name bte
```

The `Snakefile` can be run using:

```bash
snakemake --use-conda --conda-frontend conda --cores 1
```

## Repository organisation

