# 

## Setup & installation
We use conda to define an environment containing all the required depencencies to run the Snakemake pipeline.

The `Snakefile` can be run using:

```bash
snakemake --use-conda --conda-frontend conda --cores 1
```

### NCBI API Key Setup

The pipeline queries NCBI databases (BioSample and SRA) for metadata extraction. It is possible to add the API key to one's configuration by editing `config/config.yaml` and add your API key:
```yaml
ncbi_api_key: "your_api_key_here"
```

**Note**: The NCBI EDirect tools use the `NCBI_API_KEY` environment variable internally. The pipeline automatically sets this from the config file.

## Repository organisation

