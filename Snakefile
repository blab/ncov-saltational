"""
Adaptive branch identification and analysis pipeline for SARS-CoV-2.

This workflow identifies adaptive branches in the SARS-CoV-2 phylogeny using Bayesian
modeling and performs proximity-based subsampling of background sequences.

Workflow stages:
    1. viridian: Tree preparation and temporal inference
    2. sra_metadata: Extract sample metadata from NCBI SRA/BioSample databases
    3. mutations: Mutation analysis and tree topology
    4. bayesian_model: Selection analysis with Bayesian modeling
    5. adaptive_branches: Identification and preparation of adaptive branches
    6. mutation_enrichment: Enrichment analysis of mutations in adaptive branches
"""

# Load configuration
configfile: "config/config.yaml"

# Global configuration
GENES = ['ORF1a', 'ORF1b', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b', 'Genome']
CV_PRIORS = [0.95, 1.5, 2.0, 3.0, 5.0]
BF_THRESHOLDS = [31.62277660]

rule all:
    input:
        # Viridian analysis
        df_nuc_mut = "results/mutations/df_nuc_mut.tsv",
        df_aa_muts_coded = "results/mutations/nodes_aa_muts_summary_coded.tsv",
        df_sites = "results/bte/df_sites_syn_non_syn_truncated_at_stop.tsv",
        df_parent_relationships = "results/bte/df_parent_child_relationship.tsv",
        df_n_descendants = "results/bte/df_n_descendants_per_node.tsv",
        chronumental = "results/chronumental/chronumental_inferred_dates.tsv",

        # Download metadata from SRA and BioSample
        biosample_metadata = "results/sra_metadata/biosample_metadata.tsv",

        # Bayesian model run
        data_viridian = expand("results/input_model/data_viridian_{gene}.tsv", gene=GENES),
        data_viridian_grouped = expand("results/input_model/data_viridian_{gene}_grouped.tsv", gene=GENES),
        df_prior_parametrization = "results/input_model/df_prior_parametrization.tsv",
        #posterior_simul = expand("results/posterior_predictive_checks/posterior_predictive_check_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),
        #df_param_estim = expand("results/output_model/df_param_estim_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),
        #df_BF = expand("results/output_model/df_BF_{gene}_CV_{cv_prior}.tsv", gene=GENES, cv_prior = CV_PRIORS),

        # Adaptive branches characterization
        #name_adaptive_branches = expand("results/adaptive_branches/name_adaptive_branches_{gene}_CV_{cv_prior}_BF_{BF_threshold}.tsv", gene=GENES, cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #concat_name_adaptive_branches = expand("results/adaptive_branches/concat_name_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

        # 
        #ancestral_fasta = expand("results/alignment_adaptive_branches/alignment_adaptive_internal_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #alignment_adaptive_final_nodes = expand("results/alignment_adaptive_branches/alignment_adaptive_final_nodes_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #alignment_all_adaptive_branches = expand("results/alignment_adaptive_branches/alignment_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #metadata = expand("results/alignment_adaptive_branches/metadata_all_adaptive_branches_CV_{cv_prior}_BF_{BF_threshold}.tsv", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        
        #proximal_samples = expand("results/proximity/proximal_samples_all_CV_{cv_prior}_BF_{BF_threshold}.txt", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #alignment_proximal = expand("results/proximity/alignment_proximal_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #temporal_alignment = expand("results/proximity/temporal_background_sequences_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #tree_sequences = expand("results/tree_context/sequences_CV_{cv_prior}_BF_{BF_threshold}.fasta", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

        # Phylogenetic tree building
        #auspice_trees = expand("results/tree/auspice_tree_CV_{cv_prior}_BF_{BF_threshold}.json", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

        # Enrichment analysis
        #df_AA_pos_enrichment = expand("results/mutation_enrichment/df_AA_pos_enrichment_CV_{cv_prior}_BF_{BF_threshold}.tsv", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),
        #df_AA_mut_enrichment = expand("results/mutation_enrichment/df_AA_mut_enrichment_CV_{cv_prior}_BF_{BF_threshold}.tsv", cv_prior = CV_PRIORS, BF_threshold = BF_THRESHOLDS),

rule params:
    input:
        usher_reference_gtf = "params/reference/reference_seq_all.gtf",
        usher_reference_fasta = "params/reference/reference_seq.fasta",
        viridian_tree = "data/viridian/tree.all_viridian.202409.pb.gz",
        viridian_taxonium = "data/viridian/tree.all_viridian.202409.jsonl.gz"
    params:
        ncbi_sra_url = "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab",
        sra_metadata_download_batch_size = 1000

# Include workflow modules
include: "workflow/rules/viridian_analysis.smk"
include: "workflow/rules/sra_metadata.smk"
include: "workflow/rules/bayesian_model.smk"
include: "workflow/rules/adaptive_branches.smk"
include: "workflow/rules/mutation_enrichment.smk"
include: "workflow/rules/tree_building.smk"
