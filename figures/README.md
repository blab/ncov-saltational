# Scripts to generate figures and analyze model results

## Scripts to reproduce figures
### Figure 1
- [illustrate_bayesian_model.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/illustrate_bayesian_model.R) contains the script to reproduce Figure 1B and 1C.
- [summary_mut_per_branch_type_saltational.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/summary_mut_per_branch_type_saltational.R) contains the script to reproduce Figure 1D and 1E.

### Figure 2
- [plot_saltational_branches_through_time.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_saltational_branches_through_time.R) contains the script to reproduce Figure 2A and 2B.
- [plot_saltational_branches_through_space.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_saltational_branches_through_space.R) contains the script to reproduce Figure 2C, 2D and 2E.

### Figure 3
- [plot_AA_pos_enrichment_saltational.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/) contains the script to reproduce Figure 3A, 3B, 3C and 3D.
- [plot_OR_AA_pos_antigenic_spike.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_OR_AA_pos_antigenic_spike.R) contains the script to reproduce Figure 3E.
- [plot_enriched_mut_on_spike_structure.pml](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_enriched_mut_on_spike_structure.pml) contains the script to reproduce Figure 3F (to be run with Pymol).

### Figure 4
- [full_figure_chronic.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/full_figure_chronic.R) contains the script to reproduce Figure 4A-C.

### Figure 5
- [plot_saltational_branches_with_descendants.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_saltational_branches_with_descendants.R) 
contains the script to reproduce Figure 5A and 5B.
- [plot_fold_change_nb_desc_nb_muts_by_gene.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_fold_change_nb_desc_nb_muts_by_gene.R) contains the script to reproduce Figure 5C.
- [plot_impact_mutation_descendant_proba.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_impact_mutation_descendant_proba.R) contains the script to reproduce Figure 5D and 5E.

### Supplementary figures
- [plot_dist_nonsyn_mut_count.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_dist_nonsyn_mut_count.R) contains the script to reproduce Figure S1.
- [plot_saltational_branches_through_time_sensitivity.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_saltational_branches_through_time_sensitivity.R) contains the script to reproduce Figure S2.
- [plot_saltational_branches_through_space.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_saltational_branches_through_space.R) contains the script to reproduce Figure S3.
- [plot_timeline_cluster_from_saltational_branches.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_timeline_cluster_from_saltational_branches.R) contains the script to reproduce Figure S4.
- [impact_thresholds_number_branches.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/impact_thresholds_number_branches.R) contains the script to reproduce Figure S5.
- [plot_posterior_predictive_checks.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_posterior_predictive_checks.R) contains the script to reproduce Figure S6, S7 and S8.
- [compare_OR_AA_mut_accounting_background.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/compare_OR_AA_mut_accounting_background.R) contains the script to reproduce Figure S9.
- [plot_fold_change_nb_desc_nb_muts_by_gene.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/plot_fold_change_nb_desc_nb_muts_by_gene.R) contains the script to reproduce Figure S10.

## Utility scripts
The following scripts are sourced by the figure scripts above and do not produce figures directly:
- [config_plot.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/config_plot.R): color palettes and shared visual configuration used by most figure scripts.
- [utils_plot_mutation_profile.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/utils_plot_mutation_profile.R): plot helper functions for Figure 3.
- [utils_plot_chronic.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/utils_plot_chronic.R): plot helper functions for Figure 4.
- [utils_descendants_analysis.R](https://github.com/blab/ncov-adapt-branches/blob/main/figures/scripts/utils_descendants_analysis.R): analysis and plot helper functions for Figure 5.
