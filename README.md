# Script inventory

Run scripts from the repository root. Each script is independent unless noted.

| Script | Analysis |
|---|---|
| `run_all_figures.R` | Optional sequential runner |
| `figure_1_ttests_vs_positive_control.R` | Figure 1 marker-percentage t-tests |
| `figure_2B_3B_fpkm_venn.R` | Shared FPKM > 1 expression-set Venn implementation |
| `figure_3B_DIV7_fpkm_venn.R` | Figure 3B wrapper |
| `figure_2C_3C_marker_heatmaps.R` | Shared canonical marker heatmap implementation, including Atp1b2 |
| `figure_3C_DIV7_marker_heatmap.R` | Figure 3C wrapper |
| `figure_2D_3D_stage_associated_heatmaps.R` | Shared DESeq2 stage-associated heatmap implementation |
| `figure_3D_DIV7_stage_associated_heatmap.R` | Figure 3D wrapper |
| `figure_2E_3E_go_enrichment.R` | Shared GO enrichment implementation |
| `figure_3E_DIV7_GO_enrichment.R` | Figure 3E wrapper |
| `figure_3F_secreted_functions_stacked_barplot.R` | Secreted-function stacked barplot |
| `figure_4A_pca.R` | PCA plot |
| `figure_4B_4C_4D_core_signatures.R` | Core-signature volcano and GO plots |
| `figure_5_ttests_vs_positive_control.R` | Metabolite-experiment marker t-tests |
```

Run the scripts from the repository root. The Figure 2E and 3E GO scripts are self-contained and do not require `common.R`.
