# ==============================================================================
# Optional runner for the final figure scripts
# ==============================================================================
# Run this file from the repository root after placing all required inputs in data/.
# Each script also runs independently and is usually easier to debug separately.

scripts_to_run <- c(
  "figure_1_ttests_vs_positive_control.R",
  "figure_2B_3B_fpkm_venn.R",
  "figure_2C_3C_marker_heatmaps.R",
  "figure_2D_3D_stage_associated_heatmaps.R",
  "figure_2E_3E_go_enrichment.R",
  "figure_3F_secreted_functions_stacked_barplot.R",
  "figure_4A_pca.R",
  "figure_4B_4C_4D_core_signatures.R",
  "figure_5_ttests_vs_positive_control.R"
)

for (script_name in scripts_to_run) {
  script_path <- file.path("scripts", script_name)
  if (!file.exists(script_path)) {
    stop("Script not found: ", script_path)
  }
  message("Running ", script_path)
  source(script_path, echo = FALSE)
}

message("All final figure scripts completed.")
