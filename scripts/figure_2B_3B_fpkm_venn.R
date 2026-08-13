# ==============================================================================
# Figures 2B and 3B
# Venn diagrams of genes expressed at each developmental stage
#
# Analysis definition:
#   Original Novogene group-level FPKM values
#   Expression membership: FPKM > 1
#   Venn sets: expressed genes at P3, P10, and P14
#   This is NOT a DESeq2 DEG union or intersection.
#
# Figure 2B: directly isolated astrocytes
# Figure 3B: DIV7 cultured astrocytes
# ==============================================================================

# Run this script from the project root.
# Expected input: data/gene_fpkm_group.xls
# If the file is .xlsx, change the extension below.

# --- 1. Project paths ---------------------------------------------------------
input_file <- file.path("data", "gene_fpkm_group.xls")
output_dir <- file.path("results", "FPKM_Venn_Figures_2B_3B")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# --- 2. Required packages -----------------------------------------------------
required_packages <- c(
  "readxl",
  "writexl",
  "dplyr",
  "tibble",
  "VennDiagram",
  "grid"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tibble)
  library(VennDiagram)
  library(grid)
})

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# --- 3. Load original Novogene group-level FPKM data --------------------------
# read_excel() supports both .xls and .xlsx files.
fpkm <- readxl::read_excel(
  input_file,
  .name_repair = "unique"
)

# These are the six group-level FPKM columns used for the two Venn diagrams.
fpkm_columns <- c(
  "p3direct",
  "p10direct",
  "p14direct",
  "p3div7",
  "p10div7",
  "p14div7"
)

required_columns <- c("gene_id", fpkm_columns)
missing_columns <- setdiff(required_columns, colnames(fpkm))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing from the FPKM file: ",
    paste(missing_columns, collapse = ", ")
  )
}

fpkm <- fpkm %>%
  dplyr::mutate(
    gene_id = as.character(gene_id),
    dplyr::across(
      dplyr::all_of(fpkm_columns),
      as.numeric
    )
  )

if (anyNA(fpkm$gene_id) || any(fpkm$gene_id == "")) {
  stop("Missing gene_id values were detected in the FPKM file.")
}

if (anyDuplicated(fpkm$gene_id)) {
  stop(
    "Duplicate gene_id values were detected. Resolve duplicate gene IDs " ,
    "before generating the Venn diagrams."
  )
}

if (anyNA(as.matrix(fpkm[, fpkm_columns]))) {
  stop("Missing or non-numeric FPKM values were detected.")
}

if (any(as.matrix(fpkm[, fpkm_columns]) < 0)) {
  stop("Negative FPKM values were detected.")
}

# --- 4. Define expression sets ------------------------------------------------
# IMPORTANT: membership is defined by group-level FPKM > 1.
# No fold-change, adjusted p-value, DESeq2, or stage-associated DEG rule is used.
expression_cutoff <- 1

make_expression_set <- function(data, column, cutoff = 1) {
  data$gene_id[data[[column]] > cutoff]
}

direct_sets <- list(
  P3 = make_expression_set(fpkm, "p3direct", expression_cutoff),
  P10 = make_expression_set(fpkm, "p10direct", expression_cutoff),
  P14 = make_expression_set(fpkm, "p14direct", expression_cutoff)
)

cultured_sets <- list(
  P3 = make_expression_set(fpkm, "p3div7", expression_cutoff),
  P10 = make_expression_set(fpkm, "p10div7", expression_cutoff),
  P14 = make_expression_set(fpkm, "p14div7", expression_cutoff)
)

# --- 5. Create membership tables ----------------------------------------------
make_membership_table <- function(
  data,
  stage_sets,
  compartment_name,
  fpkm_columns_for_compartment
) {
  expressed_ids <- sort(unique(unlist(stage_sets)))

  metadata_columns <- intersect(
    c(
      "gene_id",
      "gene_name",
      "gene_chr",
      "gene_start",
      "gene_end",
      "gene_strand",
      "gene_length",
      "gene_biotype",
      "gene_description",
      "tf_family"
    ),
    colnames(data)
  )

  membership <- data %>%
    dplyr::filter(gene_id %in% expressed_ids) %>%
    dplyr::select(
      dplyr::all_of(metadata_columns),
      dplyr::all_of(fpkm_columns_for_compartment)
    ) %>%
    dplyr::mutate(
      P3_expressed = gene_id %in% stage_sets$P3,
      P10_expressed = gene_id %in% stage_sets$P10,
      P14_expressed = gene_id %in% stage_sets$P14,
      Venn_region = dplyr::case_when(
        P3_expressed & !P10_expressed & !P14_expressed ~ "P3 only",
        !P3_expressed & P10_expressed & !P14_expressed ~ "P10 only",
        !P3_expressed & !P10_expressed & P14_expressed ~ "P14 only",
        P3_expressed & P10_expressed & !P14_expressed ~ "P3 and P10 only",
        P3_expressed & !P10_expressed & P14_expressed ~ "P3 and P14 only",
        !P3_expressed & P10_expressed & P14_expressed ~ "P10 and P14 only",
        P3_expressed & P10_expressed & P14_expressed ~ "All three",
        TRUE ~ "Not expressed"
      ),
      Compartment = compartment_name
    ) %>%
    dplyr::arrange(
      factor(
        Venn_region,
        levels = c(
          "P3 only",
          "P10 only",
          "P14 only",
          "P3 and P10 only",
          "P3 and P14 only",
          "P10 and P14 only",
          "All three"
        )
      ),
      gene_id
    )

  membership
}

direct_membership <- make_membership_table(
  data = fpkm,
  stage_sets = direct_sets,
  compartment_name = "Direct",
  fpkm_columns_for_compartment = c(
    "p3direct",
    "p10direct",
    "p14direct"
  )
)

cultured_membership <- make_membership_table(
  data = fpkm,
  stage_sets = cultured_sets,
  compartment_name = "Cultured_DIV7",
  fpkm_columns_for_compartment = c(
    "p3div7",
    "p10div7",
    "p14div7"
  )
)

# --- 6. Create summaries ------------------------------------------------------
venn_region_levels <- c(
  "P3 only",
  "P10 only",
  "P14 only",
  "P3 and P10 only",
  "P3 and P14 only",
  "P10 and P14 only",
  "All three"
)

make_summary <- function(membership, compartment_name, source_columns) {
  stage_summary <- tibble::tibble(
    Compartment = compartment_name,
    Condition = c("P3", "P10", "P14"),
    Source_column = source_columns,
    Threshold_rule = "> 1 FPKM",
    Number_of_genes = c(
      sum(membership$P3_expressed),
      sum(membership$P10_expressed),
      sum(membership$P14_expressed)
    )
  )

  region_summary <- membership %>%
    dplyr::count(
      Venn_region,
      name = "Number_of_genes"
    ) %>%
    dplyr::filter(Venn_region %in% venn_region_levels) %>%
    dplyr::mutate(
      Compartment = compartment_name,
      Condition = Venn_region,
      Source_column = "Venn region",
      Threshold_rule = "> 1 FPKM"
    ) %>%
    dplyr::select(
      Compartment,
      Condition,
      Source_column,
      Threshold_rule,
      Number_of_genes
    )

  dplyr::bind_rows(stage_summary, region_summary)
}

direct_summary <- make_summary(
  direct_membership,
  "Direct",
  c("p3direct", "p10direct", "p14direct")
)

cultured_summary <- make_summary(
  cultured_membership,
  "Cultured_DIV7",
  c("p3div7", "p10div7", "p14div7")
)

summary_table <- dplyr::bind_rows(
  direct_summary,
  cultured_summary
)

analysis_settings <- tibble::tibble(
  Parameter = c(
    "Source file",
    "Source description",
    "Expression cutoff",
    "Set identifier",
    "Direct columns",
    "Cultured columns",
    "Interpretation"
  ),
  Value = c(
    basename(input_file),
    "Original Novogene group-level FPKM table",
    "Group-level FPKM > 1",
    "gene_id",
    paste(c("p3direct", "p10direct", "p14direct"), collapse = ", "),
    paste(c("p3div7", "p10div7", "p14div7"), collapse = ", "),
    "Expression-set membership; not DESeq2 DEGs"
  )
)

# --- 7. Save Venn diagrams ----------------------------------------------------
save_venn_plot <- function(
  stage_sets,
  title,
  output_stem,
  fill_colors
) {
  venn_grob <- VennDiagram::venn.diagram(
    x = stage_sets,
    filename = NULL,
    fill = fill_colors,
    alpha = 0.55,
    cex = 1.2,
    cat.cex = 1.1,
    cat.fontface = "plain",
    fontfamily = "sans",
    cat.fontfamily = "sans",
    margin = 0.08,
    scaled = FALSE,
    euler.d = FALSE,
    ind = TRUE
  )

  pdf(
    file.path(output_dir, paste0(output_stem, ".pdf")),
    width = 7,
    height = 7
  )
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  grid::grid.text(
    title,
    y = grid::unit(0.97, "npc"),
    gp = grid::gpar(
      fontsize = 14,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  dev.off()

  png(
    file.path(output_dir, paste0(output_stem, ".png")),
    width = 2100,
    height = 2100,
    res = 300
  )
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  grid::grid.text(
    title,
    y = grid::unit(0.97, "npc"),
    gp = grid::gpar(
      fontsize = 14,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  dev.off()
}

# Figure 2B: directly isolated astrocytes.
save_venn_plot(
  stage_sets = direct_sets,
  title = "Figure 2B: Genes expressed in directly isolated astrocytes (FPKM > 1)",
  output_stem = "Figure_2B_Direct_FPKM_gt1_Venn",
  fill_colors = c("#8c96c6", "#b3cde3", "#88419d")
)

# Figure 3B: DIV7 cultured astrocytes.
save_venn_plot(
  stage_sets = cultured_sets,
  title = "Figure 3B: Genes expressed in DIV7 cultured astrocytes (FPKM > 1)",
  output_stem = "Figure_3B_Cultured_DIV7_FPKM_gt1_Venn",
  fill_colors = c("#8c96c6", "#b3cde3", "#88419d")
)

# --- 8. Save the complete audit workbook -------------------------------------
output_workbook <- file.path(
  output_dir,
  "Figures_2B_3B_FPKM_Venn_Audit.xlsx"
)

writexl::write_xlsx(
  list(
    Analysis_Settings = analysis_settings,
    Summary = summary_table,
    Direct_Membership = direct_membership,
    Direct_Summary = direct_summary,
    Direct_P3 = direct_membership %>%
      dplyr::filter(P3_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p3direct, Venn_region),
    Direct_P10 = direct_membership %>%
      dplyr::filter(P10_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p10direct, Venn_region),
    Direct_P14 = direct_membership %>%
      dplyr::filter(P14_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p14direct, Venn_region),
    Cultured_Membership = cultured_membership,
    Cultured_Summary = cultured_summary,
    Cultured_P3 = cultured_membership %>%
      dplyr::filter(P3_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p3div7, Venn_region),
    Cultured_P10 = cultured_membership %>%
      dplyr::filter(P10_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p10div7, Venn_region),
    Cultured_P14 = cultured_membership %>%
      dplyr::filter(P14_expressed) %>%
      dplyr::select(gene_id, dplyr::any_of("gene_name"), p14div7, Venn_region)
  ),
  path = output_workbook
)

cat(
  "Figures 2B and 3B were generated from original Novogene group-level FPKM values.\n",
  "Expression membership threshold: FPKM > 1.\n",
  "No DESeq2 DEG rule was used for these Venn diagrams.\n",
  "Outputs saved in: ", normalizePath(output_dir), "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================

# Expected Venn-region counts for the currently supplied audit workbook:
# Direct:   P3 only 552; P10 only 67; P14 only 327;
#           P3 and P10 only 277; P3 and P14 only 95;
#           P10 and P14 only 316; All three 11159.
# Cultured: P3 only 164; P10 only 389; P14 only 285;
#           P3 and P10 only 361; P3 and P14 only 265;
#           P10 and P14 only 213; All three 10967.
# These values are a verification reference, not hard-coded analysis results.
# ==============================================================================
