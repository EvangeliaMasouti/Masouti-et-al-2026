# ==============================================================================
# Figures 2D and 3D
# Heatmaps of stage-associated upregulated genes
#
# Differential-expression definition:
#   DESeq2 on raw gene counts
#   Each target stage is compared separately with both other stages
#   A gene is stage-associated only when it is significant in BOTH comparisons
#   Adjusted p-value < 0.05 in both comparisons
#   log2 fold-change > 0.5 in both comparisons
#   Upregulated genes only
#   Ranking by mean adjusted p-value across the two comparisons
#
# Heatmap transformation:
#   DESeq2-normalized counts
#   log2(normalized counts + 1)
#   Row-wise z-score
#   Euclidean distance and complete-linkage row clustering
#
# Figure 2D, direct astrocytes:
#   P3 = 50, P10 = 5, P14 = 50, total = 105 genes
#
# Figure 3D, DIV7 cultured astrocytes:
#   P3 = 50, P10 = 50, P14 = 50, total = 150 genes
# ==============================================================================

# Run this script from the project root.
# Expected input: data/gene.count.all.xlsx

# --- 1. Project paths and analysis settings ----------------------------------
input_file <- file.path("data", "gene.count.all.xlsx")
input_sheet <- "gene_count"
output_dir <- file.path("results", "Stage_Associated_Heatmaps_Figures_2D_3D")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.5

# These are display limits, not DEG thresholds.
display_n_direct <- c(
  P3 = 50L,
  P10 = 5L,
  P14 = 50L
)

display_n_cultured <- c(
  P3 = 50L,
  P10 = 50L,
  P14 = 50L
)

# --- 2. Required packages -----------------------------------------------------
required_packages <- c(
  "readxl",
  "writexl",
  "dplyr",
  "tibble",
  "DESeq2",
  "pheatmap"
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
  library(DESeq2)
  library(pheatmap)
})

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# --- 3. Load and validate the raw count matrix -------------------------------
raw <- readxl::read_excel(
  input_file,
  sheet = input_sheet,
  .name_repair = "unique"
)

required_columns <- c(
  "gene_id",
  "gene_name",
  "p3direct1",
  "p3direct2",
  "p3direct3",
  "p10direct1",
  "p10direct2",
  "p10direct3",
  "p14direct1",
  "p14direct2",
  "p14direct3",
  "p3div7_1",
  "p3div7_2",
  "p3div7_3",
  "p10div7_1",
  "p10div7_2",
  "p10div7_3",
  "p14div7_1",
  "p14div7_2",
  "p14div7_3"
)

missing_columns <- setdiff(required_columns, colnames(raw))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

sample_columns <- setdiff(required_columns, c("gene_id", "gene_name"))

raw <- raw %>%
  dplyr::mutate(
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name),
    dplyr::across(
      dplyr::all_of(sample_columns),
      as.numeric
    )
  )

if (anyNA(raw$gene_id) || any(raw$gene_id == "")) {
  stop("Missing gene_id values were detected.")
}

if (anyDuplicated(raw$gene_id)) {
  stop("Duplicate gene_id values were detected.")
}

if (anyNA(as.matrix(raw[, sample_columns]))) {
  stop("Missing or non-numeric values were detected in the count matrix.")
}

if (any(as.matrix(raw[, sample_columns]) < 0)) {
  stop("Negative count values were detected.")
}

if (any(as.matrix(raw[, sample_columns]) != round(as.matrix(raw[, sample_columns])))) {
  stop("The DESeq2 input contains non-integer count values.")
}

# --- 4. Define the two experimental compartments -----------------------------
compartments <- list(
  Direct = list(
    samples = c(
      "p3direct1", "p3direct2", "p3direct3",
      "p10direct1", "p10direct2", "p10direct3",
      "p14direct1", "p14direct2", "p14direct3"
    ),
    sample_stage = c(
      p3direct1 = "P3",
      p3direct2 = "P3",
      p3direct3 = "P3",
      p10direct1 = "P10",
      p10direct2 = "P10",
      p10direct3 = "P10",
      p14direct1 = "P14",
      p14direct2 = "P14",
      p14direct3 = "P14"
    ),
    display_n = display_n_direct,
    figure_label = "Figure 2D",
    output_stem = "Figure_2D_Direct_Stage_Associated_Heatmap",
    title = "Figure 2D: Stage-associated genes in directly isolated astrocytes"
  ),
  Cultured_DIV7 = list(
    samples = c(
      "p3div7_1", "p3div7_2", "p3div7_3",
      "p10div7_1", "p10div7_2", "p10div7_3",
      "p14div7_1", "p14div7_2", "p14div7_3"
    ),
    sample_stage = c(
      p3div7_1 = "P3",
      p3div7_2 = "P3",
      p3div7_3 = "P3",
      p10div7_1 = "P10",
      p10div7_2 = "P10",
      p10div7_3 = "P10",
      p14div7_1 = "P14",
      p14div7_2 = "P14",
      p14div7_3 = "P14"
    ),
    display_n = display_n_cultured,
    figure_label = "Figure 3D",
    output_stem = "Figure_3D_Cultured_DIV7_Stage_Associated_Heatmap",
    title = "Figure 3D: Stage-associated genes in DIV7 cultured astrocytes"
  )
)

# --- 5. Helper functions ------------------------------------------------------
run_deseq2 <- function(data, samples, sample_stage) {
  count_data <- data %>%
    dplyr::select(
      gene_id,
      dplyr::all_of(samples)
    ) %>%
    as.data.frame()

  rownames(count_data) <- count_data$gene_id
  count_data <- count_data[, samples, drop = FALSE]
  count_data <- as.matrix(count_data)
  storage.mode(count_data) <- "integer"

  col_data <- data.frame(
    Age = factor(
      unname(sample_stage[samples]),
      levels = c("P3", "P10", "P14")
    ),
    row.names = samples
  )

  # Remove genes with zero counts in all samples of this compartment.
  keep <- rowSums(count_data) > 0
  count_data <- count_data[keep, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ Age
  )

  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  list(
    dds = dds,
    normalized_counts = DESeq2::counts(dds, normalized = TRUE)
  )
}

get_target_statistics <- function(
  dds,
  target_stage,
  reference_stages,
  data
) {
  result_1 <- DESeq2::results(
    dds,
    contrast = c("Age", target_stage, reference_stages[1]),
    alpha = PADJ_CUTOFF
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::transmute(
      gene_id,
      log2FC_vs_reference_1 = log2FoldChange,
      pvalue_vs_reference_1 = pvalue,
      padj_vs_reference_1 = padj,
      baseMean_vs_reference_1 = baseMean
    )

  result_2 <- DESeq2::results(
    dds,
    contrast = c("Age", target_stage, reference_stages[2]),
    alpha = PADJ_CUTOFF
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::transmute(
      gene_id,
      log2FC_vs_reference_2 = log2FoldChange,
      pvalue_vs_reference_2 = pvalue,
      padj_vs_reference_2 = padj,
      baseMean_vs_reference_2 = baseMean
    )

  statistics <- result_1 %>%
    dplyr::inner_join(result_2, by = "gene_id") %>%
    dplyr::left_join(
      data %>%
        dplyr::select(gene_id, gene_name) %>%
        dplyr::distinct(),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      Target_Stage = target_stage,
      Direction = "Upregulated",
      mean_adjusted_p = rowMeans(
        cbind(
          padj_vs_reference_1,
          padj_vs_reference_2
        ),
        na.rm = FALSE
      ),
      mean_log2FC = rowMeans(
        cbind(
          log2FC_vs_reference_1,
          log2FC_vs_reference_2
        ),
        na.rm = FALSE
      ),
      eligible_for_display =
        !is.na(padj_vs_reference_1) &
        !is.na(padj_vs_reference_2) &
        padj_vs_reference_1 < PADJ_CUTOFF &
        padj_vs_reference_2 < PADJ_CUTOFF &
        log2FC_vs_reference_1 > LOG2FC_CUTOFF &
        log2FC_vs_reference_2 > LOG2FC_CUTOFF
    ) %>%
    dplyr::filter(eligible_for_display) %>%
    dplyr::arrange(
      mean_adjusted_p,
      dplyr::desc(mean_log2FC),
      gene_id
    )

  statistics
}

select_display_genes <- function(
  statistics,
  target_stage,
  n_display,
  figure_label
) {
  if (nrow(statistics) < n_display) {
    stop(
      figure_label,
      ": fewer than ",
      n_display,
      " eligible genes were found for ",
      target_stage,
      ". Found ",
      nrow(statistics),
      ". Check the input counts and DESeq2 settings."
    )
  }

  statistics %>%
    dplyr::slice_head(n = n_display) %>%
    dplyr::mutate(
      Display_Rank = dplyr::row_number(),
      Requested_Display_N = n_display,
      Figure = figure_label,
      Stage_Group = paste0(target_stage, "_Unique")
    )
}

make_row_zscore <- function(normalized_counts, gene_ids) {
  log2_counts <- log2(
    normalized_counts[gene_ids, , drop = FALSE] + 1
  )

  z_matrix <- t(scale(t(log2_counts)))
  z_matrix[is.na(z_matrix)] <- 0
  z_matrix
}

make_display_gene_names <- function(statistics) {
  display_names <- statistics$gene_name
  display_names[is.na(display_names) | display_names == ""] <- statistics$gene_id[
    is.na(display_names) | display_names == ""
  ]

  # Keep gene symbols readable while ensuring unique row names.
  make.unique(display_names)
}

# --- 6. Run one complete compartment analysis --------------------------------
run_compartment <- function(data, compartment_name, settings) {
  samples <- settings$samples
  sample_stage <- settings$sample_stage
  stage_levels <- c("P3", "P10", "P14")

  deseq <- run_deseq2(
    data = data,
    samples = samples,
    sample_stage = sample_stage
  )

  all_statistics <- list()
  selected_statistics <- list()

  for (target_stage in stage_levels) {
    reference_stages <- setdiff(stage_levels, target_stage)

    stage_statistics <- get_target_statistics(
      dds = deseq$dds,
      target_stage = target_stage,
      reference_stages = reference_stages,
      data = data
    )

    stage_selected <- select_display_genes(
      statistics = stage_statistics,
      target_stage = target_stage,
      n_display = unname(settings$display_n[target_stage]),
      figure_label = settings$figure_label
    )

    all_statistics[[target_stage]] <- stage_statistics %>%
      dplyr::mutate(
        Figure = settings$figure_label,
        Stage_Group = paste0(target_stage, "_Unique")
      )

    selected_statistics[[target_stage]] <- stage_selected
  }

  selected_statistics_table <- dplyr::bind_rows(selected_statistics)
  selected_gene_ids <- selected_statistics_table$gene_id

  # Order selected genes by target stage, then by mean adjusted p-value.
  selected_statistics_table <- selected_statistics_table %>%
    dplyr::mutate(
      Stage_Group = factor(
        Stage_Group,
        levels = paste0(stage_levels, "_Unique")
      )
    ) %>%
    dplyr::arrange(Stage_Group, Display_Rank) %>%
    dplyr::mutate(Stage_Group = as.character(Stage_Group))

  selected_gene_ids <- selected_statistics_table$gene_id

  z_matrix <- make_row_zscore(
    normalized_counts = deseq$normalized_counts,
    gene_ids = selected_gene_ids
  )

  display_names <- make_display_gene_names(selected_statistics_table)
  rownames(z_matrix) <- display_names

  row_annotation <- selected_statistics_table %>%
    dplyr::select(Stage_Group) %>%
    as.data.frame()
  rownames(row_annotation) <- display_names

  # The requested clustering is explicitly Euclidean, complete linkage.
  row_hclust <- hclust(
    dist(z_matrix, method = "euclidean"),
    method = "complete"
  )

  stage_colors <- if (compartment_name == "Direct") {
    c(
      "P3_Unique" = "#8c96c6",
      "P10_Unique" = "#b3cde3",
      "P14_Unique" = "#88419d"
    )
  } else {
    c(
      "P3_Unique" = "#8c96c6",
      "P10_Unique" = "#b3cde3",
      "P14_Unique" = "#88419d"
    )
  }

  annotation_colors <- list(
    Stage_Group = stage_colors
  )

  # Keep samples ordered by developmental stage and replicate.
  ordered_samples <- samples[order(
    factor(unname(sample_stage[samples]), levels = stage_levels),
    samples
  )]
  z_matrix <- z_matrix[, ordered_samples, drop = FALSE]

  row_annotation <- row_annotation[rownames(z_matrix), , drop = FALSE]

  heatmap_pdf <- file.path(
    output_dir,
    paste0(settings$output_stem, ".pdf")
  )
  heatmap_png <- file.path(
    output_dir,
    paste0(settings$output_stem, ".png")
  )

  draw_heatmap <- function() {
    pheatmap::pheatmap(
      z_matrix,
      color = colorRampPalette(
        c("#313695", "#74add1", "#f7f7f7", "#f46d43", "#a50026")
      )(100),
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      clustering_distance_rows = "euclidean",
      clustering_method = "complete",
      annotation_row = row_annotation,
      annotation_colors = annotation_colors,
      border_color = NA,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = ifelse(nrow(z_matrix) > 120, 4, 5),
      fontsize_col = 8,
      angle_col = 45,
      main = paste0(
        settings$title,
        "\nRow-wise z-score of log2(DESeq2-normalized counts + 1)"
      )
    )
  }

  pdf(heatmap_pdf, width = 8, height = 12)
  draw_heatmap()
  dev.off()

  png(heatmap_png, width = 2400, height = 3600, res = 300)
  draw_heatmap()
  dev.off()

  # Add heatmap values and row order to the audit output.
  heatmap_data <- as.data.frame(z_matrix) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::mutate(
      Stage_Group = row_annotation$Stage_Group[match(Gene, rownames(row_annotation))],
      .before = 2
    )

  selected_statistics_table <- selected_statistics_table %>%
    dplyr::mutate(
      Heatmap_Row_Order = match(
        gene_id,
        selected_statistics_table$gene_id[row_hclust$order]
      ),
      Included_in_displayed_heatmap = TRUE
    )

  all_statistics_table <- dplyr::bind_rows(all_statistics) %>%
    dplyr::mutate(
      Included_in_displayed_heatmap = gene_id %in% selected_gene_ids
    )

  list(
    dds = deseq$dds,
    normalized_counts = deseq$normalized_counts,
    heatmap_matrix = z_matrix,
    heatmap_data = heatmap_data,
    selected_statistics = selected_statistics_table,
    all_statistics = all_statistics_table,
    row_annotation = row_annotation,
    row_hclust = row_hclust,
    heatmap_pdf = heatmap_pdf,
    heatmap_png = heatmap_png
  )
}

# --- 7. Run direct and cultured analyses -------------------------------------
direct_result <- run_compartment(
  data = raw,
  compartment_name = "Direct",
  settings = compartments$Direct
)

cultured_result <- run_compartment(
  data = raw,
  compartment_name = "Cultured_DIV7",
  settings = compartments$Cultured_DIV7
)

# Explicit final-size checks for the confirmed displayed figures.
if (nrow(direct_result$selected_statistics) != 105) {
  stop("Figure 2D must contain exactly 105 displayed genes.")
}

if (nrow(cultured_result$selected_statistics) != 150) {
  stop("Figure 3D must contain exactly 150 displayed genes.")
}

# --- 8. Save audit workbooks --------------------------------------------------
write_compartment_audit <- function(result, compartment_name, file_stem) {
  output_file <- file.path(
    output_dir,
    paste0(file_stem, ".xlsx")
  )

  writexl::write_xlsx(
    list(
      Heatmap_Data = result$heatmap_data,
      Selected_Gene_Statistics = result$selected_statistics,
      All_Eligible_Gene_Statistics = result$all_statistics
    ),
    path = output_file
  )
}

write_compartment_audit(
  direct_result,
  compartment_name = "Direct",
  file_stem = "Figure_2D_Direct_Heatmap_Audit"
)

write_compartment_audit(
  cultured_result,
  compartment_name = "Cultured_DIV7",
  file_stem = "Figure_3D_Cultured_DIV7_Heatmap_Audit"
)

analysis_settings <- tibble::tibble(
  Parameter = c(
    "Input file",
    "Input sheet",
    "Differential-expression method",
    "Stage comparison",
    "Adjusted p-value cutoff",
    "Log2 fold-change cutoff",
    "Direction",
    "Ranking variable",
    "Heatmap transformation",
    "Row clustering",
    "Direct display",
    "Cultured display"
  ),
  Value = c(
    basename(input_file),
    input_sheet,
    "DESeq2 using raw counts",
    "Each target stage versus both other stages; significant in both comparisons",
    "padj < 0.05 in both comparisons",
    "log2FC > 0.5 in both comparisons",
    "Upregulated genes only",
    "Mean adjusted p-value across the two target-versus-other-stage comparisons",
    "log2(DESeq2-normalized counts + 1), followed by row-wise z-score",
    "Euclidean distance with complete-linkage clustering",
    "P3 = 50; P10 = 5; P14 = 50; total = 105",
    "P3 = 50; P10 = 50; P14 = 50; total = 150"
  )
)

writexl::write_xlsx(
  list(
    Analysis_Settings = analysis_settings,
    Direct_Selected_Genes = direct_result$selected_statistics,
    Cultured_Selected_Genes = cultured_result$selected_statistics
  ),
  path = file.path(
    output_dir,
    "Figures_2D_3D_Stage_Associated_Heatmaps_Audit.xlsx"
  )
)

cat(
  "Figures 2D and 3D were generated using the specified DESeq2 stage-associated-gene workflow.\n",
  "Direct display: 105 genes (P3 = 50, P10 = 5, P14 = 50).\n",
  "Cultured display: 150 genes (P3 = 50, P10 = 50, P14 = 50).\n",
  "Outputs saved in: ", normalizePath(output_dir), "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================
