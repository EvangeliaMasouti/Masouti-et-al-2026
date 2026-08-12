# ==============================================================================
# Figure 3D
# Top unique stage-associated genes in DIV7 cultured astrocytes
# ==============================================================================

# ==============================================================================
# Expected project structure
# ==============================================================================
#
# project/
# ├── data/
# │   └── gene.count.all.xlsx
# ├── scripts/
# │   └── figure_3D_top_unique_stage_genes.R
# └── results/
#
# ==============================================================================


# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================

project_dir <- getwd()

data_dir <- file.path(
  project_dir,
  "data"
)

results_dir <- file.path(
  project_dir,
  "results",
  "figure_3D"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_file <- file.path(
  data_dir,
  "gene.count.all.xlsx"
)

if (!file.exists(input_file)) {
  stop(
    "The input file was not found:\n",
    input_file
  )
}

library(readxl)
library(writexl)
library(DESeq2)
library(dplyr)
library(pheatmap)


# ==============================================================================
# 2. LOAD COMPLETE COUNT MATRIX
# ==============================================================================

df <- read_excel(
  input_file,
  .name_repair = "unique"
)

if (!"gene_name" %in% colnames(df)) {
  stop(
    "The input file must contain a column named 'gene_name'."
  )
}


# ==============================================================================
# 3. IDENTIFY AND ORDER ALL SAMPLE COLUMNS
# ==============================================================================

condition_patterns <- c(
  p3direct  = "^p3direct",
  p10direct = "^p10direct",
  p14direct = "^p14direct",
  p3div7    = "^p3div7",
  p10div7   = "^p10div7",
  p14div7   = "^p14div7"
)

set_order <- names(
  condition_patterns
)

sample_cols <- unlist(
  lapply(
    set_order,
    function(condition) {
      grep(
        condition_patterns[[condition]],
        colnames(df),
        value = TRUE
      )
    }
  )
)

sample_cols <- unique(
  sample_cols
)

if (length(sample_cols) == 0) {
  stop(
    "No sample columns were detected in the input file."
  )
}

counts_data <- df %>%
  dplyr::select(
    gene_name,
    all_of(sample_cols)
  ) %>%
  mutate(
    across(
      all_of(sample_cols),
      as.numeric
    )
  )

if (anyNA(counts_data[, sample_cols])) {
  stop(
    "Missing or non-numeric values were detected in the count matrix."
  )
}

if (any(counts_data[, sample_cols] < 0, na.rm = TRUE)) {
  stop(
    "Negative values were detected in the count matrix."
  )
}


# ==============================================================================
# 4. AGGREGATE DUPLICATE GENE ROWS
# ==============================================================================

gene_ids <- as.character(
  counts_data$gene_name
)

missing_gene_ids <- is.na(gene_ids) |
  trimws(gene_ids) == ""

gene_ids[missing_gene_ids] <- paste0(
  "unnamed_gene_",
  which(missing_gene_ids)
)

gene_ids <- trimws(
  gene_ids
)

# Duplicate gene rows are summed before DESeq2 analysis.
counts_df <- counts_data %>%
  mutate(
    Gene = gene_ids
  ) %>%
  dplyr::select(
    Gene,
    all_of(sample_cols)
  ) %>%
  group_by(
    Gene
  ) %>%
  summarise(
    across(
      all_of(sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

counts_matrix <- as.matrix(
  counts_df[, sample_cols, drop = FALSE]
)

rownames(counts_matrix) <- counts_df$Gene

counts_matrix <- round(
  counts_matrix
)

storage.mode(counts_matrix) <- "integer"

# Remove genes with no counts in any sample
counts_matrix <- counts_matrix[
  rowSums(counts_matrix) > 0,
  ,
  drop = FALSE
]


# ==============================================================================
# 5. PREPARE SAMPLE METADATA
# ==============================================================================

sample_names <- colnames(
  counts_matrix
)

condition_code <- vapply(
  sample_names,
  function(sample_name) {
    
    matches <- names(condition_patterns)[
      vapply(
        condition_patterns,
        function(pattern) {
          grepl(
            pattern,
            sample_name
          )
        },
        logical(1)
      )
    ]
    
    if (length(matches) != 1) {
      stop(
        "Sample could not be assigned to exactly one condition: ",
        sample_name
      )
    }
    
    matches
  },
  character(1)
)

coldata <- data.frame(
  condition = factor(
    condition_code,
    levels = set_order
  ),
  row.names = sample_names,
  check.names = FALSE
)

ordered_sample_names <- unlist(
  lapply(
    set_order,
    function(condition) {
      sample_names[
        condition_code == condition
      ]
    }
  )
)

counts_matrix <- counts_matrix[
  ,
  ordered_sample_names,
  drop = FALSE
]

coldata <- coldata[
  ordered_sample_names,
  ,
  drop = FALSE
]


# ==============================================================================
# 6. RUN DESEQ2
# ==============================================================================

cat("Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(
  dds,
  quiet = TRUE
)


# ==============================================================================
# 7. IDENTIFY TOP UNIQUE DIV7 STAGE-ASSOCIATED GENES
# ==============================================================================

padj_threshold <- 0.05
log2fc_threshold <- 0.5
top_n <- 50

get_top_unique_div7 <- function(
    target,
    reference_1,
    reference_2,
    top_n = 50
) {
  
  result_reference_1 <- results(
    dds,
    contrast = c(
      "condition",
      target,
      reference_1
    )
  )
  
  result_reference_2 <- results(
    dds,
    contrast = c(
      "condition",
      target,
      reference_2
    )
  )
  
  genes_reference_1 <- rownames(
    result_reference_1
  )[
    !is.na(result_reference_1$padj) &
      result_reference_1$padj < padj_threshold &
      result_reference_1$log2FoldChange > log2fc_threshold
  ]
  
  genes_reference_2 <- rownames(
    result_reference_2
  )[
    !is.na(result_reference_2$padj) &
      result_reference_2$padj < padj_threshold &
      result_reference_2$log2FoldChange > log2fc_threshold
  ]
  
  # A unique stage-associated gene must be significant against
  # both other developmental stages.
  significant_genes <- intersect(
    genes_reference_1,
    genes_reference_2
  )
  
  if (length(significant_genes) == 0) {
    return(character(0))
  }
  
  average_adjusted_p <- (
    result_reference_1[
      significant_genes,
      "padj"
    ] +
      result_reference_2[
        significant_genes,
        "padj"
      ]
  ) / 2
  
  names(average_adjusted_p) <- significant_genes
  
  ranked_genes <- names(
    sort(
      average_adjusted_p,
      decreasing = FALSE,
      na.last = NA
    )
  )
  
  head(
    ranked_genes,
    top_n
  )
}

p3_div7_unique_top <- get_top_unique_div7(
  target = "p3div7",
  reference_1 = "p10div7",
  reference_2 = "p14div7",
  top_n = top_n
)

p10_div7_unique_top <- get_top_unique_div7(
  target = "p10div7",
  reference_1 = "p3div7",
  reference_2 = "p14div7",
  top_n = top_n
)

p14_div7_unique_top <- get_top_unique_div7(
  target = "p14div7",
  reference_1 = "p3div7",
  reference_2 = "p10div7",
  top_n = top_n
)

display_genes_div7 <- unique(
  c(
    p3_div7_unique_top,
    p10_div7_unique_top,
    p14_div7_unique_top
  )
)

display_genes_div7 <- display_genes_div7[
  !is.na(display_genes_div7) &
    display_genes_div7 != ""
]

if (length(display_genes_div7) == 0) {
  stop(
    "No DIV7 stage-associated genes passed the specified thresholds."
  )
}


# ==============================================================================
# 8. VST TRANSFORMATION
# ==============================================================================

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vsd
)


# ==============================================================================
# 9. EXTRACT DIV7 CONDITIONS FOR THE HEATMAP
# ==============================================================================

div7_cols <- rownames(
  coldata
)[
  coldata$condition %in% c(
    "p3div7",
    "p10div7",
    "p14div7"
  )
]

heatmap_annotation_div7 <- coldata[
  div7_cols,
  ,
  drop = FALSE
]

mat_unq_div7 <- vst_matrix[
  display_genes_div7,
  div7_cols,
  drop = FALSE
]


# ==============================================================================
# 10. ROW-WISE Z-SCORE SCALING
# ==============================================================================

row_sd_div7 <- apply(
  mat_unq_div7,
  1,
  sd,
  na.rm = TRUE
)

valid_genes_div7 <- is.finite(row_sd_div7) &
  row_sd_div7 > 0

mat_unq_div7 <- mat_unq_div7[
  valid_genes_div7,
  ,
  drop = FALSE
]

display_genes_div7 <- rownames(
  mat_unq_div7
)

if (nrow(mat_unq_div7) == 0) {
  stop(
    "All selected DIV7 genes had zero or non-finite variance."
  )
}

mat_unq_scaled_div7 <- t(
  scale(
    t(mat_unq_div7)
  )
)

rownames(mat_unq_scaled_div7) <- display_genes_div7


# ==============================================================================
# 11. PREPARE GENE-GROUP ANNOTATIONS
# ==============================================================================

gene_groups_div7 <- data.frame(
  Gene = display_genes_div7,
  Stage_Group = case_when(
    display_genes_div7 %in% p3_div7_unique_top ~
      "P3_DIV7_Unique",
    
    display_genes_div7 %in% p10_div7_unique_top ~
      "P10_DIV7_Unique",
    
    display_genes_div7 %in% p14_div7_unique_top ~
      "P14_DIV7_Unique",
    
    TRUE ~ "Other"
  ),
  check.names = FALSE
)

rownames(gene_groups_div7) <- gene_groups_div7$Gene

row_annotation_div7 <- gene_groups_div7 %>%
  dplyr::select(
    Stage_Group
  )

row_annotation_div7$Stage_Group <- factor(
  row_annotation_div7$Stage_Group,
  levels = c(
    "P3_DIV7_Unique",
    "P10_DIV7_Unique",
    "P14_DIV7_Unique",
    "Other"
  )
)

annotation_colours_div7 <- list(
  condition = c(
    "p3direct" = "lightcyan",
    "p10direct" = "orchid",
    "p14direct" = "pink",
    "p3div7" = "#008B8B",
    "p10div7" = "#4B0082",
    "p14div7" = "#C71585"
  ),
  Stage_Group = c(
    "P3_DIV7_Unique" = "#008B8B",
    "P10_DIV7_Unique" = "#4B0082",
    "P14_DIV7_Unique" = "#C71585",
    "Other" = "grey80"
  )
)


# ==============================================================================
# 12. GENERATE FIGURE 3D
# ==============================================================================

heatmap_colours <- colorRampPalette(
  c(
    "#000033",
    "blue",
    "cyan",
    "green",
    "yellow"
  )
)(100)

heatmap_object_div7 <- pheatmap(
  mat_unq_scaled_div7,
  annotation_col = heatmap_annotation_div7,
  annotation_row = row_annotation_div7,
  annotation_colors = annotation_colours_div7,
  show_rownames = TRUE,
  fontsize_row = 6,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = heatmap_colours,
  border_color = NA,
  main = "Top Unique Stage Genes (DIV7 Cultured, n=150)",
  silent = TRUE
)

pdf_file <- file.path(
  results_dir,
  "Figure_3D_Top_Unique_Stage_Genes_Heatmap_DIV7.pdf"
)

png_file <- file.path(
  results_dir,
  "Figure_3D_Top_Unique_Stage_Genes_Heatmap_DIV7.png"
)

pdf(
  pdf_file,
  width = 8,
  height = 14,
  useDingbats = FALSE
)

grid::grid.newpage()
grid::grid.draw(
  heatmap_object_div7$gtable
)

dev.off()

png(
  png_file,
  width = 2400,
  height = 4200,
  res = 300
)

grid::grid.newpage()
grid::grid.draw(
  heatmap_object_div7$gtable
)

dev.off()


# ==============================================================================
# 13. PREPARE EXPORT TABLES
# ==============================================================================

df_scaled_div7 <- data.frame(
  Gene = rownames(mat_unq_scaled_div7),
  mat_unq_scaled_div7,
  check.names = FALSE
)

df_vst_div7 <- data.frame(
  Gene = rownames(mat_unq_div7),
  mat_unq_div7,
  check.names = FALSE
)

df_p3_div7 <- data.frame(
  Rank = seq_along(p3_div7_unique_top),
  Gene = p3_div7_unique_top,
  check.names = FALSE
)

df_p10_div7 <- data.frame(
  Rank = seq_along(p10_div7_unique_top),
  Gene = p10_div7_unique_top,
  check.names = FALSE
)

df_p14_div7 <- data.frame(
  Rank = seq_along(p14_div7_unique_top),
  Gene = p14_div7_unique_top,
  check.names = FALSE
)

analysis_settings <- data.frame(
  Parameter = c(
    "Adjusted p-value threshold",
    "log2 fold-change threshold",
    "Top genes per stage",
    "Number of genes displayed"
  ),
  Value = c(
    padj_threshold,
    log2fc_threshold,
    top_n,
    nrow(mat_unq_scaled_div7)
  ),
  check.names = FALSE
)


# ==============================================================================
# 14. EXPORT EXCEL WORKBOOK
# ==============================================================================

excel_file <- file.path(
  results_dir,
  "Figure_3D_Top_Unique_Stage_Genes_Data_DIV7.xlsx"
)

write_xlsx(
  list(
    "Heatmap_Scaled_Values" = df_scaled_div7,
    "Normalized_Counts_VST" = df_vst_div7,
    "Gene_Group_Metadata" = gene_groups_div7,
    "P3_DIV7_Unique_Top50" = df_p3_div7,
    "P10_DIV7_Unique_Top50" = df_p10_div7,
    "P14_DIV7_Unique_Top50" = df_p14_div7,
    "Sample_Metadata" = data.frame(
      Sample = rownames(heatmap_annotation_div7),
      heatmap_annotation_div7,
      check.names = FALSE
    ),
    "Analysis_Settings" = analysis_settings
  ),
  path = excel_file
)


# ==============================================================================
# 15. SAVE SESSION INFORMATION
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "sessionInfo.txt"
  )
)


# ==============================================================================
# 16. COMPLETION MESSAGE
# ==============================================================================

cat(
  "\nFigure 3D analysis completed successfully.\n",
  "PDF: ", pdf_file, "\n",
  "PNG: ", png_file, "\n",
  "Excel data: ", excel_file, "\n",
  "Number of displayed genes: ",
  nrow(mat_unq_scaled_div7),
  "\n",
  sep = ""
)
