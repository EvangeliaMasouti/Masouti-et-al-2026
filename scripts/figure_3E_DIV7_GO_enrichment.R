# Figure 3E: self-contained DIV7 cultured GO enrichment
# This reproduces the Figure 2E-style workflow and does not source common.R.
# Gene sets are upregulated-only unions across the two target-versus-other-stage comparisons.

# =========================================================
# Figure 3E
# GO enrichment of cultured DIV7 stage-specific genes
#
# Exact Figure 2E-style workflow:
#   - cultured DIV7 samples only
#   - union of significant upregulated genes from either comparison
#   - padj < 0.05
#   - log2FoldChange > 0.5
#   - SYMBOL to ENTREZID mapping
#   - compareCluster() with enrichGO(), BP ontology
#   - top 10 GO terms per stage in the plot
#   - Excel sheets corresponding directly to the plot
#
# This is the upregulated-only Figure 2E definition.
# Downregulated genes are not included in this plot or workbook.
# =========================================================

rm(list = ls(all.names = TRUE))
graphics.off()

# =========================================================
# 1. LOAD LIBRARIES
# =========================================================

library(readxl)
library(writexl)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

# =========================================================
# 2. PATHS AND SETTINGS
# =========================================================

project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "results", "Figure_3E_DIV7_GO")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "gene.count.all.xlsx")
input_sheet <- 1

output_png <- file.path(output_dir, "Figure3E_GO_Dotplot.png")
output_pdf <- file.path(output_dir, "Figure3E_GO_Dotplot.pdf")
output_excel <- file.path(output_dir, "Figure3E_GO_Results.xlsx")

PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.5
GO_ONTOLOGY <- "BP"
GO_PVALUE_CUTOFF <- 0.05
TOP_TERMS_PER_STAGE <- 10

if (!file.exists(input_file)) {
  stop(
    input_file
  )
}

# =========================================================
# 3. LOAD CULTURED DIV7 COUNT DATA
# =========================================================

# Reading the first worksheet matches the original Figure 2E script.
df <- readxl::read_excel(
  input_file,
  sheet = input_sheet,
  .name_repair = "unique"
)

sample_columns <- c(
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

required_columns <- c(
  "gene_name",
  sample_columns
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    "The following cultured-sample columns are missing:\n",
    paste(missing_columns, collapse = ", ")
  )
}

counts_data <- df %>%
  dplyr::select(
    gene_name,
    dplyr::all_of(sample_columns)
  ) %>%
  dplyr::mutate(
    gene_name = as.character(gene_name),
    dplyr::across(
      dplyr::all_of(sample_columns),
      as.numeric
    )
  )

if (anyNA(counts_data[, sample_columns])) {
  stop(
    "Missing or non-numeric values were detected in the cultured count matrix."
  )
}

if (any(as.matrix(counts_data[, sample_columns]) < 0)) {
  stop(
    "Negative counts were detected in the cultured count matrix."
  )
}

# =========================================================
# 4. PREPARE THE COUNT MATRIX
# =========================================================

# This matches the original Figure 2E workflow:
# duplicate symbols receive unique DESeq2 row names and are not summed.

gene_symbol_original <- trimws(
  as.character(counts_data$gene_name)
)

missing_symbols <- is.na(gene_symbol_original) |
  gene_symbol_original == ""

if (any(missing_symbols)) {
  stop(
    "Missing gene symbols were found. Resolve them before running GO analysis."
  )
}

deseq2_row_names <- make.unique(
  gene_symbol_original
)

gene_lookup <- data.frame(
  DESeq2_row = deseq2_row_names,
  SYMBOL = gene_symbol_original,
  stringsAsFactors = FALSE
)

counts_matrix <- as.matrix(
  counts_data[, sample_columns, drop = FALSE]
)

counts_matrix <- round(
  counts_matrix
)

storage.mode(counts_matrix) <- "integer"
rownames(counts_matrix) <- deseq2_row_names

counts_matrix <- counts_matrix[
  rowSums(counts_matrix) > 0,
  ,
  drop = FALSE
]

if (anyDuplicated(rownames(counts_matrix))) {
  stop("Duplicate DESeq2 row names remain.")
}

# =========================================================
# 5. PREPARE SAMPLE METADATA AND RUN DESEQ2
# =========================================================

sample_names <- colnames(
  counts_matrix
)

conditions_vec <- dplyr::case_when(
  grepl("^p3div7", sample_names) ~ "p3div7",
  grepl("^p10div7", sample_names) ~ "p10div7",
  grepl("^p14div7", sample_names) ~ "p14div7",
  TRUE ~ NA_character_
)

if (anyNA(conditions_vec)) {
  stop(
    "At least one sample could not be assigned to a cultured DIV7 stage."
  )
}

coldata <- data.frame(
  condition = factor(
    conditions_vec,
    levels = c(
      "p3div7",
      "p10div7",
      "p14div7"
    )
  ),
  row.names = sample_names,
  check.names = FALSE
)

if (!identical(rownames(coldata), colnames(counts_matrix))) {
  stop(
    "The sample metadata and count matrix are not aligned."
  )
}

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- dds[
  rowSums(DESeq2::counts(dds)) > 0,
  ,
  drop = FALSE
]

dds <- DESeq2::DESeq(
  dds,
  quiet = TRUE
)

# =========================================================
# 6. ORIGINAL UNION-BASED UPREGULATED GENE SETS
# =========================================================

get_total_up <- function(
  target,
  ref1,
  ref2
) {

  res1 <- DESeq2::results(
    dds,
    contrast = c(
      "condition",
      target,
      ref1
    ),
    alpha = PADJ_CUTOFF
  )

  res2 <- DESeq2::results(
    dds,
    contrast = c(
      "condition",
      target,
      ref2
    ),
    alpha = PADJ_CUTOFF
  )

  genes_ref1 <- rownames(
    res1[
      which(
        !is.na(res1$padj) &
          !is.na(res1$log2FoldChange) &
          res1$padj < PADJ_CUTOFF &
          res1$log2FoldChange > LOG2FC_CUTOFF
      ),
      ,
      drop = FALSE
    ]
  )

  genes_ref2 <- rownames(
    res2[
      which(
        !is.na(res2$padj) &
          !is.na(res2$log2FoldChange) &
          res2$padj < PADJ_CUTOFF &
          res2$log2FoldChange > LOG2FC_CUTOFF
      ),
      ,
      drop = FALSE
    ]
  )

  # Original union rule from Figure 2E
  unique(
    c(
      genes_ref1,
      genes_ref2
    )
  )
}

p3_cultured_total <- get_total_up(
  "p3div7",
  "p10div7",
  "p14div7"
)

p10_cultured_total <- get_total_up(
  "p10div7",
  "p3div7",
  "p14div7"
)

p14_cultured_total <- get_total_up(
  "p14div7",
  "p3div7",
  "p10div7"
)

union_gene_sets <- list(
  "p3 cultured" = p3_cultured_total,
  "p10 cultured" = p10_cultured_total,
  "p14 cultured" = p14_cultured_total
)

message("Union upregulated gene counts before mapping:")
for (stage_name in names(union_gene_sets)) {
  message(
    stage_name,
    ": ",
    length(union_gene_sets[[stage_name]])
  )
}

# =========================================================
# 7. SYMBOL TO ENTREZID MAPPING
# =========================================================

get_original_symbols <- function(
  deseq2_rows
) {

  symbols <- gene_lookup$SYMBOL[
    match(
      deseq2_rows,
      gene_lookup$DESeq2_row
    )
  ]

  unique(
    symbols[
      !is.na(symbols) &
        symbols != ""
    ]
  )
}

get_entrez_table <- function(
  deseq2_rows,
  stage_name
) {

  if (length(deseq2_rows) == 0) {
    return(data.frame(
      Condition = character(),
      DESeq2_row = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      stringsAsFactors = FALSE
    ))
  }

  symbols <- get_original_symbols(
    deseq2_rows
  )

  mapped <- suppressWarnings(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db,
      drop = TRUE
    )
  )

  if (is.null(mapped) || nrow(mapped) == 0) {
    return(data.frame(
      Condition = character(),
      DESeq2_row = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      stringsAsFactors = FALSE
    ))
  }

  mapped <- mapped %>%
    dplyr::mutate(
      SYMBOL = as.character(SYMBOL),
      ENTREZID = as.character(ENTREZID),
      Condition = stage_name
    ) %>%
    dplyr::left_join(
      gene_lookup,
      by = "SYMBOL"
    ) %>%
    dplyr::select(
      Condition,
      DESeq2_row,
      SYMBOL,
      ENTREZID
    ) %>%
    dplyr::distinct()

  mapped
}

mapped_entrez_table <- dplyr::bind_rows(
  lapply(
    names(union_gene_sets),
    function(stage_name) {
      get_entrez_table(
        union_gene_sets[[stage_name]],
        stage_name
      )
    }
  )
)

mapped_gene_lists <- lapply(
  names(union_gene_sets),
  function(stage_name) {

    mapped_entrez_table %>%
      dplyr::filter(
        Condition == stage_name
      ) %>%
      dplyr::pull(ENTREZID) %>%
      unique()
  }
)

names(mapped_gene_lists) <- names(union_gene_sets)

mapped_gene_lists <- mapped_gene_lists[
  lengths(mapped_gene_lists) > 0
]

if (length(mapped_gene_lists) == 0) {
  stop(
    "No union upregulated genes could be mapped to Entrez IDs."
  )
}

# =========================================================
# 8. RUN compareCluster() GO ENRICHMENT
# =========================================================

ck <- clusterProfiler::compareCluster(
  geneCluster = mapped_gene_lists,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = GO_ONTOLOGY,
  pAdjustMethod = "BH",
  pvalueCutoff = GO_PVALUE_CUTOFF,
  readable = TRUE
)

df_go <- as.data.frame(
  ck
)

if (nrow(df_go) == 0) {
  stop(
    "compareCluster() returned no significant GO terms."
  )
}

# =========================================================
# 9. REBUILD THE FIGURE 2E-STYLE PLOT FOR FIGURE 3E
# =========================================================

# Select the exact rows used in the plot, stage by stage.
# The GO ID is retained so that identical descriptions from different stages
# are not accidentally mixed between the plot and the Excel workbook.
selected_plot_rows <- df_go %>%
  dplyr::filter(
    !is.na(p.adjust)
  ) %>%
  dplyr::group_by(
    Cluster
  ) %>%
  dplyr::arrange(
    p.adjust,
    ID,
    Description,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = TOP_TERMS_PER_STAGE
  ) %>%
  dplyr::mutate(
    Stage_Top_Rank = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

if (!nrow(selected_plot_rows)) {
  stop("No GO terms were selected for the Figure 3E plot.")
}

plot_data <- selected_plot_rows %>%
  dplyr::mutate(
    Cluster = factor(
      Cluster,
      levels = c(
        "p3 cultured",
        "p10 cultured",
        "p14 cultured"
      )
    ),
    Gene_Count = as.numeric(
      Count
    ),
    padj = as.numeric(
      p.adjust
    )
  )

term_order <- plot_data %>%
  dplyr::group_by(
    Description
  ) %>%
  dplyr::summarize(
    min_padj = min(padj),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    desc(min_padj)
  ) %>%
  dplyr::pull(
    Description
  )

plot_data$Description <- factor(
  plot_data$Description,
  levels = term_order
)

p_final <- ggplot(
  plot_data,
  aes(
    x = Cluster,
    y = Description
  )
) +
  geom_point(
    aes(
      size = Gene_Count,
      color = padj
    )
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  scale_color_gradientn(
    colors = c(
      "#8A2BE2",
      "#C71585",
      "#E04880",
      "#F28E2B",
      "#F1CE63"
    ),
    limits = c(
      0,
      GO_PVALUE_CUTOFF
    ),
    name = "adj. p-value"
  ) +
  scale_size_continuous(
    range = c(
      2.5,
      8.5
    ),
    name = "Gene Count"
  ) +
  labs(
    title = "GO enrichment of cultured stage-specific genes",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(
      size = 11,
      color = "black",
      family = "sans"
    ),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      face = "bold",
      family = "sans"
    ),
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.box = "vertical",
    panel.grid.major.x = element_line(
      color = "grey92"
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

print(
  p_final
)

ggsave(
  output_png,
  p_final,
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  output_pdf,
  p_final,
  width = 10,
  height = 12,
  device = cairo_pdf
)

# =========================================================
# 10. PREPARE EXACT EXCEL TABLES
# =========================================================

# Full GO table
all_go_terms <- df_go %>%
  dplyr::select(
    Cluster,
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) %>%
  dplyr::rename(
    Condition = Cluster,
    GO_ID = ID,
    Term_Description = Description,
    adj_pvalue = p.adjust,
    Gene_Count = Count,
    Gene_Symbols = geneID
  )

# Exact GO rows represented in the plot, matched by stage and GO ID.
selected_keys <- selected_plot_rows %>%
  dplyr::transmute(
    Condition = as.character(Cluster),
    GO_ID = as.character(ID)
  ) %>%
  dplyr::distinct()

plot_terms <- all_go_terms %>%
  dplyr::semi_join(
    selected_keys,
    by = c("Condition", "GO_ID")
  )

plot_data_export <- plot_data %>%
  dplyr::mutate(
    Cluster = as.character(Cluster),
    Description = as.character(Description)
  ) %>%
  dplyr::select(
    Cluster,
    Description,
    Gene_Count,
    padj
  )

input_gene_table <- dplyr::bind_rows(
  lapply(
    names(union_gene_sets),
    function(stage_name) {
      rows <- union_gene_sets[[stage_name]]
      data.frame(
        Condition = stage_name,
        DESeq2_row = rows,
        Original_SYMBOL = gene_lookup$SYMBOL[
          match(rows, gene_lookup$DESeq2_row)
        ],
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::distinct()

p3_go <- all_go_terms %>%
  dplyr::filter(
    Condition == "p3 cultured"
  )

p10_go <- all_go_terms %>%
  dplyr::filter(
    Condition == "p10 cultured"
  )

p14_go <- all_go_terms %>%
  dplyr::filter(
    Condition == "p14 cultured"
  )

analysis_settings <- data.frame(
  Parameter = c(
    "Input file",
    "Input worksheet",
    "Samples",
    "Analysis direction",
    "Stage-specific gene rule",
    "Adjusted p-value cutoff",
    "log2FoldChange cutoff",
    "GO ontology",
    "GO p-value cutoff",
    "Top GO terms per stage in plot",
    "Gene identifier used for DESeq2",
    "Gene identifier used for GO mapping",
    "Plot file",
    "Excel file"
  ),
  Value = c(
    input_file,
    "First worksheet",
    paste(sample_columns, collapse = ", "),
    "Upregulated only",
    "Union of genes significant against either comparator",
    PADJ_CUTOFF,
    LOG2FC_CUTOFF,
    GO_ONTOLOGY,
    GO_PVALUE_CUTOFF,
    TOP_TERMS_PER_STAGE,
    "make.unique(gene_name)",
    "Original SYMBOL to ENTREZID",
    output_png,
    output_excel
  ),
  stringsAsFactors = FALSE
)

validation_checks <- data.frame(
  Check = c(
    "Input file exists",
    "Nine cultured DIV7 samples present",
    "Three P3 replicates",
    "Three P10 replicates",
    "Three P14 replicates",
    "No missing count values",
    "No negative counts",
    "No duplicate DESeq2 row names",
    "All selected plot rows occur in full GO table",
    "All plot clusters occur in input gene sets",
    "Excel plot-term rows match plot data rows"
  ),
  Result = c(
    file.exists(input_file),
    length(sample_columns) == 9,
    sum(grepl("^p3div7", sample_columns)) == 3,
    sum(grepl("^p10div7", sample_columns)) == 3,
    sum(grepl("^p14div7", sample_columns)) == 3,
    !anyNA(counts_matrix),
    !any(counts_matrix < 0),
    !anyDuplicated(rownames(counts_matrix)),
    all(paste(plot_terms$Condition, plot_terms$GO_ID) %in%
          paste(all_go_terms$Condition, all_go_terms$GO_ID)),
    all(unique(plot_data_export$Cluster) %in% names(union_gene_sets)),
    nrow(plot_terms) == nrow(plot_data_export)
  ),
  stringsAsFactors = FALSE
)

# =========================================================
# 11. EXPORT THE PI-READY EXCEL WORKBOOK
# =========================================================

writexl::write_xlsx(
  list(
    All_GO_Terms = all_go_terms,
    Selected_Plot_Rows = selected_plot_rows,
    Top_Terms_In_Plot = plot_terms,
    Plot_Data = plot_data_export,
    P3_Cultured_GO = p3_go,
    P10_Cultured_GO = p10_go,
    P14_Cultured_GO = p14_go,
    Union_Input_Genes = input_gene_table,
    Mapped_Entrez_IDs = mapped_entrez_table,
    Analysis_Settings = analysis_settings,
    Validation_Checks = validation_checks
  ),
  path = output_excel
)

# =========================================================
# 12. COMPLETION MESSAGE
# =========================================================

message(
  "Done! Figure 3E union-based cultured plot and Excel workbook exported."
)

cat(
  "\nFiles saved in:\n",
  getwd(),
  "\n\nPlot:\n",
  file.path(getwd(), output_png),
  "\n\nExcel workbook:\n",
  file.path(getwd(), output_excel),
  "\n"
)

# =========================================================
# End of script
# =========================================================


writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Figure 3E outputs were saved in: ", normalizePath(output_dir))
