# ==============================================================================
# Figure 4B, 4C and 4D
# Core transcriptional signatures distinguishing cultured and direct astrocytes
#
# Figure 4B: P10 cultured-versus-direct volcano plot with core signatures marked
# Figure 4C: GO enrichment of the core direct signature
# Figure 4D: GO enrichment of the core cultured signature
#
# Core cultured signature:
#   Genes upregulated in cultured versus direct astrocytes at P3, P10 and P14.
#
# Core direct signature:
#   Genes upregulated in direct versus cultured astrocytes at P3, P10 and P14.
#
# The three stage-matched contrasts are analyzed separately using DESeq2.
# ==============================================================================

rm(list = ls(all.names = TRUE))
graphics.off()

# Run from the project root.
# Expected input: data/gene.count.all.xlsx
project_dir <- getwd()
input_file <- file.path(
  project_dir,
  "data",
  "gene.count.all.xlsx"
)
input_sheet <- "gene_count"
output_dir <- file.path(
  project_dir,
  "results",
  "Figure_4B_4C_4D_Core_Signatures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# --- 1. Load libraries --------------------------------------------------------
required_packages <- c(
  "readxl",
  "writexl",
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "dplyr",
  "ggplot2",
  "ggrepel"
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
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# --- 2. Analysis settings -----------------------------------------------------
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.58
GO_ONTOLOGY <- "BP"
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.20

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# --- 3. Load and validate raw counts -----------------------------------------
df <- readxl::read_excel(
  input_file,
  sheet = input_sheet,
  .name_repair = "unique"
)

sample_columns <- c(
  "p3div7_1", "p3div7_2", "p3div7_3",
  "p10div7_1", "p10div7_2", "p10div7_3",
  "p14div7_1", "p14div7_2", "p14div7_3",
  "p3direct1", "p3direct2", "p3direct3",
  "p10direct1", "p10direct2", "p10direct3",
  "p14direct1", "p14direct2", "p14direct3"
)

required_columns <- c(
  "gene_id",
  "gene_name",
  sample_columns
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

raw <- df %>%
  dplyr::select(
    gene_id,
    gene_name,
    dplyr::all_of(sample_columns)
  ) %>%
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
  stop("Missing or non-numeric count values were detected.")
}

if (any(as.matrix(raw[, sample_columns]) < 0)) {
  stop("Negative count values were detected.")
}

if (any(as.matrix(raw[, sample_columns]) != round(as.matrix(raw[, sample_columns])))) {
  stop("Non-integer count values were detected in the DESeq2 input.")
}

# Use gene_id for DESeq2 and retain gene_name only for annotation and labels.
count_matrix <- raw %>%
  dplyr::select(
    gene_id,
    dplyr::all_of(sample_columns)
  ) %>%
  as.data.frame()

rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[, sample_columns, drop = FALSE]
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

# --- 4. Metadata and DESeq2 model ---------------------------------------------
coldata <- data.frame(
  sample = sample_columns,
  group = dplyr::case_when(
    grepl("^p3div7", sample_columns) ~ "P3Cultured",
    grepl("^p10div7", sample_columns) ~ "P10Cultured",
    grepl("^p14div7", sample_columns) ~ "P14Cultured",
    grepl("^p3direct", sample_columns) ~ "P3Direct",
    grepl("^p10direct", sample_columns) ~ "P10Direct",
    grepl("^p14direct", sample_columns) ~ "P14Direct",
    TRUE ~ NA_character_
  ),
  row.names = sample_columns,
  stringsAsFactors = FALSE
)

if (anyNA(coldata$group)) {
  stop("One or more samples could not be assigned to a DESeq2 group.")
}

coldata$group <- factor(
  coldata$group,
  levels = c(
    "P3Direct",
    "P3Cultured",
    "P10Direct",
    "P10Cultured",
    "P14Direct",
    "P14Cultured"
  )
)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ group
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

# --- 5. Extract stage-matched direct-versus-cultured contrasts ----------------
extract_contrast <- function(
  dds,
  cultured_group,
  direct_group
) {
  DESeq2::results(
    dds,
    contrast = c("group", cultured_group, direct_group),
    alpha = PADJ_CUTOFF
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(
      raw %>%
        dplyr::select(gene_id, gene_name),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      Contrast = paste(cultured_group, "vs", direct_group)
    )
}

res_p3 <- extract_contrast(
  dds,
  cultured_group = "P3Cultured",
  direct_group = "P3Direct"
)

res_p10 <- extract_contrast(
  dds,
  cultured_group = "P10Cultured",
  direct_group = "P10Direct"
)

res_p14 <- extract_contrast(
  dds,
  cultured_group = "P14Cultured",
  direct_group = "P14Direct"
)

# --- 6. Identify core consensus signatures -----------------------------------
# Core cultured genes are upregulated in cultured versus direct at all stages.
up_cultured_p3 <- res_p3 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange > LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

up_cultured_p10 <- res_p10 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange > LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

up_cultured_p14 <- res_p14 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange > LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

core_cultured_genes <- Reduce(
  intersect,
  list(up_cultured_p3, up_cultured_p10, up_cultured_p14)
)

# Core direct genes are upregulated in direct versus cultured at all stages,
# represented as negative log2FC in the cultured-versus-direct contrasts.
up_direct_p3 <- res_p3 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange < -LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

up_direct_p10 <- res_p10 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange < -LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

up_direct_p14 <- res_p14 %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(log2FoldChange),
    padj < PADJ_CUTOFF,
    log2FoldChange < -LOG2FC_CUTOFF
  ) %>%
  dplyr::pull(gene_id)

core_direct_genes <- Reduce(
  intersect,
  list(up_direct_p3, up_direct_p10, up_direct_p14)
)

if (length(core_cultured_genes) < 5) {
  warning("Fewer than five core cultured genes were detected.")
}

if (length(core_direct_genes) < 5) {
  warning("Fewer than five core direct genes were detected.")
}

core_consensus <- dplyr::bind_rows(
  raw %>%
    dplyr::filter(gene_id %in% core_cultured_genes) %>%
    dplyr::transmute(
      Gene_ID = gene_id,
      Gene_Symbol = gene_name,
      Category = "Core Cultured Signature"
    ),
  raw %>%
    dplyr::filter(gene_id %in% core_direct_genes) %>%
    dplyr::transmute(
      Gene_ID = gene_id,
      Gene_Symbol = gene_name,
      Category = "Core Direct Signature"
    )
)

# --- 7. GO enrichment ---------------------------------------------------------
run_go <- function(gene_ids, group_name) {
  if (length(gene_ids) < 5) {
    return(tibble::tibble())
  }

  mapped <- clusterProfiler::bitr(
    gene_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  ) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  if (!nrow(mapped)) {
    return(tibble::tibble())
  }

  ego <- clusterProfiler::enrichGO(
    gene = mapped$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = GO_ONTOLOGY,
    pAdjustMethod = "BH",
    pvalueCutoff = GO_PVALUE_CUTOFF,
    qvalueCutoff = GO_QVALUE_CUTOFF,
    readable = TRUE
  )

  if (is.null(ego)) {
    return(tibble::tibble())
  }

  as.data.frame(ego) %>%
    dplyr::mutate(Group = group_name)
}

go_core_cultured <- run_go(
  core_cultured_genes,
  "Core Cultured Signature"
)

go_core_direct <- run_go(
  core_direct_genes,
  "Core Direct Signature"
)

all_go_consensus <- dplyr::bind_rows(
  go_core_cultured,
  go_core_direct
)

if (!nrow(all_go_consensus)) {
  stop("No enriched GO terms were found for the core signatures.")
}

# --- 8. GO dot-plot helper ----------------------------------------------------
build_dotplot <- function(go_data, title_text) {
  if (!nrow(go_data)) {
    stop("No GO terms available for: ", title_text)
  }

  top_data <- go_data %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust <= GO_PVALUE_CUTOFF
    ) %>%
    dplyr::slice_min(
      order_by = p.adjust,
      n = 20,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(
      neg_log10_p = -log10(p.adjust)
    )

  if (!nrow(top_data)) {
    stop("No GO terms passed the adjusted p-value cutoff for: ", title_text)
  }

  ggplot(
    top_data,
    aes(
      x = neg_log10_p,
      y = reorder(Description, neg_log10_p)
    )
  ) +
    geom_point(
      aes(
        size = Count,
        color = p.adjust
      )
    ) +
    scale_color_gradient(
      low = "red",
      high = "blue",
      name = "p.adjust"
    ) +
    scale_size_continuous(
      range = c(3, 9),
      name = "Gene Count"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.18))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_text,
      x = expression(-log[10](p.adjust)),
      y = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_line(
        color = "grey85",
        linewidth = 0.5
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(
        size = 11,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 10,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 11,
        face = "bold"
      ),
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      legend.position = "left",
      legend.box = "vertical",
      plot.margin = margin(
        t = 10,
        r = 25,
        b = 10,
        l = 10
      )
    )
}

p_dotplot_direct <- build_dotplot(
  go_core_direct,
  "GO direct astrocytes"
)

p_dotplot_cultured <- build_dotplot(
  go_core_cultured,
  "GO cultured astrocytes"
)

# --- 9. Figure 4B volcano plot -----------------------------------------------
# This volcano plot uses the P10 stage-matched contrast. The colored genes are
# core signatures defined consistently across P3, P10 and P14.
volcano_df <- res_p10 %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    Significance = dplyr::case_when(
      gene_id %in% core_cultured_genes ~ "Core Cultured Gene",
      gene_id %in% core_direct_genes ~ "Core Direct Gene",
      TRUE ~ "Not Core / Not Sig"
    )
  )

volcano_df$Significance <- factor(
  volcano_df$Significance,
  levels = c(
    "Core Cultured Gene",
    "Core Direct Gene",
    "Not Core / Not Sig"
  )
)

top_label_cult <- volcano_df %>%
  dplyr::filter(Significance == "Core Cultured Gene") %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene_name)

top_label_dir <- volcano_df %>%
  dplyr::filter(Significance == "Core Direct Gene") %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene_name)

top_core_labels <- unique(
  c(top_label_cult, top_label_dir)
)

volcano_df <- volcano_df %>%
  dplyr::mutate(
    Label = ifelse(
      gene_name %in% top_core_labels,
      gene_name,
      NA_character_
    )
  )

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(
    data = dplyr::filter(
      volcano_df,
      Significance == "Not Core / Not Sig"
    ),
    color = "grey88",
    alpha = 0.4,
    size = 1.2
  ) +
  geom_point(
    data = dplyr::filter(
      volcano_df,
      Significance != "Not Core / Not Sig"
    ),
    aes(color = Significance),
    alpha = 0.85,
    size = 2
  ) +
  geom_hline(
    yintercept = -log10(PADJ_CUTOFF),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "Core Cultured Gene" = "#4B0082",
      "Core Direct Gene" = "#F3B5FF"
    ),
    drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3.5,
    fontface = "bold.italic",
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 30,
    segment.color = "grey40",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  labs(
    title = "DEG: direct versus cultured astrocytes",
    subtitle = "P10 stage-matched contrast; core signatures are conserved across P3, P10 and P14",
    x = "Log2 Fold Change (P10 Cultured vs P10 Direct)",
    y = "-log10 (Adjusted p-value)",
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(
      size = 13,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 10,
      face = "bold"
    )
  )

# --- 10. Save Figure 4B, 4C and 4D ------------------------------------------
ggsave(
  file.path(output_dir, "Figure4B_Volcano_P10_Cultured_vs_Direct.png"),
  p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure4B_Volcano_P10_Cultured_vs_Direct.pdf"),
  p_volcano,
  width = 8,
  height = 7
)

ggsave(
  file.path(output_dir, "Figure4C_GO_Direct_Core_Signature.png"),
  p_dotplot_direct,
  width = 9.5,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure4C_GO_Direct_Core_Signature.pdf"),
  p_dotplot_direct,
  width = 9.5,
  height = 8
)

ggsave(
  file.path(output_dir, "Figure4D_GO_Cultured_Core_Signature.png"),
  p_dotplot_cultured,
  width = 9.5,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure4D_GO_Cultured_Core_Signature.pdf"),
  p_dotplot_cultured,
  width = 9.5,
  height = 8
)

# --- 11. Save audit tables ----------------------------------------------------
write_xlsx(
  list(
    Core_Consensus_Genes = core_consensus,
    P3_Cultured_vs_Direct = res_p3,
    P10_Cultured_vs_Direct = res_p10,
    P14_Cultured_vs_Direct = res_p14
  ),
  path = file.path(
    output_dir,
    "Figure4B_Core_Consensus_and_DESeq2_Results.xlsx"
  )
)

write_xlsx(
  list(
    Core_Direct_GO = go_core_direct,
    Core_Cultured_GO = go_core_cultured
  ),
  path = file.path(
    output_dir,
    "Figure4C_4D_Core_GO_Results.xlsx"
  )
)

analysis_settings <- data.frame(
  Parameter = c(
    "Input file",
    "Input sheet",
    "DESeq2 design",
    "Stage-matched contrasts",
    "Adjusted p-value cutoff",
    "Absolute log2FC cutoff",
    "Core cultured definition",
    "Core direct definition",
    "Figure 4B volcano contrast",
    "Figure 4C",
    "Figure 4D"
  ),
  Value = c(
    basename(input_file),
    input_sheet,
    "~ group",
    "P3, P10 and P14 cultured versus corresponding direct samples",
    PADJ_CUTOFF,
    LOG2FC_CUTOFF,
    "Upregulated in cultured versus direct at all three stages",
    "Upregulated in direct versus cultured at all three stages",
    "P10 cultured versus P10 direct",
    "GO enrichment of the core direct signature",
    "GO enrichment of the core cultured signature"
  ),
  stringsAsFactors = FALSE
)

write_xlsx(
  list(Analysis_Settings = analysis_settings),
  path = file.path(
    output_dir,
    "Figure4_Analysis_Settings.xlsx"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

cat(
  "Figure 4B, 4C and 4D outputs were saved in: ",
  normalizePath(output_dir),
  "\n",
  "Figure 4B uses the P10 stage-matched volcano contrast.\n",
  "Figures 4C and 4D use the core direct and core cultured signatures, respectively.\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================
