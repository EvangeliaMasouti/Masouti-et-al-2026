# Figure 2E: self-contained direct GO enrichment
# Uses the original Figure 2E-style union of upregulated genes from either comparator.
# This script does not source common.R.

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
library(purrr)

# =========================================================
# 2. DATA LOADING & DESeq2 ANALYSIS (DIRECT SAMPLES)
# =========================================================

# New working directory: iCloud Drive/Desktop
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "results", "Figure_2E_Direct_GO")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "gene.count.all.xlsx")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

df <- readxl::read_excel(input_file, sheet = 1, .name_repair = "unique")

# Explicit dplyr::select targeting direct samples
counts_data <- df %>%
  dplyr::select(
    gene_name,
    starts_with("p3direct"),
    starts_with("p10direct"),
    starts_with("p14direct")
  ) %>%
  as.data.frame()

rownames(counts_data) <- make.unique(
  as.character(counts_data$gene_name)
)

counts_matrix <- as.matrix(
  counts_data[, -1]
)

sample_names <- colnames(
  counts_matrix
)

conditions_vec <- case_when(
  grepl("p3direct", sample_names) ~ "p3direct",
  grepl("p10direct", sample_names) ~ "p10direct",
  grepl("p14direct", sample_names) ~ "p14direct"
)

coldata <- data.frame(
  condition = factor(
    conditions_vec,
    levels = c(
      "p3direct",
      "p10direct",
      "p14direct"
    )
  )
)

rownames(coldata) <- sample_names

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(
  dds
)

# =========================================================
# 3. ORIGINAL UNION-BASED UPREGULATED GENE SETS
# =========================================================

# Helper function to get unique upregulated genes per stage in DIRECT
get_total_up <- function(
  target,
  ref1,
  ref2
) {

  res1 <- results(
    dds,
    contrast = c(
      "condition",
      target,
      ref1
    )
  )

  res2 <- results(
    dds,
    contrast = c(
      "condition",
      target,
      ref2
    )
  )

  genes_ref1 <- rownames(
    res1[
      which(
        res1$padj < 0.05 &
          res1$log2FoldChange > 0.5
      ),
      ,
      drop = FALSE
    ]
  )

  genes_ref2 <- rownames(
    res2[
      which(
        res2$padj < 0.05 &
          res2$log2FoldChange > 0.5
      ),
      ,
      drop = FALSE
    ]
  )

  # Original union rule
  unique(
    c(
      genes_ref1,
      genes_ref2
    )
  )
}

p3_direct_total <- get_total_up(
  "p3direct",
  "p10direct",
  "p14direct"
)

p10_direct_total <- get_total_up(
  "p10direct",
  "p3direct",
  "p14direct"
)

p14_direct_total <- get_total_up(
  "p14direct",
  "p3direct",
  "p10direct"
)

# =========================================================
# 4. MAP GENES TO ENTREZ IDS
# =========================================================

get_entrez <- function(
  genes
) {

  if (length(genes) == 0) {
    return(NULL)
  }

  ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db"
  )

  ids$ENTREZID
}

gene_list_cluster <- list(
  "p3 Direct" = get_entrez(
    p3_direct_total
  ),
  "p10 Direct" = get_entrez(
    p10_direct_total
  ),
  "p14 Direct" = get_entrez(
    p14_direct_total
  )
)

# Remove any NULL or empty lists if a stage has no specific genes
gene_list_cluster <- compact(
  gene_list_cluster
)

# =========================================================
# 5. RUN ENRICHMENT
# =========================================================

ck <- compareCluster(
  geneCluster = gene_list_cluster,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

df_go <- as.data.frame(
  ck
)

# =========================================================
# 6. REBUILD GRAPHIC FOR DIRECT SAMPLES
# =========================================================

# Extract top 10 GO terms per stage dynamically for Direct
top_terms_direct <- df_go %>%
  group_by(Cluster) %>%
  slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  pull(Description) %>%
  unique()

plot_data <- df_go %>%
  filter(
    Description %in% top_terms_direct
  ) %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = c(
        "p3 Direct",
        "p10 Direct",
        "p14 Direct"
      )
    ),
    Gene_Count = as.numeric(
      Count
    ),
    padj = as.numeric(
      p.adjust
    )
  )

# Order Y-axis terms by overall minimum adjusted p-value
term_order <- plot_data %>%
  group_by(Description) %>%
  summarize(
    min_padj = min(padj),
    .groups = "drop"
  ) %>%
  arrange(
    desc(min_padj)
  ) %>%
  pull(Description)

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

  # Bright Magenta -> Pink -> Orange -> Yellow Color Scale
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
      0.05
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
    title = "GO enrichment of Direct stage-specific genes",
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

# Print and export Figure 2E
print(
  p_final
)

ggsave(
  file.path(output_dir, "Figure2E_GO_Dotplot.png"),
  p_final,
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  file.path(output_dir, "Figure2E_GO_Dotplot.pdf"),
  p_final,
  width = 10,
  height = 12,
  device = cairo_pdf
)

# =========================================================
# 7. EXPORT GO ENRICHMENT RESULTS TO EXCEL
# =========================================================

df_go_export <- df_go %>%
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

# Filter individual condition datasets
df_p3_direct_go <- df_go_export %>%
  filter(
    Condition == "p3 Direct"
  )

df_p10_direct_go <- df_go_export %>%
  filter(
    Condition == "p10 Direct"
  )

df_p14_direct_go <- df_go_export %>%
  filter(
    Condition == "p14 Direct"
  )

# Terms shown in plot
df_plot_export <- df_go_export %>%
  filter(
    Term_Description %in% top_terms_direct
  )

# Export multi-tab Excel file
write_xlsx(
  list(
    "All_GO_Terms" = df_go_export,
    "Top_Terms_In_Plot" = df_plot_export,
    "P3_Direct_GO" = df_p3_direct_go,
    "P10_Direct_GO" = df_p10_direct_go,
    "P14_Direct_GO" = df_p14_direct_go
  ),
  path = file.path(output_dir, "Figure2E_GO_Results.xlsx")
)

message(
  "Done! Direct plot and Excel file exported successfully."
)

cat(
  "\nFiles saved in:\n",
  getwd(),
  "\n"
)


# Reproducibility metadata
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Figure 2E outputs were saved in: ", normalizePath(output_dir))
