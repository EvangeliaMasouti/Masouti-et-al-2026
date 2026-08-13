# ==============================================================================
# Cell Marker Expression Analysis and Annotated Heatmap Generation
# Figures 2C and 3C
# ==============================================================================

# Run this script from the project root.
# Expected input:  data/gene.count.all.xlsx
# Outputs:         results/figures/Cell_Markers_Direct.{pdf,png}
#                  results/figures/Cell_Markers_DIV7.{pdf,png}

# --- 1. Project paths ---------------------------------------------------------
input_file <- file.path("data", "gene.count.all.xlsx")
output_dir <- file.path("results", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load libraries --------------------------------------------------------
library(readxl)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(stringr)

# --- 3. Load counts and metadata ---------------------------------------------
df <- read_excel(input_file)

if (!"gene_name" %in% names(df)) {
  stop("The input file must contain a column named 'gene_name'.")
}

# Fail explicitly rather than silently omitting the required Atp1b2 marker.
if (!"Atp1b2" %in% df$gene_name) {
  stop("Atp1b2 was not found in the gene_name column of the input file.")
}

counts_data <- df %>%
  dplyr::select(
    gene_name,
    starts_with("p3direct"),
    starts_with("p10direct"),
    starts_with("p14direct"),
    starts_with("p3div7"),
    starts_with("p10div7"),
    starts_with("p14div7")
  ) %>%
  as.data.frame()

# --- 4. Define categorized markers and explicit order ------------------------
# Atp1b2 is intentionally the first astrocyte marker.
astro_markers <- c(
  "Atp1b2",
  "Gfap",
  "Aqp4",
  "S100b",
  "Aldh1l1",
  "Aldoc",
  "Apoe",
  "Kcnj10",
  "Slc1a2",
  "Slc1a3"
)

stem_markers <- c(
  "Nes",
  "Vim",
  "Sox2",
  "Pax6",
  "Prom1"
)

cycle_markers <- c(
  "Mki67",
  "Pcna",
  "Mcm2",
  "Cdk1",
  "Top2a",
  "Ccnb1"
)

# Master order: astrocytes, stem cells, then cell-cycle markers.
ordered_markers <- c(
  astro_markers,
  stem_markers,
  cycle_markers
)

counts_filtered <- counts_data %>%
  filter(gene_name %in% ordered_markers)

# Deduplicate rows while retaining the highest-expression entry.
counts_unique <- counts_filtered %>%
  mutate(
    total_expression = rowSums(
      dplyr::select(., -gene_name)
    )
  ) %>%
  group_by(gene_name) %>%
  slice_max(
    total_expression,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  dplyr::select(-total_expression) %>%
  as.data.frame()

rownames(counts_unique) <- counts_unique$gene_name
counts_matrix <- as.matrix(counts_unique[, -1, drop = FALSE])

# --- 5. Build row annotations -------------------------------------------------
marker_annotation <- data.frame(
  Gene = rownames(counts_matrix)
) %>%
  mutate(
    Category = case_when(
      Gene %in% astro_markers ~ "Astrocyte",
      Gene %in% stem_markers  ~ "Stem Cell",
      Gene %in% cycle_markers ~ "Cell Cycle",
      TRUE                    ~ "Other"
    )
  ) %>%
  tibble::column_to_rownames("Gene")

annotation_colors <- list(
  Category = c(
    "Astrocyte"  = "#88419d",
    "Stem Cell"  = "#8c96c6",
    "Cell Cycle" = "#b3cde3"
  )
)

sample_info <- data.frame(
  sample = colnames(counts_matrix),
  row.names = colnames(counts_matrix)
) %>%
  mutate(
    Condition = factor(
      if_else(str_detect(sample, "direct"), "Direct", "DIV7"),
      levels = c("Direct", "DIV7")
    ),
    Age = factor(
      case_when(
        str_detect(sample, "p3")  ~ "P3",
        str_detect(sample, "p10") ~ "P10",
        str_detect(sample, "p14") ~ "P14"
      ),
      levels = c("P3", "P10", "P14")
    )
  )

# --- 6. Normalize counts using DESeq2 ----------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info,
  design = ~ Condition + Age
)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Log2 transform absolute normalized expression values using a +1 pseudocount.
mat_log2 <- log2(norm_counts + 1)

# --- 7. Apply strict row ordering and display labels --------------------------
present_genes_ordered <- ordered_markers[
  ordered_markers %in% rownames(mat_log2)
]

mat_log2 <- mat_log2[
  present_genes_ordered,
  ,
  drop = FALSE
]

marker_annotation <- marker_annotation[
  present_genes_ordered,
  ,
  drop = FALSE
]

if (nrow(mat_log2) == 0 || !"Atp1b2" %in% rownames(mat_log2)) {
  stop("Atp1b2 is absent from the heatmap matrix after filtering.")
}

if (rownames(mat_log2)[1] != "Atp1b2") {
  stop("Atp1b2 is not the first marker in the heatmap matrix.")
}

# Custom display labels only. The underlying gene symbols remain unchanged.
custom_gene_labels <- rownames(mat_log2)
custom_gene_labels[custom_gene_labels == "Nes"] <- "Nestin"
custom_gene_labels[custom_gene_labels == "Mki67"] <- "Ki-67"

rownames(mat_log2) <- custom_gene_labels
rownames(marker_annotation) <- custom_gene_labels

# --- 8. Establish palette and split matrices ---------------------------------
pastel_palette <- colorRampPalette(
  c(
    "#f4f6f7",
    "#b3cde3",
    "#8c96c6",
    "#88419d"
  )
)(100)

mat_direct <- mat_log2[
  ,
  grep("direct", colnames(mat_log2), value = TRUE),
  drop = FALSE
]

mat_div7 <- mat_log2[
  ,
  grep("div7", colnames(mat_log2), value = TRUE),
  drop = FALSE
]

# --- 9. Heatmap helper --------------------------------------------------------
make_heatmap <- function(mat, output_stem, title) {
  pdf(
    file.path(output_dir, paste0(output_stem, ".pdf")),
    width = 7,
    height = 8
  )

  pheatmap(
    mat,
    color = pastel_palette,
    scale = "none",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = marker_annotation,
    annotation_colors = annotation_colors,
    border_color = "white",
    main = title,
    angle_col = 45
  )

  dev.off()

  png(
    file.path(output_dir, paste0(output_stem, ".png")),
    width = 2100,
    height = 2400,
    res = 300
  )

  pheatmap(
    mat,
    color = pastel_palette,
    scale = "none",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = marker_annotation,
    annotation_colors = annotation_colors,
    border_color = "white",
    main = title,
    angle_col = 45
  )

  dev.off()
}

# --- 10. Generate Figures 2C and 3C ------------------------------------------
make_heatmap(
  mat_direct,
  "Cell_Markers_Direct",
  "Cell Markers: Direct Isolation (Log2 Counts)"
)

make_heatmap(
  mat_div7,
  "Cell_Markers_DIV7",
  "Cell Markers: DIV 7 Culture (Log2 Counts)"
)

cat(
  "Figures 2C and 3C were generated with Atp1b2 as the first astrocyte marker.\n",
  "Files were saved in: ", normalizePath(output_dir), "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================

