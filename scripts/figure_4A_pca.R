# ==============================================================================
# Figure 4A
# PCA of directly isolated and DIV7 cultured astrocytes
# ==============================================================================

rm(list = ls(all.names = TRUE))
graphics.off()

# Run from the project root.
# Expected input: data/pca.data.xlsx
project_dir <- getwd()
input_file <- file.path(
  project_dir,
  "data",
  "pca.data.xlsx"
)
output_dir <- file.path(
  project_dir,
  "results",
  "Figure_4A_PCA"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# --- 1. Load libraries --------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# --- 2. Read and validate PCA data -------------------------------------------
pca_data <- readxl::read_excel(
  input_file,
  .name_repair = "unique"
)

required_columns <- c(
  "sample",
  "group",
  "pc1",
  "pc2"
)

missing_columns <- setdiff(required_columns, colnames(pca_data))

if (length(missing_columns) > 0) {
  stop(
    "The following required PCA columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

pca_data <- pca_data %>%
  dplyr::mutate(
    sample = as.character(sample),
    group = as.character(group),
    pc1 = as.numeric(pc1),
    pc2 = as.numeric(pc2)
  )

if (anyNA(pca_data[, required_columns])) {
  stop("Missing or non-numeric values were detected in the PCA data.")
}

# --- 3. Lock group order and define the Figure 4A palette --------------------
group_levels <- c(
  "p3div7",
  "p10div7",
  "p14div7",
  "p3direct",
  "p10direct",
  "p14direct"
)

unknown_groups <- setdiff(unique(pca_data$group), group_levels)

if (length(unknown_groups) > 0) {
  stop(
    "Unexpected group labels detected: ",
    paste(unknown_groups, collapse = ", ")
  )
}

pca_data$group <- factor(
  pca_data$group,
  levels = group_levels
)

custom_colors <- c(
  "p3div7" = "#008B8B",
  "p10div7" = "#4B0082",
  "p14div7" = "#C71585",
  "p3direct" = "powderblue",
  "p10direct" = "orchid",
  "p14direct" = "pink"
)

# --- 4. Generate Figure 4A ----------------------------------------------------
p <- ggplot(
  pca_data,
  aes(
    x = pc1,
    y = pc2,
    color = group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 0.5
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    aes(label = sample),
    size = 2.8,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = custom_colors,
    drop = FALSE,
    name = "Group"
  ) +
  labs(
    x = "PC1 (60.10%)",
    y = "PC2 (13.70%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f5f5f5"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = grid::unit(0.8, "lines"),
    plot.margin = ggplot2::margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5,
      unit = "mm"
    )
  )

print(p)

# --- 5. Save Figure 4A -------------------------------------------------------
ggsave(
  filename = file.path(
    output_dir,
    "Figure4A_PCA.png"
  ),
  plot = p,
  width = 5,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure4A_PCA.pdf"
  ),
  plot = p,
  width = 5,
  height = 5
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

cat(
  "Figure 4A PCA plots were saved in: ",
  normalizePath(output_dir),
  "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================
