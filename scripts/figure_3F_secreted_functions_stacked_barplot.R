# ==============================================================================
# Figure 3F
# Secreted-gene function analysis and stacked barplot
#
# Input:
#   unique genes div7.xlsx
#   secreted_genes_summary_with_functions.xlsx
#
# Secreted-gene classification:
#   A gene is retained if it is present in the manual secreted-gene list OR
#   annotated to one of the extracellular cellular-component GO terms below.
#
# Displayed stages:
#   p3div7, p10div7 and p14div7
# ==============================================================================

rm(list = ls(all.names = TRUE))
graphics.off()

# Run from the project root.
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(
  project_dir,
  "results",
  "Figure_3F_Secreted_Functions"
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
  library(tidyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
  library(ggplot2)
  library(writexl)
})

# --- 2. Input files -----------------------------------------------------------
unique_genes_file <- file.path(
  data_dir,
  "unique genes div7.xlsx"
)

manual_secreted_file <- file.path(
  data_dir,
  "secreted_genes_summary_with_functions.xlsx"
)

if (!file.exists(unique_genes_file)) {
  stop("Input file not found: ", unique_genes_file)
}

if (!file.exists(manual_secreted_file)) {
  stop("Manual secreted-gene file not found: ", manual_secreted_file)
}

# --- 3. Load and reshape unique stage-gene lists -----------------------------
message("Loading unique stage-gene file...")

unique_div7_raw <- readxl::read_excel(
  unique_genes_file,
  .name_repair = "unique"
)

unique_div7 <- unique_div7_raw %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Group",
    values_to = "Gene"
  ) %>%
  dplyr::mutate(
    Group = as.character(Group),
    Gene = trimws(as.character(Gene))
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::distinct(Group, Gene)

genes_vec <- sort(unique(unique_div7$Gene))

if (!length(genes_vec)) {
  stop("No genes were found in the unique stage-gene file.")
}

# --- 4. Retrieve GO cellular-component annotations ---------------------------
go_data <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = genes_vec,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GO", "ONTOLOGY")
) %>%
  dplyr::filter(
    !is.na(SYMBOL),
    !is.na(GO),
    ONTOLOGY == "CC"
  ) %>%
  dplyr::distinct(SYMBOL, GO, ONTOLOGY)

go_details <- AnnotationDbi::select(
  GO.db,
  keys = sort(unique(go_data$GO)),
  columns = "TERM",
  keytype = "GOID"
) %>%
  dplyr::filter(
    !is.na(GOID),
    !is.na(TERM)
  ) %>%
  dplyr::distinct(GOID, .keep_all = TRUE)

gene_functions_go <- go_data %>%
  dplyr::left_join(
    go_details,
    by = c("GO" = "GOID")
  ) %>%
  dplyr::filter(!is.na(TERM)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    # Preserve one fallback function per gene, as in the original analysis.
    GO_Function = first(TERM),
    .groups = "drop"
  )

# Extracellular cellular-component terms used to identify predicted secreted
# genes. These are annotation-based candidates, not experimentally confirmed
# secreted proteins.
secreted_GO_ids <- c(
  "GO:0005576", # extracellular region
  "GO:0005615", # extracellular space
  "GO:0044421"  # extracellular region part
)

predicted_secreted_genes <- go_data %>%
  dplyr::filter(GO %in% secreted_GO_ids) %>%
  dplyr::pull(SYMBOL) %>%
  unique()

# --- 5. Load and validate manual secretome annotations -----------------------
message("Loading manual secreted-gene reference file...")

manual <- readxl::read_excel(
  manual_secreted_file,
  .name_repair = "unique"
) %>%
  as.data.frame()

colnames(manual)[tolower(colnames(manual)) == "gene"] <- "Gene"
colnames(manual)[tolower(colnames(manual)) == "function"] <- "Function"

if (!all(c("Gene", "Function") %in% colnames(manual))) {
  stop(
    "The manual secreted-gene file must contain columns named Gene and Function."
  )
}

manual <- manual %>%
  dplyr::transmute(
    Gene = trimws(as.character(Gene)),
    Function = trimws(as.character(Function))
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::mutate(
    Function = dplyr::na_if(Function, "")
  ) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

# --- 6. Build the final secreted-gene table ----------------------------------
final_table <- unique_div7 %>%
  dplyr::mutate(
    is_secreted =
      Gene %in% predicted_secreted_genes |
      Gene %in% manual$Gene
  ) %>%
  dplyr::filter(is_secreted) %>%
  dplyr::left_join(
    manual,
    by = "Gene"
  ) %>%
  dplyr::left_join(
    gene_functions_go,
    by = c("Gene" = "SYMBOL")
  ) %>%
  dplyr::mutate(
    Function = dplyr::coalesce(
      Function,
      GO_Function,
      "Unclassified secreted gene"
    ),
    Status = "TRUE"
  ) %>%
  dplyr::select(
    Gene,
    Function,
    Group,
    Status,
    is_secreted
  ) %>%
  dplyr::distinct()

conditions_to_plot <- c(
  "p3div7",
  "p10div7",
  "p14div7"
)

missing_conditions <- setdiff(
  conditions_to_plot,
  unique(final_table$Group)
)

if (length(missing_conditions) > 0) {
  stop(
    "The following plotting groups were not found: ",
    paste(missing_conditions, collapse = ", ")
  )
}

# --- 7. Aggregate functions for the stacked barplot --------------------------
plot_data_long <- final_table %>%
  dplyr::filter(
    Group %in% conditions_to_plot,
    Status == "TRUE"
  ) %>%
  dplyr::count(
    Group,
    Function,
    name = "Gene_Count"
  ) %>%
  dplyr::mutate(
    Condition = factor(
      Group,
      levels = conditions_to_plot
    )
  ) %>%
  dplyr::select(
    Condition,
    Function,
    Gene_Count
  )

function_order <- plot_data_long %>%
  dplyr::group_by(Function) %>%
  dplyr::summarise(
    Total = sum(Gene_Count),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Total),
    Function
  ) %>%
  dplyr::pull(Function)

plot_data_long <- plot_data_long %>%
  dplyr::mutate(
    Function = factor(
      Function,
      levels = function_order
    )
  )

# --- 8. Define the color palette ---------------------------------------------
distinct_colors <- c(
  "#b3cde3",
  "#8c96c6",
  "#88419d",
  "#cbd5e8",
  "#e0f3db",
  "#a8ddb5",
  "#4eb3d3",
  "#fdcdac",
  "#cbd5e1",
  "#fbcfe8",
  "#c084fc",
  "#fed976",
  "#feb24c",
  "#94a3b8",
  "#64748b"
)

color_palette <- colorRampPalette(distinct_colors)(
  length(function_order)
)

names(color_palette) <- function_order

# --- 9. Generate Figure 3F ----------------------------------------------------
stacked_barplot <- ggplot(
  plot_data_long,
  aes(
    x = Condition,
    y = Gene_Count,
    fill = Function
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = NA,
    width = 0.55
  ) +
  scale_fill_manual(
    values = color_palette,
    name = "Biological Functions",
    drop = FALSE
  ) +
  labs(
    title = "Distribution of Secreted Functions Across Stages",
    x = "Developmental Stage (DIV7)",
    y = "Number of Secreted Genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15,
      margin = margin(b = 15)
    ),
    axis.text = element_text(
      color = "black",
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

print(stacked_barplot)

# --- 10. Save outputs ---------------------------------------------------------
ggsave(
  filename = file.path(
    output_dir,
    "Figure3F_Secreted_Functions_Stacked_Barplot.png"
  ),
  plot = stacked_barplot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure3F_Secreted_Functions_Stacked_Barplot.pdf"
  ),
  plot = stacked_barplot,
  width = 11,
  height = 7,
  bg = "white"
)

writexl::write_xlsx(
  list(
    Secreted_Gene_Table = final_table,
    Plot_Data = plot_data_long,
    GO_Function_Annotations = gene_functions_go,
    Predicted_Secreted_Genes = data.frame(
      Gene = predicted_secreted_genes,
      stringsAsFactors = FALSE
    ),
    Manual_Secreted_Genes = manual
  ),
  path = file.path(
    output_dir,
    "Figure3F_Secreted_Functions_Audit.xlsx"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

cat(
  "Figure 3F secreted-function stacked barplot was saved in: ",
  normalizePath(output_dir),
  "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================
