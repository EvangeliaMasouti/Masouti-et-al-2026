# ==============================================================================
# Figure 1
# GFAP, MAP2 and DCX percentages
# Treatment versus positive-control two-tailed unpaired t-tests
# ==============================================================================

rm(list = ls(all.names = TRUE))
graphics.off()

# Run from the project root.
# Expected input: data/EM124 tests.xlsx
project_dir <- getwd()
input_file <- file.path(
  project_dir,
  "data",
  "EM124 tests.xlsx"
)
output_dir <- file.path(
  project_dir,
  "results",
  "Figure_1_ttests_vs_positive_control"
)

# --- 1. Load libraries --------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(writexl)
})

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# --- 2. Read and clean data ---------------------------------------------------
raw_data <- readxl::read_excel(
  input_file,
  .name_repair = "unique"
)

if (ncol(raw_data) < 4) {
  stop("The input file must contain a sample column and three marker columns.")
}

# The first column contains the experiment/sample labels. The remaining columns
# are identified by their marker names, for example GFAP+, MAP2+ and DCX+.
data_clean <- raw_data %>%
  dplyr::rename(sample_id = 1) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id)
  ) %>%
  dplyr::filter(
    !is.na(sample_id),
    sample_id != ""
  ) %>%
  dplyr::mutate(
    # Remove replicate numbers at the end of the sample label.
    # Example: "positive control 1" becomes "positive control".
    treatment = trimws(
      gsub("\\s*\\d+$", "", sample_id)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = -c(sample_id, treatment),
    names_to = "raw_marker",
    values_to = "percentage"
  ) %>%
  dplyr::mutate(
    marker = dplyr::case_when(
      stringr::str_detect(
        raw_marker,
        stringr::regex("GFAP", ignore_case = TRUE)
      ) ~ "GFAP",
      stringr::str_detect(
        raw_marker,
        stringr::regex("MAP2", ignore_case = TRUE)
      ) ~ "MAP2",
      stringr::str_detect(
        raw_marker,
        stringr::regex("DCX", ignore_case = TRUE)
      ) ~ "DCX",
      TRUE ~ NA_character_
    ),
    percentage = suppressWarnings(
      as.numeric(percentage)
    )
  ) %>%
  dplyr::filter(
    !is.na(marker),
    !is.na(percentage)
  )

if (!nrow(data_clean)) {
  stop("No valid GFAP, MAP2 or DCX measurements were detected.")
}

# --- 3. Define the positive-control pattern ----------------------------------
control_pattern <- stringr::regex(
  "positive\\s*control|\\+\\s*control",
  ignore_case = TRUE
)

target_markers <- c(
  "GFAP",
  "MAP2",
  "DCX"
)

# --- 4. Run t-tests for every marker and treatment ---------------------------
all_ttests <- list()

for (marker_name in target_markers) {
  marker_data <- data_clean %>%
    dplyr::filter(marker == marker_name)

  control_values <- marker_data %>%
    dplyr::filter(
      stringr::str_detect(treatment, control_pattern)
    ) %>%
    dplyr::pull(percentage)

  treatment_data <- marker_data %>%
    dplyr::filter(
      !stringr::str_detect(treatment, control_pattern)
    )

  if (length(control_values) < 2) {
    warning(
      "Fewer than two positive-control replicates were found for ",
      marker_name,
      "."
    )
    next
  }

  if (!nrow(treatment_data)) {
    warning(
      "No treatment data were found for ",
      marker_name,
      "."
    )
    next
  }

  result <- treatment_data %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(
      marker = marker_name,
      n_replicates = dplyr::n(),
      n_control = length(control_values),
      control_mean = mean(control_values),
      treatment_mean = mean(percentage),
      mean_difference = treatment_mean - control_mean,
      p_value = if (dplyr::n() >= 2) {
        stats::t.test(
          x = percentage,
          y = control_values,
          alternative = "two.sided",
          var.equal = TRUE
        )$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # BH correction is applied within each marker across its treatment
      # comparisons, matching the original Figure 1 analysis.
      p_adj_fdr = stats::p.adjust(
        p_value,
        method = "BH"
      ),
      is_significant = dplyr::case_when(
        is.na(p_value) ~ "Insufficient data",
        p_adj_fdr < 0.05 ~ "YES (*)",
        TRUE ~ "NO"
      ),
      statistical_test = "Two-tailed unpaired Student's t-test",
      variance_assumption = "Equal variances assumed",
      fdr_scope = "Within marker across treatment comparisons"
    ) %>%
    dplyr::select(
      marker,
      treatment,
      n_replicates,
      n_control,
      control_mean,
      treatment_mean,
      mean_difference,
      p_value,
      p_adj_fdr,
      is_significant,
      statistical_test,
      variance_assumption,
      fdr_scope
    ) %>%
    dplyr::arrange(
      p_adj_fdr,
      p_value
    )

  all_ttests[[marker_name]] <- result
}

if (!length(all_ttests)) {
  stop("No t-tests could be performed.")
}

combined_ttests <- dplyr::bind_rows(all_ttests)

# --- 5. Export results --------------------------------------------------------
excel_output <- c(
  list(
    Combined_All_Markers = combined_ttests,
    Cleaned_Long_Data = data_clean
  ),
  all_ttests
)

output_path <- file.path(
  output_dir,
  "Figure1_ttests_vs_positive_control.xlsx"
)

writexl::write_xlsx(
  excel_output,
  path = output_path
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

cat(
  "Figure 1 t-test analysis completed.\n",
  "Results saved to: ",
  normalizePath(output_path),
  "\n",
  sep = ""
)

# ==============================================================================
# End of script
# ==============================================================================
