# -----------------------------
# 0. Parse Command Line Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R <input_csv_folder> <lipid_classification_file> <output_folder>")
} 
input_folder <- args[1]
classification_file <- args[2]
output_dir <- args[3]

# Create the output directory if it doesn't exist.
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------
# 1. Load required libraries
# -----------------------------
if (!require("fgsea")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("fgsea")
}
library(fgsea)

if (!require("dplyr")) {
  install.packages("dplyr")
}
library(dplyr)

if (!require("tibble")) {
  install.packages("tibble")
  library(tibble)
} else {
  library(tibble)
}

if (!require("ggplot2")) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!require("gridExtra")) {
  install.packages("gridExtra")
}
library(gridExtra)

# -----------------------------
# 2. Function to convert list columns to comma-separated strings
# -----------------------------
convert_list_columns <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if (is.null(x)) {
          return(NA)
        } else {
          paste(x, collapse = ",")
        }
      })
    }
  }
  return(df)
}

# -----------------------------
# 3. Read the Lipid Classification File
# -----------------------------
# This file must have columns: "Lipid", "Saturation", "LipidGroup", and "LipidType"
classification <- read.delim(classification_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Lipid classification file loaded. Sample:\n")
print(head(classification))

# -----------------------------
# 4. List Differential Analysis Files from the Input Folder
# -----------------------------
diff_files <- list.files(path = input_folder, pattern = "\\.csv$", full.names = TRUE)
if (length(diff_files) == 0) {
  stop("No CSV files found in the input folder.")
}
cat("Found", length(diff_files), "CSV files in", input_folder, "\n")

# -----------------------------
# 5. Set fgsea parameters
# -----------------------------
minSetSize <- 5    # minimum set size (adjust as needed)
maxSetSize <- 500  # maximum set size (adjust as needed)
nperm <- 10000     # number of permutations

# -----------------------------
# 6. Process Each Differential Analysis File
# -----------------------------
for (file in diff_files) {
  cat("Processing file:", file, "\n")
  
  # Read the differential analysis file.
  # We assume the first column (with lipid names) is unnamed, so we use row.names = 1.
  diff_results <- read.delim(file, header = TRUE, sep = ",", stringsAsFactors = FALSE,
                             check.names = FALSE, row.names = 1)
  
  # Add lipid names as a proper column for merging.
  diff_results$Lipid <- rownames(diff_results)
  
  # Merge with the classification file.
  merged_data <- merge(diff_results, classification, by = "Lipid", all.x = TRUE)
  
  # Create a ranked vector using the "logFC" column.
  ranked_vector <- merged_data$logFC
  names(ranked_vector) <- merged_data$Lipid
  ranked_vector <- sort(ranked_vector, decreasing = TRUE)
  
  # -----------------------------
  # 7. Define Lipid Sets by Classification
  # -----------------------------
  # Lipid sets by Saturation
  saturation_sets <- merged_data %>%
    group_by(Saturation) %>%
    summarise(Lipids = list(Lipid)) %>%
    deframe()
  
  # Lipid sets by LipidGroup
  group_sets <- merged_data %>%
    group_by(LipidGroup) %>%
    summarise(Lipids = list(Lipid)) %>%
    deframe()
  
  # Lipid sets by LipidType
  type_sets <- merged_data %>%
    group_by(LipidType) %>%
    summarise(Lipids = list(Lipid)) %>%
    deframe()
  
  # -----------------------------
  # 8. Run fgsea on Each Classification
  # -----------------------------
  fgsea_saturation <- fgsea(
    pathways = saturation_sets,
    stats = ranked_vector,
    minSize = minSetSize,
    maxSize = maxSetSize,
    nperm = nperm
  )
  
  fgsea_group <- fgsea(
    pathways = group_sets,
    stats = ranked_vector,
    minSize = minSetSize,
    maxSize = maxSetSize,
    nperm = nperm
  )
  
  fgsea_type <- fgsea(
    pathways = type_sets,
    stats = ranked_vector,
    minSize = minSetSize,
    maxSize = maxSetSize,
    nperm = nperm
  )
  
  # Sort results by adjusted p-value for clarity.
  fgsea_saturation <- fgsea_saturation %>% arrange(padj)
  fgsea_group <- fgsea_group %>% arrange(padj)
  fgsea_type <- fgsea_type %>% arrange(padj)
  
  # Convert list columns (e.g., leadingEdge) to comma-separated strings.
  fgsea_saturation <- convert_list_columns(fgsea_saturation)
  fgsea_group <- convert_list_columns(fgsea_group)
  fgsea_type <- convert_list_columns(fgsea_type)
  
  # -----------------------------
  # 9. Write the Results to Output Files
  # -----------------------------
  base_name <- tools::file_path_sans_ext(basename(file))
  out_sat <- file.path(output_dir, paste0("fgsea_saturation_", base_name, ".tsv"))
  out_group <- file.path(output_dir, paste0("fgsea_group_", base_name, ".tsv"))
  out_type <- file.path(output_dir, paste0("fgsea_type_", base_name, ".tsv"))
  
  write.table(fgsea_saturation, out_sat, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fgsea_group, out_group, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fgsea_type, out_type, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("GSEA outputs written:\n", out_sat, "\n", out_group, "\n", out_type, "\n")
  
  # -----------------------------
  # 10. Plot the Results
  # -----------------------------
  # 10a. Enrichment Curves: For each classification, select the top (most significant) pathway,
  # and plot its enrichment curve.
  
  # Saturation enrichment curve
  if(nrow(fgsea_saturation) > 0 && !is.null(saturation_sets[[fgsea_saturation$pathway[1]]])) {
    top_sat <- fgsea_saturation$pathway[1]
    p_sat <- plotEnrichment(saturation_sets[[top_sat]], ranked_vector) +
             ggtitle(paste("Saturation:", top_sat))
  } else {
    p_sat <- ggplot() + ggtitle("No significant saturation pathway")
  }
  
  # LipidGroup enrichment curve
  if(nrow(fgsea_group) > 0 && !is.null(group_sets[[fgsea_group$pathway[1]]])) {
    top_group <- fgsea_group$pathway[1]
    p_group <- plotEnrichment(group_sets[[top_group]], ranked_vector) +
               ggtitle(paste("LipidGroup:", top_group))
  } else {
    p_group <- ggplot() + ggtitle("No significant lipid group pathway")
  }
  
  # LipidType enrichment curve
  if(nrow(fgsea_type) > 0 && !is.null(type_sets[[fgsea_type$pathway[1]]])) {
    top_type <- fgsea_type$pathway[1]
    p_type <- plotEnrichment(type_sets[[top_type]], ranked_vector) +
              ggtitle(paste("LipidType:", top_type))
  } else {
    p_type <- ggplot() + ggtitle("No significant lipid type pathway")
  }
  
  # Arrange the three enrichment curve plots side by side
  p_enrichment <- grid.arrange(p_sat, p_group, p_type, ncol = 3)
  out_enrich_pdf <- file.path(output_dir, paste0("Enrichment_Curves_", base_name, ".pdf"))
  ggsave(out_enrich_pdf, p_enrichment, width = 15, height = 5)
  
  # 10b. NES Bar Plots: Create horizontal bar plots for significant pathways (padj < 0.05)
  # for each classification.
  padj_cutoff <- 0.05
  
  fgsea_sat_sig <- fgsea_saturation %>% filter(padj < padj_cutoff) %>% arrange(desc(NES))
  p_bar_sat <- ggplot(fgsea_sat_sig, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Significant Saturation Pathways (NES)",
         x = "Saturation", y = "Normalized Enrichment Score") +
    theme_minimal()
  
  fgsea_group_sig <- fgsea_group %>% filter(padj < padj_cutoff) %>% arrange(desc(NES))
  p_bar_group <- ggplot(fgsea_group_sig, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Significant LipidGroup Pathways (NES)",
         x = "Lipid Group", y = "Normalized Enrichment Score") +
    theme_minimal()
  
  fgsea_type_sig <- fgsea_type %>% filter(padj < padj_cutoff) %>% arrange(desc(NES))
  p_bar_type <- ggplot(fgsea_type_sig, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Significant LipidType Pathways (NES)",
         x = "Lipid Type", y = "Normalized Enrichment Score") +
    theme_minimal()
  
  # Arrange the three bar plots vertically
  p_barplots <- grid.arrange(p_bar_sat, p_bar_group, p_bar_type, ncol = 1)
  out_bar_pdf <- file.path(output_dir, paste0("NES_Barplots_", base_name, ".pdf"))
  ggsave(out_bar_pdf, p_barplots, width = 8, height = 12)
  
  cat("Plots saved:\n", out_enrich_pdf, "\n", out_bar_pdf, "\n\n")
}
