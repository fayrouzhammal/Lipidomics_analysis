# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
library(pheatmap)
library(RColorBrewer) 
library(factoextra)  # For PCA visualization

# Function to normalize data based on user choice
normalize_data <- function(data, method) {
  if (method == "TIC") {
    tic <- colSums(data, na.rm = TRUE)
    mean_tic <- mean(tic, na.rm = TRUE)
    normalized_data <- sweep(data, 2, tic, "/") * mean_tic
  } else if (method == "Median") {
    median_intensity <- apply(data, 2, median, na.rm = TRUE)
    normalized_data <- sweep(data, 2, median_intensity, "/")
  } else if (method == "PQN") {
    reference <- apply(data, 1, median, na.rm = TRUE)
    quotients <- sweep(data, 1, reference, "/")
    normalization_factors <- apply(quotients, 2, median, na.rm = TRUE)
    normalized_data <- sweep(data, 2, normalization_factors, "/")
  } else {
    stop("Unsupported normalization method. Choose from: TIC, Median, PQN.")
  }
  return(normalized_data)
}

# Function to generate RLE plots
generate_rle_plot <- function(data, title, output_file) {
  rle_data <- apply(data, 1, function(x) x - median(x, na.rm = TRUE))
  boxplot(rle_data, main = title, outline = FALSE, col = brewer.pal(9, "Set1"))
  dev.copy(pdf, output_file, width = 10, height = 6)
  dev.off()
}

# Function to generate PCA plots
generate_pca_plot <- function(data, metadata, title, output_file) {
  pca_result <- prcomp(t(data), scale. = TRUE)  # Transpose data for PCA
  pca_data <- as.data.frame(pca_result$x)
  pca_data$Sample_Group <- metadata$`Sample Group`  # Add sample groups for coloring

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample_Group)) +
    geom_point(size = 3) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(legend.position = "top")
  
  ggsave(output_file, plot = p, width = 8, height = 6)
}

# Main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript lipidomics_analysis.R <input_data.csv> <contrasts_file.csv> <output_dir> <normalization_method>")
}

# Input arguments
input_data_file <- args[1]
contrasts_file <- args[2]
output_dir <- args[3]
normalization_method <- args[4]

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
normalization_dir <- file.path(output_dir, "Normalization_Results")
dir.create(normalization_dir, showWarnings = FALSE, recursive = TRUE)
diff_analysis_dir <- file.path(output_dir, "Differential_Analysis_Results")
dir.create(diff_analysis_dir, showWarnings = FALSE, recursive = TRUE)
volcano_dir <- file.path(diff_analysis_dir, "Volcano_Plots")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
heatmap_dir <- file.path(diff_analysis_dir, "Heatmaps")
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Load the raw data
data <- read.csv(input_data_file, check.names = FALSE)
metadata_columns <- c('Sample Name', 'Sample Group')
numeric_data <- data[, !names(data) %in% metadata_columns]

# Generate RLE plot before normalization
generate_rle_plot(numeric_data, "RLE Plot Before Normalization", file.path(normalization_dir, "RLE_Before_Normalization.pdf"))

# Generate PCA plot before normalization
generate_pca_plot(numeric_data, data, "PCA Plot Before Normalization", file.path(normalization_dir, "PCA_Before_Normalization.pdf"))

# Normalize the data
normalized_data <- normalize_data(numeric_data, normalization_method)

# Generate RLE plot after normalization
generate_rle_plot(normalized_data, "RLE Plot After Normalization", file.path(normalization_dir, "RLE_After_Normalization.pdf"))

# Generate PCA plot after normalization
generate_pca_plot(normalized_data, data, "PCA Plot After Normalization", file.path(normalization_dir, "PCA_After_Normalization.pdf"))

# Log2 transformation
log2_transformed_data <- log2(normalized_data + 1)

# Save normalized data
write.csv(cbind(data[, metadata_columns], log2_transformed_data), file.path(normalization_dir, "Normalized_Data.csv"), row.names = FALSE)

# Load contrasts
contrasts <- read.csv(contrasts_file, check.names = FALSE)

# Differential analysis using LIMMA
design <- model.matrix(~ 0 + data$`Sample Group`)
colnames(design) <- gsub("data\\$`Sample Group`", "", colnames(design))
fit <- lmFit(t(log2_transformed_data), design)

# Apply contrasts
contrast_matrix <- makeContrasts(contrasts = contrasts$Contrast, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Function to create volcano plots
create_volcano_plot <- function(results, title, output_file) {
  results$Significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Significant", "Not Significant")
  p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(title = title, x = "Log Fold Change", y = "-log10 Adjusted P-Value") +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")
  ggsave(output_file, plot = p, width = 8, height = 6)
}

# Function to create heatmaps
create_heatmap <- function(data, metadata, title, output_file) {
  pheatmap(
    mat = data,
    annotation_col = metadata,
    main = title,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(50),
    fontsize = 10,
    filename = output_file
  )
}

# Process each contrast
for (i in 1:nrow(contrasts)) {
  contrast_name <- contrasts$Contrast[i]
  results <- topTable(fit2, coef = i, number = Inf, adjust = "fdr")
  
  # Save results
  write.csv(results, file.path(diff_analysis_dir, paste0(contrast_name, "_Results.csv")), row.names = TRUE)
  
  # Volcano plot
  create_volcano_plot(results, contrast_name, file.path(volcano_dir, paste0(contrast_name, "_Volcano.pdf")))
  
  # Heatmap
  significant_genes <- rownames(results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ])
  significant_data <- log2_transformed_data[significant_genes, ]
  create_heatmap(significant_data, data[, metadata_columns, drop = FALSE], contrast_name, file.path(heatmap_dir, paste0(contrast_name, "_Heatmap.pdf")))
}

print("Analysis completed successfully!")
