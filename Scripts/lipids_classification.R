# Load required packages
library(stringr)
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input file path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input) || is.null(opt$output)) {
  stop("Both input and output file paths must be provided.")
}

# Read the input file
lipids <- read.delim(opt$input, header = FALSE, stringsAsFactors = FALSE, sep='\t')

# Rename the column to "Lipid" to standardize
names(lipids) <- c("Lipid")

# Define the lipid classification function
classify_lipid <- function(lipid) {
  lipid_clean <- str_trim(lipid)
  
  group <- "Other/Unclassified"
  type <- "Other"
  
  if (str_detect(lipid_clean, regex("^PC", ignore_case = TRUE))) {
    group <- "Glycerophospholipids"
    type <- "Phosphatidylcholine"
  } else if (str_detect(lipid_clean, regex("^PE", ignore_case = TRUE))) {
    group <- "Glycerophospholipids"
    type <- "Phosphatidylethanolamine"
  } else if (str_detect(lipid_clean, regex("^PG", ignore_case = TRUE))) {
    group <- "Glycerophospholipids"
    type <- "Phosphatidylglycerol"
  } else if (str_detect(lipid_clean, regex("^PI", ignore_case = TRUE))) {
    group <- "Glycerophospholipids"
    type <- "Phosphatidylinositol"
  } else if (str_detect(lipid_clean, regex("^PS", ignore_case = TRUE))) {
    group <- "Glycerophospholipids"
    type <- "Phosphatidylserine"
  } else if (str_detect(lipid_clean, regex("^TG", ignore_case = TRUE))) {
    group <- "Glycerolipids"
    type <- "Triacylglycerol"
  } else if (str_detect(lipid_clean, regex("^Cer", ignore_case = TRUE))) {
    group <- "Sphingolipids"
    type <- "Ceramide"
  } else if (str_detect(lipid_clean, regex("^SM", ignore_case = TRUE))) {
    group <- "Sphingolipids"
    type <- "Sphingomyelin"
  } else if (str_detect(lipid_clean, regex("^CE", ignore_case = TRUE))) {
    group <- "Sterol Lipids"
    type <- "Cholesteryl Ester"
  }
  
  return(list(group = group, type = type))
}

# Define the saturation classification function
classify_saturation <- function(lipid) {
  matches <- str_match_all(lipid, "(\\d+):(\\d+)")[[1]]
  
  if (length(matches) == 0) {
    return("Unclassified")
  }
  
  double_bonds <- as.numeric(matches[,2])
  
  if (all(double_bonds == 0)) {
    return("SFA")
  } else if (any(double_bonds == 1)) {
    return("MUFA")
  } else if (any(double_bonds >= 2)) {
    return("PUFA")
  }
  
  return("Unclassified")
}

# Apply the classification functions
results <- lapply(lipids$Lipid, classify_lipid)
lipids$LipidGroup <- sapply(results, function(x) x$group)
lipids$LipidType <- sapply(results, function(x) x$type)
lipids$Saturation <- sapply(lipids$Lipid, classify_saturation)

# Show first few rows to check
head(lipids)

# Write the annotated table to file
write.table(lipids, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
