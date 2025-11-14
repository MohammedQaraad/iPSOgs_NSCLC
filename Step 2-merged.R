# Load required libraries
library(tidyverse)

# Read the LUAD and LUSC expression data
luad_data <- read.csv("TCGA-LUAD.star_counts_updated.csv", row.names = 1, check.names = FALSE)
lusc_data <- read.csv("TCGA-LUSC.star_counts_updated.csv", row.names = 1, check.names = FALSE)

# Extract sample IDs (column names) from both datasets
luad_samples <- colnames(luad_data)
lusc_samples <- colnames(lusc_data)


# Find common genes (row names)
common_genes <- intersect(rownames(luad_data), rownames(lusc_data))

# Subset both datasets to include only common genes
luad_data <- luad_data[common_genes, , drop = FALSE]
lusc_data <- lusc_data[common_genes, , drop = FALSE]

# Check if column names (sample IDs) are consistent; if not, proceed with all columns
if (ncol(luad_data) != ncol(lusc_data) || !all(colnames(luad_data) == colnames(lusc_data))) {
  # Combine by binding columns after transposing to align genes, then transpose back
  combined_data <- cbind(luad_data, lusc_data)
} else {
  # If by chance columns match, rbind would work, but this is unlikely given your dims
  combined_data <- rbind(luad_data, lusc_data)
}
dim(luad_data)
dim(lusc_data)
dim(combined_data)
# Save the combined expression data to a new CSV file
write.csv(combined_data, "TCGA-LUAD_LUSC_expression_combined.csv", row.names = TRUE, quote = FALSE)


##########################################################
# Combined Labels
##########################################################

# Load required libraries

# Read the LUAD and LUSC label files
luad_labels <- read.csv("TCGA-LUAD_matched_samples.csv")
lusc_labels <- read.csv("TCGA-LUSC_matched_samples.csv")

# Extract sample IDs and labels
luad_labels <- luad_labels[, c("sample", "label")]
lusc_labels <- lusc_labels[, c("sample", "label")]

# Combine the label data frames
combined_labels <- rbind(luad_labels, lusc_labels)

# Save the combined labels to a new CSV file
write.csv(combined_labels, "TCGA-LUAD_LUSC_labels_combined.csv", row.names = FALSE, quote = FALSE)
