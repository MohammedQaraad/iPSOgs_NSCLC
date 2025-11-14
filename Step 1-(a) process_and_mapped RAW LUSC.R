# Load required libraries
library(tidyverse)
library(reshape2)
library(ggplot2)

# Read the TSV file
data <- read.table("TCGA-LUSC.star_counts.tsv", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Inspect the first few rows and structure
cat("First few rows of the data:\n")
print(head(data))
cat("\nStructure of the data:\n")
str(data)

# Check if data is likely raw counts or normalized
# Raw counts are typically integers, while normalized data (e.g., TPM, FPKM) are decimals
is_integer <- function(x) all(x == floor(x), na.rm = TRUE)
numeric_cols <- sapply(data, is.numeric)  # Check all columns (samples)
integer_check <- sapply(data[, numeric_cols], is_integer)
cat("\nAre columns integers (suggesting raw counts)?\n")
print(integer_check)

# Calculate summary statistics for numeric columns
summary_stats <- summary(data[, numeric_cols])
cat("\nSummary statistics of expression values:\n")
print(summary_stats)

# Check range of values per sample to infer normalization
range_values <- apply(data[, numeric_cols], 2, range, na.rm = TRUE)
cat("\nRange of values per sample:\n")
print(range_values)

# Melt data for visualization (convert to long format)
data_long <- melt(as.matrix(data), id.vars = NULL, variable.name = "Sample", value.name = "Expression")
data_long$Gene <- rownames(data)[row(data)[,1]]

# Plot overall distribution of expression values
# Histogram of all expression values
ggplot(data_long, aes(x = Expression)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Expression Values (All Genes and Samples)", x = "Expression", y = "Frequency") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

# Save the histogram
ggsave("expression_histogram.png", width = 8, height = 6)

# Boxplot of expression values per sample (first 10 samples for clarity)
subset_samples <- unique(data_long$Var2)[1:min(10, length(unique(data_long$Var2)))]
ggplot(data_long %>% filter(Var2 %in% subset_samples), aes(x = Var2, y = Expression)) +
  geom_boxplot(fill = "lightgreen", outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "Expression Distribution Across Samples", x = "Sample ID", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the boxplot
ggsave("expression_boxplot.png", width = 12, height = 6)

# If data appears to be raw counts, plot log-transformed distribution
if (any(integer_check)) {
  cat("\nData appears to contain raw counts. Plotting log-transformed distribution...\n")
  data_long$LogExpression <- log2(data_long$Expression + 1)  # Add 1 to avoid log(0)
  ggplot(data_long, aes(x = LogExpression)) +
    geom_histogram(bins = 100, fill = "coral", color = "black") +
    theme_minimal() +
    labs(title = "Log2(Expression + 1) Distribution", x = "Log2(Expression + 1)", y = "Frequency") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  ggsave("log_expression_histogram.png", width = 8, height = 6)
}




##############################################
# Replace with real gene ID
##############################################

# Load required libraries
library(tidyverse)

# Read the TSV file with expression data
data <- read.table("TCGA-LUAD.star_counts.tsv", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Read the gene ID mapping file
gene_map <- read.csv("gene_id.csv", header = TRUE, stringsAsFactors = FALSE)

# Create a mapping dictionary using full Ensembl IDs
gene_mapping <- setNames(gene_map$gene, gene_map$id)




# Count duplicate gene symbols
duplicate_counts <- table(gene_map$gene)
duplicates <- duplicate_counts[duplicate_counts > 1]
num_duplicates <- length(duplicates)
cat("Number of unique gene symbols with duplicates:", num_duplicates, "\n")
cat("Details of duplicate gene symbols:\n")
print(duplicates)

# Total number of genes affected by duplicates
total_duplicate_instances <- sum(duplicates - 1)  # Subtract 1 to count extra occurrences
cat("Total number of duplicate gene instances:", total_duplicate_instances, "\n")



######################### fix duplicate

# Create a mapping dictionary using full Ensembl IDs
gene_mapping <- setNames(gene_map$gene, gene_map$id)

# Replace row names with gene symbols
new_row_names <- gene_mapping[rownames(data)]

# Handle duplicates by appending a suffix
row_names_df <- data.frame(ensembl = rownames(data), gene = new_row_names, stringsAsFactors = FALSE)
row_names_df <- row_names_df %>%
  group_by(gene) %>%
  mutate(new_gene = ifelse(is.na(gene), ensembl, 
                           ifelse(duplicated(gene) | duplicated(gene, fromLast = TRUE), 
                                  paste(gene, ave(gene, gene, FUN = seq_along), sep = "_"), gene)))
new_row_names <- row_names_df$new_gene

# Assign new row names to data
rownames(data) <- new_row_names

# Inspect the first few rows to verify the change
cat("First few rows with updated gene names:\n")
print(head(data))

# Save the updated data to a new file
write.table(data, "TCGA-LUSC.star_counts_with_Gene_ID.tsv", sep = "\t", quote = FALSE)

write.csv(data, "TCGA-LUSC.star_counts_updated.csv", row.names = TRUE, quote = FALSE)



##################################################################
# process clinical data
##################################################################

# Load required libraries
library(tidyverse)

# Read the clinical sample IDs
clinical_samples <- read.csv("TCGA-LUSC_clinical.csv", header = FALSE, col.names = "sample")

# Read the expression data (assuming it’s the updated CSV from the previous step)
expression_data <- read.csv("TCGA-LUAD.star_counts_updated.csv", row.names = 1, check.names = FALSE)

# Extract column names (sample IDs) from the expression data
data_samples <- colnames(expression_data)

# Find matching samples, preserving the order from expression_data
matching_samples <- data_samples[data_samples %in% clinical_samples$sample]

# Create a data frame with the matching samples in the original order
final_samples <- data.frame(sample = matching_samples)

# Save the matching samples to a CSV file
write.csv(final_samples, "TCGA-LUSC_matched_samples.csv", row.names = FALSE, quote = FALSE)

