
#########################################
# DEG DESeq2
#########################################
# 1. Install and load necessary packages
# If you haven't installed them, uncomment the following lines:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# install.packages("tibble")

library(DESeq2)
library(dplyr)
library(tibble) # For easier data frame manipulation, especially tibble::rownames_to_column

# --- User-provided file paths ---
expression_file <- "TCGA-LUAD_LUSC_expression_combined.csv"
survival_file <- "TCGA-LUAD_LUSC_labels_combined.csv"

# 2. Load the expression data
cat("Loading expression data from:", expression_file, "\n")
combined_expression <- read.csv(expression_file, row.names = 1)
cat("Dimensions of combined expression dataset:", dim(combined_expression), "\n")

# Ensure expression data contains only numeric values (counts)
# It's crucial for DESeq2 that counts are non-negative integers.
# If your data are not raw counts (e.g., TPM, FPKM, or log-transformed), DESeq2 is not appropriate.

# --- NEW CHECK FOR NEGATIVE VALUES ---
if (any(combined_expression < 0)) {
  stop("Error: Negative values detected in your expression data. DESeq2 requires raw count data, which must be non-negative integers. Please provide untransformed raw counts or consider using a different differential expression tool for normalized data.")
}

# Handle any NA values by setting to 0 (should be done AFTER negative value check if you expect raw counts)
combined_expression[is.na(combined_expression)] <- 0

# Round to integers (if not already integers, and if they are indeed counts)
# This step is only appropriate if your data are raw counts that might have been read as floats.
# If your data are already normalized and contain decimals, rounding them will not make them suitable for DESeq2.
combined_expression <- round(combined_expression)

combined_expression <- as.matrix(combined_expression)


# 3. Load the survival/label data
cat("Loading survival data from:", survival_file, "\n")
combined_survival <- read.csv(survival_file)
cat("Dimensions of combined survival dataset:", dim(combined_survival), "\n")

# 4. Prepare metadata (coldata) for DESeq2
# The 'dataset' column in combined_survival.csv is assumed to contain the group labels ("LUAD", "LUSC").
# The 'sample_id' column in combined_survival.csv should match the column names of combined_expression.
# Adjust 'sample_id_column_name' if your sample ID column has a different name.
sample_id_column_name <- "sample" # <--- IMPORTANT: Adjust this if your sample ID column is named differently

# Check if the sample_id_column_name exists in combined_survival
if (!sample_id_column_name %in% colnames(combined_survival)) {
  stop(paste0("Error: Column '", sample_id_column_name, "' not found in combined_survival.csv. Please check your file or adjust 'sample_id_column_name' variable."))
}

# Ensure the order of samples in coldata matches the order of columns in the expression matrix
# This is CRUCIAL for DESeq2.
# First, get sample names from expression data
colnames(combined_expression) <- gsub("\\.", "-", colnames(combined_expression))
expression_sample_names <- colnames(combined_expression)

# Filter and order combined_survival to match expression_sample_names
coldata <- combined_survival %>%
  filter(!!sym(sample_id_column_name) %in% expression_sample_names) %>%
  column_to_rownames(var = sample_id_column_name) %>%
  as.data.frame() # Ensure it's a data.frame for DESeq2

# Reorder coldata rows to match expression matrix columns
coldata <- coldata[expression_sample_names, , drop = FALSE]

# Check if 'label' column exists in coldata
if (!"label" %in% colnames(coldata)) {
  stop("Error: Column 'dataset' not found in the processed coldata. Please ensure combined_survival.csv has a 'dataset' column.")
}

# Convert 'dataset' to a factor
coldata$dataset <- factor(coldata$label)

# Verify that column names of expression matrix match row names of coldata
if (!all(colnames(combined_expression) == rownames(coldata))) {
  stop("Sample names in expression matrix and coldata do not match or are not in the same order. Please check your data preparation.")
} else {
  cat("Sample names in expression matrix and coldata match and are in order.\n")
}

# 5. Create a DESeqDataSet object
cat("Creating DESeqDataSet object...\n")
dds <- DESeqDataSetFromMatrix(countData = combined_expression,
                              colData = coldata,
                              design = ~ dataset)

# Pre-filtering: Remove genes with very low counts across all samples
# This helps to reduce the number of tests and improve statistical power.
# A common threshold is to keep genes that have at least 10 counts in total across all samples.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Removed", nrow(combined_expression) - nrow(dds), "genes with low counts.\n")
cat("Remaining genes for analysis:", nrow(dds), "\n")

# 6. Run the DESeq analysis
cat("Running DESeq analysis (this may take a while for large datasets)...\n")
dds <- DESeq(dds)
cat("DESeq analysis complete.\n")

# 7. Extract results
# You can specify the contrast for comparison (e.g., LUAD vs LUSC)
# Replace 'LUAD' and 'LUSC' with the exact names of your levels in the 'dataset' column if they are different.
cat("Extracting results for 'dataset' comparison...\n")
res <- results(dds, contrast = c("dataset", "LUAD", "LUSC")) # Example: Comparing LUAD vs LUSC

# Order results by adjusted p-value
res_ordered <- res[order(res$padj),]

# Summarize results
cat("\nSummary of DESeq2 results:\n")
summary(res)

# 8. Display top differentially expressed genes
cat("\nTop 20 differentially expressed genes (by adjusted p-value):\n")
print(head(res_ordered, 20))

# You can save the results to a CSV file
write.csv(as.data.frame(res_ordered), file = "DEG_results_LUAD_LUSC.csv")
cat("\nFull DEG results saved to DEG_results_LUAD_LUSC.csv\n")

# Further analysis: Volcano plot, Heatmap, etc.
# For example, to get genes with padj < 0.05 and |log2FoldChange| > 1
significant_genes <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nNumber of significant genes (padj < 0.05, |log2FoldChange| > 1):", nrow(significant_genes), "\n")
cat("\nFirst few significant genes:\n")
print(head(significant_genes))


###########################
# Visulation
###########################

# BiocManager::install("EnhancedVolcano")
library(DESeq2)
library(dplyr)
library(tibble) # For easier data frame manipulation, especially tibble::rownames_to_column
library(EnhancedVolcano) # For volcano plot
library(ggplot2) # For general plotting
library(pheatmap) # For heatmap visualization


# Slice combined_expression for only significant genes
if (nrow(significant_genes) > 0) {
  significant_expression <- combined_expression[rownames(significant_genes), ]
  write.csv(significant_expression, file = "significant_genes_expression_LUAD_LUSC.csv")
  cat("\nExpression data for significant genes saved to significant_genes_expression_LUAD_LUSC.csv\n")
} else {
  cat("\nNo significant genes found based on the criteria (padj < 0.05, |log2FoldChange| > 1).\n")
}

# --- NEW SECTION: Volcano Plot for Journal Publication ---
cat("\nGenerating Volcano Plot...\n")

# Define custom colors for the volcano plot
# You can adjust these colors as needed for your publication
my_colors <- c(
  "red" = "red",       # Upregulated
  "blue" = "blue",     # Downregulated
  "grey" = "grey"      # Not significant
)

volcano_plot <- EnhancedVolcano(
  res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano Plot: LUAD vs LUSC Differential Gene Expression',
  subtitle = 'Adjusted p-value < 0.05, |Log2 Fold Change| > 1',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 1.5,
  labSize = 3.0,
  colAlpha = 0.8,
  legendLabels = c('Not significant', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
  col = c('grey', 'forestgreen', 'royalblue', 'red2'), # Default colors: not significant, FC only, p-value only, both
  boxedLabels = FALSE, # Avoid overlapping labels
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  # Removed 'border', 'borderWidth', 'borderColor' as they are not supported in all versions
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  # Customize axis titles and font sizes for publication
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~italic(adj.P.Val)),
  caption = paste0('Total genes = ', nrow(res), '\n',
                   'Up-regulated = ', sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE), '\n',
                   'Down-regulated = ', sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)),
  # Adjust text sizes for journal quality
  axisLabSize = 12,
  titleLabSize = 14,
  subtitleLabSize = 10,
  captionLabSize = 8,
  legendLabSize = 10,
  legendIconSize = 4.0
)

# Save the volcano plot with high quality
ggsave(
  filename = "volcano_plot_LUAD_vs_LUSC.tiff",
  plot = volcano_plot,
  device = "tiff",
  width = 7, # inches
  height = 8, # inches
  dpi = 300 # dots per inch (high resolution)
)
cat("Volcano plot saved as volcano_plot_LUAD_vs_LUSC.tiff\n")


# --- NEW SECTION: Heatmap of Top 20 Differentially Expressed Genes ---
cat("\nGenerating Heatmap for Top 20 Differentially Expressed Genes...\n")

# Get the top 20 genes based on adjusted p-value
top20_genes <- head(res_ordered, 20)
top20_gene_ids <- rownames(top20_genes)

# Extract expression data for these top 20 genes
# Ensure you use the *normalized* counts for visualization, not raw counts.
# DESeq2 provides variance-stabilized transformation (vst) or rlog transformation for this.
# Let's use vst for visualization.
vsd <- vst(dds, blind = FALSE) # Variance stabilizing transformation
expression_for_heatmap <- assay(vsd)[top20_gene_ids, ]

# Scale the expression data for better visualization in heatmap (row-wise scaling)
scaled_expression_for_heatmap <- t(scale(t(expression_for_heatmap)))

# Prepare annotation for columns (samples)
# We'll use the 'dataset' column from coldata
annotation_col <- data.frame(Dataset = coldata$dataset)
rownames(annotation_col) <- rownames(coldata)

# Define colors for the annotation (optional, but good for publication)
ann_colors = list(
  Dataset = c(LUAD = "orange", LUSC = "purple") # Adjust colors as needed
)

# Generate the heatmap
pheatmap(
  scaled_expression_for_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE, # Usually too many samples to show names
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Top 20 Differentially Expressed Genes (LUAD vs LUSC)",
  fontsize_row = 8, # Adjust font size for gene names
  fontsize_col = 6, # Adjust font size for sample names if shown
  fontsize = 10, # General font size for title/legend
  color = colorRampPalette(c("blue", "white", "red"))(100), # Blue-white-red color scale
  filename = "heatmap_top20_genes_LUAD_vs_LUSC.tiff", # Save heatmap directly
  width = 8, # inches
  height = 10, # inches
  dpi = 300 # high resolution
)
cat("Heatmap for top 20 genes saved as heatmap_top20_genes_LUAD_vs_LUSC.tiff\n")


# install.packages("patchwork") # For combining ggplot2 plots
# install.packages("pheatmap")
library(DESeq2)
library(dplyr)
library(tibble) # For easier data frame manipulation, especially tibble::rownames_to_column
library(EnhancedVolcano) # For volcano plot
library(ggplot2) # For general plotting
library(pheatmap) # For heatmap visualization
library(reshape2) # For melting data for ggplot2
library(patchwork) # For combining ggplot2 plots

# --- NEW SECTION: Heatmap of Top 20 Differentially Expressed Genes ---
cat("\nGenerating Heatmap for Top 20 Differentially Expressed Genes...\n")

# Get the top 20 genes based on adjusted p-value
top20_genes <- head(res_ordered, 20)
top20_gene_ids <- rownames(top20_genes)

# Extract expression data for these top 20 genes
# Ensure you use the *normalized* counts for visualization, not raw counts.
# DESeq2 provides variance-stabilized transformation (vst) or rlog transformation for this.
# Let's use vst for visualization.
vsd <- vst(dds, blind = FALSE) # Variance stabilizing transformation
expression_for_heatmap <- assay(vsd)[top20_gene_ids, ]

# Scale the expression data for better visualization in heatmap (row-wise scaling)
scaled_expression_for_heatmap <- t(scale(t(expression_for_heatmap)))

# Prepare annotation for columns (samples)
# We'll use the 'dataset' column from coldata
annotation_col <- data.frame(Dataset = coldata$dataset)
rownames(annotation_col) <- rownames(coldata)

# Define colors for the annotation (optional, but good for publication)
ann_colors = list(
  Dataset = c(LUAD = "orange", LUSC = "purple") # Adjust colors as needed
)

# Generate the heatmap
pheatmap(
  scaled_expression_for_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE, # Usually too many samples to show names
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Top 20 Differentially Expressed Genes (LUAD vs LUSC)",
  fontsize_row = 8, # Adjust font size for gene names
  fontsize_col = 6, # Adjust font size for sample names if shown
  fontsize = 10, # General font size for title/legend
  color = colorRampPalette(c("blue", "white", "red"))(100), # Blue-white-red color scale
  filename = "heatmap_top20_genes_LUAD_vs_LUSC.tiff", # Save heatmap directly
  width = 8, # inches
  height = 10, # inches
  dpi = 300 # high resolution
)
cat("Heatmap for top 20 genes saved as heatmap_top20_genes_LUAD_vs_LUSC.tiff\n")

# --- NEW SECTION: Boxplots for Top 10 Differentially Expressed Genes ---
cat("\nGenerating Boxplots for Top 10 Differentially Expressed Genes...\n")

# Get the top 10 genes based on adjusted p-value
top10_genes <- head(res_ordered, 10)
top10_gene_ids <- rownames(top10_genes)

# Extract variance-stabilized transformed (VST) expression data for these top 10 genes
# This is suitable for visualization as it accounts for library size and variance.
expression_for_boxplot <- assay(vsd)[top10_gene_ids, ]

# Convert to a data frame and add gene names as a column
expression_df <- as.data.frame(expression_for_boxplot) %>%
  rownames_to_column(var = "Gene")

# Melt the data frame to long format for ggplot2
# This creates columns for 'Gene', 'Sample', and 'Expression_Level'
melted_expression <- melt(expression_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression_Level")

# Merge with coldata to get the 'dataset' information for each sample
merged_data_for_boxplot <- left_join(melted_expression,
                                     coldata %>% rownames_to_column(var = "Sample"),
                                     by = "Sample")

# Create a list to store individual plots
plot_list <- list()

# Loop through each of the top 10 genes to create a boxplot
for (gene in top10_gene_ids) {
  p <- ggplot(subset(merged_data_for_boxplot, Gene == gene),
              aes(x = dataset, y = Expression_Level, fill = dataset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Hide outliers from boxplot, show as jittered points
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) + # Add jittered points
    scale_fill_manual(values = c("LUAD" = "#FF9933", "LUSC" = "#009966")) + # Custom colors for LUAD and LUSC
    labs(title = gene, x = "", y = "") + # Gene name as title, no x/y labels for individual plots
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none", # No legend for individual plots
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.border = element_rect(colour = "black", fill=NA, size=0.5) # Add a border
    )
  plot_list[[gene]] <- p
}

# Combine all plots into a single figure using patchwork
# Adjust layout (e.g., 2 rows, 5 columns) based on the number of plots (10 genes)
combined_boxplot_plot <- wrap_plots(plot_list, ncol = 5, nrow = 2) +
  plot_annotation(
    title = 'Top 10 Differentially Expressed Genes (LUAD vs. LUSC)',
    caption = 'Log2 Expression Level',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.0) # Align caption to left
    )
  ) & theme(axis.title.y = element_text(size = 10, angle = 90)) # Apply common y-axis label

# Save the combined boxplot plot with high quality
ggsave(
  filename = "boxplot_top10_genes_LUAD_vs_LUSC.tiff",
  plot = combined_boxplot_plot,
  device = "tiff",
  width = 12, # inches (adjust as needed for 10 plots)
  height = 8, # inches (adjust as needed)
  dpi = 300 # dots per inch (high resolution)
)
cat("Boxplots for top 10 genes saved as boxplot_top10_genes_LUAD_vs_LUSC.tiff\n")






########################################################################################
# Slice only Significant Genes from GEO (Test data)
########################################################################################
library(data.table)
library(biomaRt)
library(dplyr)

# --- Step 1: Load expression matrix ---
fpkm_data <- read.csv("GSE81089_log2FPKM_with_GeneSymbols.csv", row.names = 1)
DESeq_significant_genes <- read.csv("significant_genes_expression_LUAD_LUSC.csv", row.names = 1)

# Extract gene IDs (row names)
DESeq_significant_genes <- rownames(DESeq_significant_genes)
fpkm_data_genes <- rownames(fpkm_data)

# Find matched (common) genes
matched_genes <- intersect(DESeq_significant_genes, fpkm_data_genes)

# Find genes unique to LUAD (not in LUSC)
unique_luad <- setdiff(DESeq_significant_genes, fpkm_data_genes)

# Find genes unique to LUSC (not in LUAD)
unique_lusc <- setdiff(fpkm_data_genes, DESeq_significant_genes)

# Print results
cat("Number of matched genes:", length(matched_genes), "\n")
cat("Matched genes:\n")
print(matched_genes)
cat("\nNumber of genes unique to LUAD:", length(unique_luad), "\n")
cat("Genes unique to LUAD:\n")
print(unique_luad)
cat("\nNumber of genes unique to LUSC:", length(unique_lusc), "\n")
cat("Genes unique to LUSC:\n")
print(unique_lusc)
library(dplyr)

# Optionally, save results to files
write.csv(matched_genes, "matched_genes.csv", row.names = FALSE)
write.csv(unique_luad, "DESeq_significant_genes_Not_in_GEO.csv", row.names = FALSE)
write.csv(unique_lusc, "GEO_genes_Not_in_DESeq.csv", row.names = FALSE)

################################################
# Slice Only Matched Genes
###############################################

# Slice the datasets to keep only matched genes
# luad_matched <- luad_data[matched_genes, , drop = FALSE]
lusc_matched <- fpkm_data[matched_genes, , drop = FALSE]
lusc_matched <- round(lusc_matched)
# Print number of matched genes and dimensions of sliced datasets
cat("Number of matched genes:", length(matched_genes), "\n")
# cat("Dimensions of sliced LUAD dataset:", dim(luad_matched), "\n")
cat("Dimensions of sliced LUSC dataset:", dim(lusc_matched), "\n")

# Optionally, save the sliced datasets to new CSV files
# write.csv(luad_matched, "LUAD_matched_genes.csv")
write.csv(lusc_matched, "GEO_matched_Genes_expression_data_test.csv")


###############################################################
# Slice MAtched Genes from DESeq_significant_genes to equivelnt 
#################################################################
DEseq_matched <- DESeq_significant_genes[matched_genes, , drop = FALSE]
DEseq_matched <- round(DEseq_matched)
# Print number of matched genes and dimensions of sliced datasets
cat("Number of matched genes:", length(matched_genes), "\n")
# cat("Dimensions of sliced LUAD dataset:", dim(luad_matched), "\n")
cat("Dimensions of sliced LUSC dataset:", dim(DEseq_matched), "\n")

# Optionally, save the sliced datasets to new CSV files
# write.csv(luad_matched, "LUAD_matched_genes.csv")
write.csv(DEseq_matched, "DEseq_matched_Genes_expression_data_trian_validation.csv")
