# Clear workspace and load required library
rm(list = ls())
library(data.table)

# Define directories
data_dir <- "/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/"
output_dir <- "/Users/tianyuzhang/Documents/HMC/tutorial/"
input_data_path <- paste0(data_dir, "data/")
intermediate_path <- paste0(input_data_path, "intermediate_data/")

# Load gene clustering info
clustering_path <- paste0(input_data_path, "module_list_df_2000_genes.csv")
clustering <- fread(clustering_path)
clustering[, gene_name := gsub("\\.", "-", gene_name)]

# Identify and remove mitochondrial genes (cluster 42)
mt_genes <- clustering[cluster_index == 42, gene_name]
clustering <- clustering[cluster_index != 42]

# Load residual expression matrix
resid_path <- paste0(intermediate_path, "residual_matrix_small.rds")
residual_subset <- readRDS(resid_path)
residual_subset <- data.table(residual_subset)[, -c("ID", "Cell_cycle_phase")]

# Standardize column names
setnames(residual_subset, gsub("\\.", "-", names(residual_subset)))

# Remove mitochondrial genes
residual_subset <- residual_subset[, !..mt_genes]
cat("Residual matrix after removing mt genes:", dim(residual_subset), "\n")

# Subsample 500 non-targeting cells
control_cells <- which(residual_subset$Guides_collapsed_by_gene == "non-targeting")
n_control <- 1000
set.seed(2)
remove_idx <- sample(control_cells, length(control_cells) - n_control)
residual_subset <- residual_subset[-remove_idx]
cat("Residual matrix after subsampling controls:", dim(residual_subset), "\n")

# Check distribution of guide groups
print(table(residual_subset$Guides_collapsed_by_gene))

# Save cleaned dataset and gene clustering
saveRDS(residual_subset, file = paste0(output_dir, "/data/tiny_cleary_for_HMC.rds"))
saveRDS(clustering, file = paste0(output_dir, "/data/gene_clustering_for_HMC.rds"))

