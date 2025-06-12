library(data.table)
library(HMC)

# Load data
output_dir <- "~/Documents/HMC/docs/"
residual_subset <- readRDS(paste0(output_dir, "/data/tiny_cleary_for_HMC.rds"))
clustering <- readRDS(paste0(output_dir, "/data/gene_clustering_for_HMC.rds"))

# Subset to control and treatment groups, then remove guide labels
control <- residual_subset[
  Guides_collapsed_by_gene == "non-targeting",
  !"Guides_collapsed_by_gene"
]
STAT1 <- residual_subset[
  Guides_collapsed_by_gene == "STAT1",
  !"Guides_collapsed_by_gene"
]

STAT2 <- residual_subset[
  Guides_collapsed_by_gene == "STAT2",
  !"Guides_collapsed_by_gene"
]


clustering <- readRDS(paste0(output_dir, "/data/gene_clustering_for_HMC.rds"))
gene_to_keep <- which(clustering$cluster_index %in% c(31, 34, 37))
control_subset <- control[, ..gene_to_keep]
STAT1_subset <- STAT1[, ..gene_to_keep]
STAT2_subset <- STAT2[, ..gene_to_keep]

set.seed(123)
test_result <- convergence_testing(
  control = control_subset,
  treatment1 = STAT1_subset,
  treatment2 = STAT2_subset,
  pca_method = "sparse_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  verbose = FALSE
)
debug(convergence_testing)
test_result <- mean_comparison_anchor(
  control = control_subset,
  treatment = STAT1_subset,
  pca_method = "sparse_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = FALSE
)
