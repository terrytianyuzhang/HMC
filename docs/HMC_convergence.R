# THIS FILE CAN BE DELETED, IT WAS FOR DEVELOPING THE TUTORIAL
# devtools::install_github("terrytianyuzhang/HMC/HMC_package")
rm(list = ls())
library(data.table)
library(HMC)
library(glmnet)

source('~/Documents/HMC/HMC_package/R/31_pre_processing_functions.R')
source('~/Documents/HMC/HMC_package/R/51_convergence.R')
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
gene_to_keep <- which(clustering$cluster_index %in% 1:41)
# gene_to_keep <- which(clustering$cluster_index %in% c(5))
control_subset <- control[, ..gene_to_keep]
STAT1_subset <- STAT1[, ..gene_to_keep]
STAT2_subset <- STAT2[, ..gene_to_keep]


test_result <- convergence_testing(
  control = control_subset,
  treatment1 = STAT1_subset,
  treatment2 = STAT2_subset,
  pca_method = "dense_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = TRUE
)



sapply(test_result$fold_data, function(x) x$variance)
control_scores <- sapply(test_result$fold_data, function(x) x$control_score, simplify = FALSE)
boxplot(control_scores, main = "Control Score Distribution per Fold",
        xlab = "Fold", ylab = "Control Score")

tr2_scores <- sapply(test_result$fold_data, function(x) x$tr2_score, simplify = FALSE)
boxplot(tr2_scores, main = "Tretment Score Distribution per Fold",
        xlab = "Fold", ylab = "Control Score")


test_result$p_value
collect_active_features(test_result)

saveRDS(test_result, file = paste0(output_dir, "data/HMC_convergence_control_STAT1_STAT2.rds"))

control_subset <- control[1:100, ]
psudo_control_subset <- control[101:200,]

test_result_null <- convergence_testing(
  control = control_subset,
  treatment1 = STAT1,
  treatment2 = psudo_control_subset,
  pca_method = "dense_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = TRUE
)
test_result_null$p_value
collect_active_features(test_result_null)


control_subset <- control[1:100, ]
psudo_control_subset <- control[101:200,]
psudo_control_subset[, 1:5] <- 10 * psudo_control_subset[, 1:5] 

test_result_null <- convergence_testing(
  control = control_subset,
  treatment1 =  STAT1,
  treatment2 = psudo_control_subset,
  pca_method = "dense_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = TRUE
)

test_result_null$p_value
collect_active_features(test_result_null)
psudo_control_subset[1:5, 1:10]
