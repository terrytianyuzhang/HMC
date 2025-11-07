# THIS FILE CAN BE DELETED, IT WAS FOR DEVELOPING THE TUTORIAL
# devtools::install_github("terrytianyuzhang/HMC/HMC_package")
rm(list = ls())
library(data.table)
library(HMC)
# library(glmnet)


# source('~/Documents/HMC/HMC_package/R/31_pre_processing_functions.R')
# source('~/Documents/HMC/HMC_package/R/51_convergence.R')
# Load data
output_dir <- "~/Documents/HMC/docs/"
gp_lasso_test_result <- readRDS(file = paste0(output_dir, "/data/HMC_convergence_control_STAT1_STAT2_grLasso.rds"))
gp_lasso_test_result

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
# Function to extract top nonzero genes per fold

combined_df <- visualize_convergence_top_genes(
  fold_data = gp_lasso_test_result$fold_data[1],
  top_n = 20,
  save_path = "top_genes_by_fold_and_type.pdf"
)
