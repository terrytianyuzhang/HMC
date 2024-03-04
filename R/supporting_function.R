library(dplyr)

from_covariance_to_correlation <- function(simulation_covariance){
  simulation_correlation <- simulation_covariance
  for(row_index in 1:nrow(simulation_correlation)){
    for(column_index in 1:ncol(simulation_correlation)){
      simulation_correlation[row_index, column_index] <- simulation_covariance[row_index, column_index] / (sqrt(simulation_covariance[row_index, row_index]) * sqrt(simulation_covariance[column_index, column_index]))
    }
  }
  return(simulation_correlation)
}

generate_zero_inflated_normal <- function(sample_size, true_mean, normal_covariance,
                                          zero_prob){
  
  sample <- data.frame(MASS::mvrnorm(sample_size, 
                               mu = true_mean, 
                               Sigma = normal_covariance))
  num_feature <- length(true_mean)
  
  mask_vector <- rbinom(number_feature * sample_size, 1, zero_prob)
  mask_matrix <- matrix(mask_vector, nrow = sample_size)
  
  sample <- sample * mask_matrix
  
  inflated_covariance <- normal_covariance * zero_prob^2
  diag(inflated_covariance) <- diag(normal_covariance) * zero_prob
  
  return(list(sample = sample,
              inflated_covariance = inflated_covariance))
  
  ###test the above formula is correct when true_mean is zero
  # inflated_list <- generate_zero_inflated_normal(sample_size_1,
  #                                                true_mean_1, simulation_covariance, mask_probability)
  # print(inflated_list$inflated_covariance[1:5,1:5])
  # 
  # sample_matrix <- as.matrix(inflated_list$sample)
  # est_covariance <- t(sample_matrix) %*% (sample_matrix)/sample_size_1
  # est_covariance[1:5,1:5]
  # 
  # inflated_eigen_result <- eigen(inflated_list$inflated_covariance)
  # plot(inflated_eigen_result$vectors[,2])
  # inflated_eigen_result$values
}

# generate_zero_inflated_normal(sample_size_1, true_mean_1, simulation_covariance,
#                               0.5)

funny_lasso_plug_in <- function(cross_fitting_sample_1,
                                       cross_fitting_sample_2,
                                       nuisance_collection){
  
  # cross_fitting_sample_1 <<- cross_fitting_sample_1
  # cross_fitting_sample_2 <<- cross_fitting_sample_2
  # nuisance_collection <<- nuisance_collection
  
  estimate_projection_direction <- nuisance_collection$estimate_lasso_beta
  cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
  
  #####inner product statistics
  inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_projection_direction
  
  influence_each_subject_1 <- inner_product_1
  
  
  
  cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)
  
  #####inner product statistics
  inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_projection_direction
  
  influence_each_subject_2 <- inner_product_2 
  
  
  
  return(list(influence_each_subject_1 = influence_each_subject_1,
              influence_each_subject_2 = influence_each_subject_2))
}

funny_lasso <- function(sample_1, sample_2,
                        pca_method = "dense",
                        mean_method = "lasso",
                        n_folds = 5){
  set.seed(1)
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                        n_folds = n_folds)
  
  for(split_index in 1:n_folds){
    print(paste0("processiong fold ", split_index))
    sample_1_cross <- sample_1[sample_1_split_index[[split_index]], ]
    sample_1_nuisance <- sample_1[-sample_1_split_index[[split_index]], ]
    
    sample_2_cross <- sample_2[sample_2_split_index[[split_index]], ]
    sample_2_nuisance <- sample_2[-sample_2_split_index[[split_index]], ]
    
    split_data[[split_index]]$sample_1_cross_index <- sample_1_split_index[[split_index]]
    split_data[[split_index]]$sample_2_cross_index <- sample_2_split_index[[split_index]]
    
    split_data[[split_index]]$nuisance_collection <- estimate_nuisance_parameter_lasso(sample_1_nuisance, 
                                                                                       sample_2_nuisance,
                                                                                       pca_method = pca_method,
                                                                                       mean_method = mean_method)
    split_data[[split_index]]$influence_function_value <- funny_lasso_plug_in(sample_1_cross,
                                                                               sample_2_cross,
                                                                               split_data[[split_index]]$nuisance_collection)
    
    split_data[[split_index]]$test_statistic <- mean(split_data[[split_index]]$influence_function_value$influence_each_subject_1) - mean(split_data[[split_index]]$influence_function_value$influence_each_subject_2)
    split_data[[split_index]]$variance_sample_1 <- var(split_data[[split_index]]$influence_function_value$influence_each_subject_1)
    split_data[[split_index]]$variance_sample_2 <- var(split_data[[split_index]]$influence_function_value$influence_each_subject_2)
  }
  
  ####now combine the folds
  test_statistics_before_studentization <- 0
  variance_sample_1 <- variance_sample_2 <- 0
  
  for(split_index in 1:n_folds){
    inner_product_projection_direction <- crossprod(split_data[[1]]$nuisance_collection$estimate_projection_direction, 
                                                    split_data[[split_index]]$nuisance_collection$estimate_projection_direction)
    
    same_sign <- sign(inner_product_projection_direction)
    print(inner_product_projection_direction)
    if(same_sign == 0){
      print('the projection directions are orthogonal')
      same_sign <- 1
    }
    
    test_statistics_before_studentization <- test_statistics_before_studentization + same_sign * split_data[[split_index]]$test_statistic
    
    variance_sample_1 <- variance_sample_1 + split_data[[split_index]]$variance_sample_1
    variance_sample_2 <- variance_sample_2 + split_data[[split_index]]$variance_sample_2
  }
  
  test_statistics_before_studentization <- test_statistics_before_studentization/n_folds
  variance_sample_1 <- variance_sample_1/n_folds
  variance_sample_2 <- variance_sample_2/n_folds
  
  standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1 + 
                           (nrow(sample_2))^(-1) *variance_sample_2)
  test_statistics <- test_statistics_before_studentization/standard_error
  
  return(list(test_statistics = test_statistics,
              standard_error = standard_error,
              test_statistics_before_studentization = test_statistics_before_studentization,
              split_data = split_data))
}
extract_pc <- function(testing_result){
  n_folds <- length(testing_result$split_data)
  pc_each_split <- vector(mode = "list", length = n_folds)
  
  for(split_index in 1:n_folds){
    pc_each_split[[split_index]] <- testing_result$split_data[[split_index]]$nuisance_collection$estimate_leading_pc
  }
  
  return(pc_each_split)
}

extract_lasso_coef <- function(testing_result){
  n_folds <- length(testing_result$split_data)
  lasso_coef_each_split <- vector(mode = "list", length = n_folds)
  
  for(split_index in 1:n_folds){
    lasso_coef_each_split[[split_index]] <- testing_result$split_data[[split_index]]$nuisance_collection$estimate_lasso_beta
  }
  
  return(lasso_coef_each_split)
}

summarize_feature_name <- function(testing_result, method = 'majority voting'){
  lasso_coef_each_split <- extract_lasso_coef(testing_result)
  
  for(split_index in 1:length(lasso_coef_each_split)){
    lasso_coef_vector <- lasso_coef_each_split[[split_index]]
    sorted_lasso_coef_vector <- sort(abs(lasso_coef_vector), decreasing = T)
    temp_gene_symbols <- names(sorted_lasso_coef_vector[which(abs(sorted_lasso_coef_vector)>0)])
    if(split_index == 1){
      gene_symbols <- temp_gene_symbols
    }else{
      if(method == 'union'){
        gene_symbols <- union(gene_symbols, temp_gene_symbols)
      }else if (method %in% c('majority voting', 'intersection')){
        gene_symbols <- c(gene_symbols, temp_gene_symbols)
      }
    }
  }
  
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  if(method == 'majority voting'){
    gene_counts <- table(gene_symbols)
    gene_symbols <- names(gene_counts[gene_counts >= ceiling(length(lasso_coef_each_split)/2)])  
  }else if(method == 'intersection'){
    gene_counts <- table(gene_symbols)
    gene_symbols <- names(gene_counts[gene_counts >= length(lasso_coef_each_split)])  
  }
  
  
  return(gene_symbols)
}

summarize_pc_name <- function(testing_result, latent_fator_index = 1, method = 'majority voting'){
  pc_each_split <- extract_pc(testing_result)
  
  for(split_index in 1:length(pc_each_split)){
    pc_vector <- pc_each_split[[split_index]][, latent_fator_index]
    sorted_pc_vector <- sort(abs(pc_vector), decreasing = T)
    temp_gene_symbols <- names(sorted_pc_vector[which(abs(sorted_pc_vector)>0)])
    if(split_index == 1){
      gene_symbols <- temp_gene_symbols
    }else{
      if(method == 'union'){
        gene_symbols <- union(gene_symbols, temp_gene_symbols)
      }else if (method %in% c('majority voting', 'intersection')){
        gene_symbols <- c(gene_symbols, temp_gene_symbols)
      }
    }
  }
  
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  if(method == 'majority voting'){
    gene_counts <- table(gene_symbols)
    gene_symbols <- names(gene_counts[gene_counts >= ceiling(length(pc_each_split)/2)])  
  }else if(method == 'intersection'){
    gene_counts <- table(gene_symbols)
    gene_symbols <- names(gene_counts[gene_counts >= length(pc_each_split)])  
  }
  return(gene_symbols)
}

combine_data_frame_of_diff_column <- function(df1, df2){
  
  # df1 <- data.frame(B = 4:6, A = 1:3)
  # df2 <- data.frame(A = 7:8, C = 9:10, D = 11:12)
  
  # Determine the full set of column names
  if(nrow(df1) > 0){
  all_columns <- union(colnames(df1), colnames(df2))
  
  # Create a function to add missing columns with NA
  add_missing_columns <- function(df, all_columns) {
    missing_columns <- setdiff(all_columns, colnames(df))
    new_columns <- as.data.frame(matrix(NA, ncol = length(missing_columns), nrow = nrow(df)))
    colnames(new_columns) <- missing_columns
    bind_cols(df, new_columns)
  }
  
  # Add missing columns with NA to both data frames
  df1_filled <- add_missing_columns(df1, all_columns)
  df2_filled <- add_missing_columns(df2, all_columns)
  
  # Combine data frames and fill missing columns with NA
  combined_df <- bind_rows(df1_filled, df2_filled)
  
  # Print the combined data frame
  return(combined_df)
  }else{
    combined_df <- bind_rows(df1, df2)
    
    # Print the combined data frame
    return(combined_df)
  }
  
}
# gene_symbols <- summarize_feature_name(simple_lasso_test_result)
# gene_ontology_analysis_pc <- function(pc_each_split, dim_to_analysis = 20,
#                                       keyType = "SYMBOL", ont = "BP"){
#   
#   library(clusterProfiler)
#   library(AnnotationDbi)
#   library(org.Hs.eg.db)
#   
#   num_latent_factor <- ncol(pc_each_split[[1]])
#   GO_result <- vector('list', num_latent_factor)
#   
#   for(latent_factor_index in 1:num_latent_factor){
#     print(paste0("performing gene ontology analysis for latent factor", latent_factor_index))
#     
#     pc_vector <- pc_each_split[[1]][,latent_factor_index]
#     gene_symbols <- names(sort(abs(pc_vector), decreasing = T)[1:dim_to_analysis])
#     
#     result <- enrichGO(gene_symbols,
#                        OrgDb = org.Hs.eg.db, keyType = keyType, ont = ont)
#     # result_df <- as.data.frame(result)
#     goplot(result)
#     plot(barplot(result, showCategory =20))
#     
#     GO_result[[latent_factor_index]] <- result
#   }
#   return(GO_result)
# }
# 
# go_result <- gene_ontology_analysis_pc(pc_each_split)
# 


