# library(PMA)
# library(irlba)
# library(glmnet)
# library(MASS)

#' Split the sample index into n_folds many groups so that we can perform cross-fitting
#'
#' @param array Sample index. Usually just an array from 1 to the number of samples in one group.
#' @param n_folds Number of splits
#' @export
#' 
index_spliter <- function(array, n_folds = 5){
  
  # array <- 1:99
  
  # Calculate the length of each part
  part_length <- length(array) %/% n_folds
  
  # Create an empty list to store the parts
  parts <- vector("list", n_folds)
  
  # Randomly shuffle the array
  shuffled_array <- sample(array)
  
  # Split the shuffled array into parts
  for (fold_index in 1:n_folds) {
    start_index <- (fold_index - 1) * part_length + 1
    
    if(fold_index < n_folds){
      end_index <- fold_index * part_length
    }else{
      end_index <- length(array)
    }
    
    parts[[fold_index]] <- shuffled_array[start_index:end_index]
  }
  
  return(parts)
}

#' The function for nuisance parameter estimation in simple_pc_testing() and debiased_pc_testing().
#'
#' @param nuisance_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the mean vector. Default is sample mean "naive". There is also a hard-thresholding sparse estiamtor "hard".
#' @param num_latent_factor Number of principle to be estimated/tested. Default is 1.
#' 
#' @export
#' 
estimate_nuisance_pc <- function(nuisance_sample_1,
                                 nuisance_sample_2 = NULL,
                                 pca_method = "sparse_pca",
                                 mean_method = "naive",
                                 num_latent_factor = 1){
  
  nuisance_sample_1 <<- nuisance_sample_1
  nuisance_sample_2 <<- nuisance_sample_2
  
  feature_number <- ncol(nuisance_sample_1)
  
  nuisance_sample_1 <- as.data.frame(nuisance_sample_1)
  nuisance_sample_size_1 <- nrow(nuisance_sample_1)
  
  ###estimate mean vector 
  estimate_mean_1 <- colMeans(nuisance_sample_1)
  # estimate_mean_1[abs(estimate_mean_1) < mean_threshold_1] <- 0
  if(mean_method == 'hard'){
    variance_feature_wise_centered <- apply(nuisance_sample_1, 2, var)
    mean_threshold <- 2*sqrt(variance_feature_wise_centered)*sqrt(2*log(2*feature_number*nuisance_sample_size_1)/nuisance_sample_size_1)
    estimate_mean_1[abs(estimate_mean_1) < mean_threshold] <- 0
  }
  
  centered_sample_1 <- sweep(nuisance_sample_1, 2, estimate_mean_1)
  
  #######
  
  if(!is.null(nuisance_sample_2)){
    
    nuisance_sample_2 <- as.data.frame(nuisance_sample_2)
    nuisance_sample_size_2 <- nrow(nuisance_sample_2)
    
    ###estimate mean vector
    estimate_mean_2 <- colMeans(nuisance_sample_2)
    # estimate_mean_2[abs(estimate_mean_2) < mean_threshold_2] <- 0
    
    if(mean_method == 'hard'){
      variance_feature_wise_centered <- apply(nuisance_sample_2, 2, var)
      mean_threshold <- 2*sqrt(variance_feature_wise_centered)*sqrt(2*log(2*feature_number*nuisance_sample_size_2)/nuisance_sample_size_2)
      estimate_mean_2[abs(estimate_mean_2) < mean_threshold] <- 0
    }
    
    centered_sample_2 <- sweep(nuisance_sample_2, 2, estimate_mean_2)
    
    sample_centered <- as.matrix(rbind(centered_sample_1, centered_sample_2))
  }else{
    estimate_mean_2 <- NULL
    sample_centered <- as.matrix(centered_sample_1)
  }
  
  ########
  
  nuisance_sample_size <- nrow(sample_centered)
  
  if(pca_method == "sparse_pca"){
    ###estimate sparse principle component

    if(hyperparameter_shared_between_folds < 0){
      # cv_result <- SPC.cv(sample_centered,
      #                     sumabsvs = seq(1, sqrt(floor(feature_number)), length = 10))
      # hyperparameter_shared_between_folds <<- cv_result$bestsumabsv
      cv_result <- PMA::SPC.cv(sample_centered)
      hyperparameter_shared_between_folds <<- max(1, 0.7*cv_result$bestsumabsv) ##I USED THIS FOR LUPAS REAL-DATA ANALYSIS
      hyperparameter_shared_between_folds <<- cv_result$bestsumabsv
    }
    print(hyperparameter_shared_between_folds)
    estimate_leading_pc <- PMA::SPC(sample_centered,
                               K = num_latent_factor,
                               sumabsv = hyperparameter_shared_between_folds)$v
    # plot(estimate_leading_pc[,1])
    # plot(estimate_leading_pc[,2])
    estimate_leading_pc <- estimate_leading_pc/apply(X = estimate_leading_pc, MARGIN = 2, FUN = norm, type = '2')
    # estimate_leading_pc <- estimate_leading_pc/norm(estimate_leading_pc, type = '2')
  }else if(pca_method == "dense_pca"){
    estimate_leading_pc <- irlba::irlba(sample_centered, nv = num_latent_factor)$v
    # plot(estimate_leading_pc[,1])
    # plot(estimate_leading_pc[,2])
  }else if(pca_method == "hard"){
    estimate_leading_pc <- dense_leading_pc <- irlba::irlba(sample_centered, nv = num_latent_factor)$v
    pc_threshold <-nuisance_sample_size^(-1/3)
    estimate_leading_pc[abs(estimate_leading_pc) < pc_threshold] <- 0
    
    for(latent_index in 1:num_latent_factor){
      estimate_vector <- estimate_leading_pc[, latent_index]
      if(norm(estimate_vector, type = '2') == 0){
        estimate_leading_pc[, latent_index] <- dense_leading_pc[, latent_index]
      }else{
        estimate_leading_pc[, latent_index] <- estimate_vector/norm(estimate_vector, type = '2')
      }
    }
    
    # if(norm(estimate_leading_pc, type = '2') == 0){
    #   estimate_leading_pc <- dense_leading_pc
    # }else{
    #   estimate_leading_pc <- estimate_leading_pc/norm(estimate_leading_pc, type = '2')
    # }
  }
  colnames(estimate_leading_pc) <- paste0('pc', 1:ncol(estimate_leading_pc))
  rownames(estimate_leading_pc) <- colnames(nuisance_sample_1)
  
  ##estimate eigenvalue
  data_times_pc <- sample_centered %*% estimate_leading_pc
  estimate_eigenvalue <- (apply(X = data_times_pc, MARGIN = 2, 
                               FUN = norm, type = '2'))^2/nuisance_sample_size
  # estimate_eigenvalue <- norm(data_times_pc, type = '2')^2/nuisance_sample_size
  
  
  ##estimate noise variance
  estimate_covariance_matrix_diagonal<- apply(sample_centered, 2, function(x) sum(x^2))/nuisance_sample_size
  estimate_covariance_matrix_trace <- sum(estimate_covariance_matrix_diagonal)
  estimate_noise_variance <- (estimate_covariance_matrix_trace - sum(estimate_eigenvalue))/(feature_number - num_latent_factor)
  
  # ##estimate the transformed eigenvector
  # sigma_times_pc <- t(sample_centered) %*% data_times_pc / nuisance_sample_size
  # 
  
  return(list(estimate_leading_pc = estimate_leading_pc,
              estimate_mean_1 = estimate_mean_1,
              estimate_mean_2 = estimate_mean_2,
              estimate_eigenvalue = estimate_eigenvalue,
              estimate_noise_variance = estimate_noise_variance))
              # estimate_eigenvalue = 1.1,
              # estimate_noise_variance = 0.1,
              # sigma_times_pc = sigma_times_pc))
}

#' Calculate the test statistics on the left-out samples. Called in simple_pc_testing().
#'
#' @param cross_fitting_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param cross_fitting_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_collection A collection of nuisance quantities estimated using "nuisance" samples. It is the output of estimate_nuisance_pc().
#' 
#' @export
#' 
#' 
evaluate_pca_plug_in <- function(cross_fitting_sample_1,
                                 cross_fitting_sample_2 = NULL,
                                 nuisance_collection){
  cross_fitting_sample_1 <<- cross_fitting_sample_1
  cross_fitting_sample_2 <<- cross_fitting_sample_2
  nuisance_collection <<- nuisance_collection
  
  cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
  
  
  estimate_leading_pc <- nuisance_collection$estimate_leading_pc
  
  #####inner product statistics
  inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
  
  influence_each_subject_1 <- inner_product_1
  colnames(influence_each_subject_1) <- paste0('pc', 1:ncol(estimate_leading_pc))
  
  if(!is.null(cross_fitting_sample_2)){
    
    cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)
    
    #####inner product statistics
    inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_leading_pc
    
    influence_each_subject_2 <- inner_product_2 
    colnames(influence_each_subject_2) <- paste0('pc', 1:ncol(estimate_leading_pc))
  }else{
    influence_each_subject_2 <- NULL
  }
  
  
  return(list(influence_each_subject_1 = influence_each_subject_1,
              influence_each_subject_2 = influence_each_subject_2))
} 

#' Calculate the test statistics on the left-out samples. Called in debiased_pc_testing().
#'
#' @param cross_fitting_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param cross_fitting_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_collection A collection of nuisance quantities estimated using "nuisance" samples. It is the output of estimate_nuisance_pc().
#' @param num_latent_factor Number of principle components to be considered.
#' @export
evaluate_influence_function_multi_factor <- function(cross_fitting_sample_1,
                                                     cross_fitting_sample_2 = NULL,
                                                     nuisance_collection,
                                                     num_latent_factor = 1){
  
  cross_fitting_sample_1 <<- cross_fitting_sample_1
  cross_fitting_sample_2 <<- cross_fitting_sample_2
  nuisance_collection <<- nuisance_collection
  
  estimate_leading_pc <- nuisance_collection$estimate_leading_pc
  estimate_eigenvalue <- nuisance_collection$estimate_eigenvalue
  estimate_noise_variance <- nuisance_collection$estimate_noise_variance
  
  if(is.null(cross_fitting_sample_2)){ # if this is a one-sample testing problem
    
    cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
    
    #####inner product statistics
    estimate_mean_1 <- nuisance_collection$estimate_mean_1
    inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
    
    ####influence function of the eigenvector
    centered_sample_1 <- as.matrix(sweep(cross_fitting_sample_1, 2, estimate_mean_1)) 
    inner_product_centered_1 <- centered_sample_1 %*% estimate_leading_pc
    
    calculate_mu_inverse_centered <- function(latent_factor_index){
      ###
      pc_interaction <- (estimate_eigenvalue[latent_factor_index] - estimate_eigenvalue)^(-1)
      pc_interaction[which(pc_interaction == Inf | pc_interaction == -Inf)] <- 0
      
      pc_noise_interaction <- (rep(estimate_eigenvalue[latent_factor_index] - estimate_noise_variance, num_latent_factor))^(-1)
      
      coefficients_in_influence_function <- pc_interaction - pc_noise_interaction
      pc_matrix <- t( coefficients_in_influence_function* t(estimate_leading_pc))
      
      mu_times_pc <- matrix(estimate_mean_1, nrow = 1) %*% pc_matrix
      mu_times_pc_times_individual <- mu_times_pc %*% t(inner_product_centered_1)
      
      mu_times_centered <- (estimate_eigenvalue[latent_factor_index] - estimate_noise_variance)^(-1) * centered_sample_1 %*% estimate_mean_1
      
      mu_inverse_centered <- mu_times_centered + t(mu_times_pc_times_individual)
      
      return(mu_inverse_centered)
    }
    
    influence_each_subject_1_all <- NULL
    for(latent_factor_index in 1:num_latent_factor){
      mu_inverse_centered <- calculate_mu_inverse_centered(latent_factor_index = latent_factor_index)
      correction_each_individual <- mu_inverse_centered * inner_product_centered_1[,latent_factor_index]
      influence_each_subject_1 <- inner_product_1[, latent_factor_index] + correction_each_individual
      influence_each_subject_1_all <- cbind(influence_each_subject_1_all, influence_each_subject_1)
    }
    
    influence_each_subject_2_all <- NULL
    colnames(influence_each_subject_1_all) <- paste0('pc', 1:num_latent_factor)
    
    return(list(influence_each_subject_1 = influence_each_subject_1_all,
                influence_each_subject_2 = influence_each_subject_2_all))
    # mu_times_centered <- centered_sample_1 %*% estimate_mean_1
    # mu_times_pc <- as.numeric(crossprod(estimate_leading_pc, estimate_mean_1))
    # 
    # influence_function_eigenvector <- (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * (mu_times_centered - mu_times_pc * inner_product_centered_1)*inner_product_centered_1
    # 
    # ####the omitted component
    # # ommited_component <- (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * (crossprod(estimate_mean_1, sigma_times_pc) - mu_times_pc * nuisance_collection$estimate_eigenvalue)
    # ommited_component <- 0
  }
  
  if(!is.null(cross_fitting_sample_2)){
    
    cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
    cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)

    #####inner product statistics, also used in the simple statistics
    estimate_mean_1 <- nuisance_collection$estimate_mean_1
    inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
    
    estimate_mean_2 <- nuisance_collection$estimate_mean_2
    inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_leading_pc
    
    
    ####influence function of the eigenvector
    centered_sample_1 <- as.matrix(sweep(cross_fitting_sample_1, 2, estimate_mean_1)) 
    centered_1_times_pc <- inner_product_centered_1 <- centered_sample_1 %*% estimate_leading_pc
    
    centered_sample_2 <- as.matrix(sweep(cross_fitting_sample_2, 2, estimate_mean_2))
    centered_2_times_pc <- inner_product_centered_2 <- centered_sample_2 %*% estimate_leading_pc
    
    mean_diff <- estimate_mean_1 - estimate_mean_2
    two_sample_weight <- nrow(cross_fitting_sample_1)/(nrow(cross_fitting_sample_1) + nrow(cross_fitting_sample_2))
    
    
    calculate_mu_inverse_centered <- function(latent_factor_index){
      
      ###the following lines are preparing the part of pseudo-inverse matrix that is related to the latent factors
      pc_interaction <- (estimate_eigenvalue[latent_factor_index] - estimate_eigenvalue)^(-1)
      pc_interaction[which(pc_interaction == Inf | pc_interaction == -Inf)] <- 0
      
      pc_noise_interaction <- (rep(estimate_eigenvalue[latent_factor_index] - estimate_noise_variance, num_latent_factor))^(-1)
      
      coefficients_in_influence_function <- pc_interaction - pc_noise_interaction
      pc_matrix <- t( coefficients_in_influence_function* t(estimate_leading_pc))
      
      #######
      
      mu_times_pc <- matrix(mean_diff, nrow = 1) %*% pc_matrix
      
      mu_times_pc_times_individual <- mu_times_pc %*% t(inner_product_centered_1)
      mu_times_centered <- (estimate_eigenvalue[latent_factor_index] - estimate_noise_variance)^(-1) * centered_sample_1 %*% mean_diff
      mu_inverse_centered_1 <- mu_times_centered + t(mu_times_pc_times_individual)
      
      mu_times_pc_times_individual <- mu_times_pc %*% t(inner_product_centered_2)
      mu_times_centered <- (estimate_eigenvalue[latent_factor_index] - estimate_noise_variance)^(-1) * centered_sample_2 %*% mean_diff
      mu_inverse_centered_2 <- mu_times_centered + t(mu_times_pc_times_individual)
      
      return(list(mu_inverse_centered_1 = mu_inverse_centered_1,
                  mu_inverse_centered_2 = mu_inverse_centered_2))
    }
    
    influence_eigenvector_each_subject_1_all <- influence_eigenvector_each_subject_2_all <- NULL
    for_variance_subject_1_all <- for_variance_subject_2_all <- NULL
    
    for(latent_factor_index in 1:num_latent_factor){
      mu_inverse_centered <- calculate_mu_inverse_centered(latent_factor_index = latent_factor_index)
      mu_inverse_centered_1 <- mu_inverse_centered$mu_inverse_centered_1
      mu_inverse_centered_2 <- mu_inverse_centered$mu_inverse_centered_2
      
      influence_eigenvector_each_subject_1 <- two_sample_weight * mu_inverse_centered_1 * inner_product_centered_1[,latent_factor_index]
      influence_eigenvector_each_subject_1_all <- cbind(influence_eigenvector_each_subject_1_all, influence_eigenvector_each_subject_1)
      
      for_variance_subject_1 <- inner_product_1[, latent_factor_index] + influence_eigenvector_each_subject_1
      for_variance_subject_1_all <- cbind(for_variance_subject_1_all, for_variance_subject_1)
      
      influence_eigenvector_each_subject_2 <- (1 - two_sample_weight) * mu_inverse_centered_2 * inner_product_centered_2[,latent_factor_index]
      influence_eigenvector_each_subject_2_all <- cbind(influence_eigenvector_each_subject_2_all, influence_eigenvector_each_subject_2)
      
      for_variance_subject_2 <- inner_product_2[, latent_factor_index] - influence_eigenvector_each_subject_2
      for_variance_subject_2_all <- cbind(for_variance_subject_2_all, for_variance_subject_2)
    }
    
    
    # centered_1_times_meandiff <- centered_sample_1 %*% mean_diff
    # centered_2_times_meandiff <- centered_sample_2 %*% mean_diff
    # 
    # meandiff_times_pc <- as.numeric(crossprod(mean_diff, estimate_leading_pc))
    # two_sample_weight <- nrow(cross_fitting_sample_1)/(nrow(cross_fitting_sample_1) + nrow(cross_fitting_sample_2))
    # 
    # influence_eigenvector_each_subject_1 <- (centered_1_times_meandiff - meandiff_times_pc * centered_1_times_pc) * centered_1_times_pc
    # influence_eigenvector_each_subject_1 <- two_sample_weight * (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * influence_eigenvector_each_subject_1
    # 
    # influence_eigenvector_each_subject_2 <- (centered_2_times_meandiff - meandiff_times_pc * centered_2_times_pc) * centered_2_times_pc
    # influence_eigenvector_each_subject_2 <- (1 - two_sample_weight) * (nuisance_collection$estimate_eigenvalue - nuisance_collection$estimate_noise_variance)^(-1) * influence_eigenvector_each_subject_2
    # 
    # for_variance_subject_1 <- inner_product_1 + influence_eigenvector_each_subject_1
    # for_variance_subject_2 <- inner_product_2 - influence_eigenvector_each_subject_2
    # 
    return(list(inner_product_1 = inner_product_1,
                inner_product_2 = inner_product_2,
                influence_eigenvector_each_subject_1 = influence_eigenvector_each_subject_1_all,
                influence_eigenvector_each_subject_2 = influence_eigenvector_each_subject_2_all,
                for_variance_subject_1 = for_variance_subject_1_all,
                for_variance_subject_2 = for_variance_subject_2_all))
  }
  
  
} 

#' Simple plug-in test for two-sample mean comparison. 
#' 
#' @param sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the mean vector. Default is sample mean "naive". There is also a hard-thresholding sparse estiamtor "hard".
#' @param num_latent_factor Number of principle to be estimated/tested. Default is 1.
#' @param n_folds Number of splits when performing cross-fitting. The default is 5, if computational time allows, you can try to set it to 10.
#' @export
#' @examples 
#' sample_size_1 <- sample_size_2 <- 300
#' 
#' true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
#' true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)
#' pc1 <- c(rep(1, 10), rep(0, 90))
#' pc1 <- pc1/norm(pc1, type = '2')
#' 
#' simulation_covariance <- 10 * pc1 %*% t(pc1)
#' simulation_covariance <- simulation_covariance + diag(1, 100)
#' 
#' sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
#'                                mu = true_mean_1,
#'                                Sigma = simulation_covariance))
#'  sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
#'                                mu = true_mean_2,
#'                                Sigma = simulation_covariance))
#'  result <- simple_pc_testing(sample_1, sample_2)
#'  result$test_statistics 
#'  ##these are test statistics. Each one of them corresponds to one PC.
#'  summarize_pc_name(result, latent_fator_index = 1) #shows which features contribute to PC1
#'  extract_pc(result) # extract the estimated leading PCs.
#'  
simple_pc_testing <- function(sample_1,
                              sample_2 = NULL,
                              pca_method = "sparse_pca",
                              mean_method = "naive",
                              num_latent_factor = 1,
                              n_folds = 5){
  
  ###split the samples into n_folds splits
  set.seed(1)
  sample_1 <- as.data.frame(sample_1)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  
  if(!is.null(sample_2)){
    sample_2 <- as.data.frame(sample_2)
    sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                          n_folds = n_folds)
  }
  
  hyperparameter_shared_between_folds <<- -1
  ###process each split
  for(split_index in 1:n_folds){
    
    print(paste0("processiong fold ", split_index))
    
    sample_1_cross <- sample_1[sample_1_split_index[[split_index]], ]
    sample_1_nuisance <- sample_1[-sample_1_split_index[[split_index]], ]
    split_data[[split_index]]$sample_1_cross_index <- sample_1_split_index[[split_index]]
    
    if(!is.null(sample_2)){
      split_data[[split_index]]$sample_2_cross_index <- sample_2_split_index[[split_index]]
      sample_2_cross <- sample_2[sample_2_split_index[[split_index]], ]
      sample_2_nuisance <- sample_2[-sample_2_split_index[[split_index]], ]
    }else{
      sample_2_cross <- sample_2_nuisance <- NULL
    }
    

    split_data[[split_index]]$nuisance_collection <- estimate_nuisance_pc(sample_1_nuisance, 
                                                                          sample_2_nuisance,
                                                                          pca_method = pca_method,
                                                                          mean_method = mean_method,
                                                                          num_latent_factor = num_latent_factor)
    split_data[[split_index]]$influence_function_value <- evaluate_pca_plug_in(sample_1_cross,
                                                                               sample_2_cross,
                                                                               split_data[[split_index]]$nuisance_collection)
    
    ##############
    
    sample_1_individual <- split_data[[split_index]]$influence_function_value$influence_each_subject_1
    split_data[[split_index]]$variance_sample_1 <- apply(X = sample_1_individual, MARGIN = 2, FUN = var)
    
    if(!is.null(sample_2)){
      sample_2_individual <- split_data[[split_index]]$influence_function_value$influence_each_subject_2
      split_data[[split_index]]$variance_sample_2 <- apply(X = sample_2_individual, MARGIN = 2, FUN = var)
      
      split_data[[split_index]]$test_statistic <- apply(X = sample_1_individual, MARGIN = 2, FUN = mean) - apply(X = sample_2_individual, MARGIN = 2, FUN = mean)
    }else{
      split_data[[split_index]]$test_statistic <- apply(X = sample_1_individual, MARGIN = 2, FUN = mean)
    }
  }
  
  ####now combine the folds
  combine_1_factor_over_splits <- function(latent_factor_index){
    test_statistics_before_studentization <- 0
    variance_sample_1 <- variance_sample_2 <- 0
    
    for(split_index in 1:n_folds){
      inner_product_projection_direction <- crossprod(split_data[[1]]$nuisance_collection$estimate_leading_pc[,latent_factor_index], 
                                                      split_data[[split_index]]$nuisance_collection$estimate_leading_pc[,latent_factor_index])
      
      same_sign <- sign(inner_product_projection_direction)
      print(inner_product_projection_direction)
      if(same_sign == 0){
        print('the projection directions are orthogonal')
        same_sign <- 1
      }
      
      test_statistics_before_studentization <- test_statistics_before_studentization + same_sign * split_data[[split_index]]$test_statistic[latent_factor_index]
      
      variance_sample_1 <- variance_sample_1 + split_data[[split_index]]$variance_sample_1[latent_factor_index]
      
      if(!is.null(sample_2)){
        variance_sample_2 <- variance_sample_2 + split_data[[split_index]]$variance_sample_2[latent_factor_index]
      }
    }
    
    test_statistics_before_studentization <- test_statistics_before_studentization/n_folds
    variance_sample_1 <- variance_sample_1/n_folds
    
    if(!is.null(sample_2)) variance_sample_2 <- variance_sample_2/n_folds
    
    if(!is.null(sample_2)){
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1 + 
                               (nrow(sample_2))^(-1) *variance_sample_2)
    }else{
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1)
    }
    
    test_statistics <- test_statistics_before_studentization/standard_error
    
    tuple_one_latent_factor <- c(test_statistics_before_studentization, standard_error, test_statistics)
    return(tuple_one_latent_factor)
  }
  
  test_statistics_before_studentization <- standard_error <- test_statistics <- rep(NULL,num_latent_factor)
  for(latent_factor_index in 1:num_latent_factor){
    tuple_one_latent_factor <- combine_1_factor_over_splits(latent_factor_index)
    
    test_statistics_before_studentization[latent_factor_index] <- tuple_one_latent_factor[1]
    standard_error[latent_factor_index] <- tuple_one_latent_factor[2]
    test_statistics[latent_factor_index] <- tuple_one_latent_factor[3]
  }
  
  
  return(list(test_statistics = test_statistics,
              standard_error = standard_error,
              test_statistics_before_studentization = test_statistics_before_studentization,
              split_data = split_data))
}

#' Debiased one-step test for two-sample mean comparison. A small p-value tells us not only there is difference in the mean vectors, but can also indicates which principle component the difference aligns with.
#' 
#' @param sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the mean vector. Default is sample mean "naive". There is also a hard-thresholding sparse estiamtor "hard".
#' @param num_latent_factor Number of principle to be estimated/tested. Default is 1.
#' @param n_folds Number of splits when performing cross-fitting. The default is 5, if computational time allows, you can try to set it to 10.
#' @export
#' @examples 
#' sample_size_1 <- sample_size_2 <- 300
#' 
#' true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
#' true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)
#' pc1 <- c(rep(1, 10), rep(0, 90))
#' pc1 <- pc1/norm(pc1, type = '2')
#' 
#' simulation_covariance <- 10 * pc1 %*% t(pc1)
#' simulation_covariance <- simulation_covariance + diag(1, 100)
#' 
#' sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
#'                                mu = true_mean_1,
#'                                Sigma = simulation_covariance))
#'  sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
#'                                mu = true_mean_2,
#'                                Sigma = simulation_covariance))
#'  result <- debiased_pc_testing(sample_1, sample_2)
#'  result$test_statistics
#'  ##these are test statistics. Each one of them corresponds to one PC.
#'  summarize_pc_name(result, latent_fator_index = 1) #shows which features contribute to PC1
#'  extract_pc(result) # extract the estimated leading PCs.
#'  
#'  
debiased_pc_testing <- function(sample_1,
                                sample_2 = NULL,
                                pca_method = "sparse_pca",
                                mean_method = "naive",
                                num_latent_factor = 1,
                                n_folds = 5){
  
  ###split the samples into n_folds splits
  set.seed(1)
  sample_1 <- as.data.frame(sample_1)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  
  if(!is.null(sample_2)){
    sample_2 <- as.data.frame(sample_2)
    sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                          n_folds = n_folds)
  }
  
  hyperparameter_shared_between_folds <<- -1
  ###process each split
  for(split_index in 1:n_folds){
    print(paste0("processiong fold ", split_index))
    
    sample_1_cross <- sample_1[sample_1_split_index[[split_index]], ]
    sample_1_nuisance <- sample_1[-sample_1_split_index[[split_index]], ]
    split_data[[split_index]]$sample_1_cross_index <- sample_1_split_index[[split_index]]
    
    if(!is.null(sample_2)){
      split_data[[split_index]]$sample_2_cross_index <- sample_2_split_index[[split_index]]
      sample_2_cross <- sample_2[sample_2_split_index[[split_index]], ]
      sample_2_nuisance <- sample_2[-sample_2_split_index[[split_index]], ]
    }else{
      sample_2_cross <- sample_2_nuisance <- NULL
    }
    
    
    split_data[[split_index]]$nuisance_collection <- estimate_nuisance_pc(sample_1_nuisance, 
                                                                          sample_2_nuisance,
                                                                          pca_method = pca_method,
                                                                          mean_method = mean_method,
                                                                          num_latent_factor = num_latent_factor)
    split_data[[split_index]]$influence_function_value <- evaluate_influence_function_multi_factor(sample_1_cross,
                                                                                                   sample_2_cross,
                                                                                                   nuisance_collection = split_data[[split_index]]$nuisance_collection,
                                                                                                   num_latent_factor = num_latent_factor)
    # testing data, should remove later
    
    # to_look_new <- split_data[[split_index]]$influence_function_value
    # to_look_new_multi <- split_data[[split_index]]$influence_function_value
    # to_look_old <- evaluate_influence_function(sample_1_cross, sample_2_cross, nuisance_collection = split_data[[split_index]]$nuisance_collection)
    # summary(to_look_new$for_variance_subject_1 - to_look_old$for_variance_subject_1)
    # to_look_new$influence_eigenvector_each_subject_1[1:10]
    # to_look_new_multi$influence_eigenvector_each_subject_1[1:10, 1]
    ##############
    
    if(!is.null(sample_2)){
      
      inner_product_1 <- split_data[[split_index]]$influence_function_value$inner_product_1
      inner_product_2 <- split_data[[split_index]]$influence_function_value$inner_product_2
      split_data[[split_index]]$test_statistic <- apply(X = inner_product_1, MARGIN = 2, FUN = mean) - apply(X = inner_product_2, MARGIN = 2, FUN = mean)
      #debias
      sample_1_correction <- split_data[[split_index]]$influence_function_value$influence_eigenvector_each_subject_1
      sample_2_correction <- split_data[[split_index]]$influence_function_value$influence_eigenvector_each_subject_2
      
      split_data[[split_index]]$test_statistic <- split_data[[split_index]]$test_statistic + apply(X = sample_1_correction, MARGIN = 2, FUN = mean) + apply(X = sample_2_correction, MARGIN = 2, FUN = mean)
      
      ##variance
      split_data[[split_index]]$variance_sample_1 <- apply(X = split_data[[split_index]]$influence_function_value$for_variance_subject_1, MARGIN = 2, FUN = var)
      split_data[[split_index]]$variance_sample_2 <- apply(X = split_data[[split_index]]$influence_function_value$for_variance_subject_2, MARGIN = 2, FUN = var)
      
    }else{
      sample_1_individual <- split_data[[split_index]]$influence_function_value$influence_each_subject_1
      split_data[[split_index]]$variance_sample_1 <- apply(X = sample_1_individual, MARGIN = 2, FUN = var)
      split_data[[split_index]]$test_statistic <- apply(X = sample_1_individual, MARGIN = 2, FUN = mean)
    }
    
}
  
  
  ####now combine the folds
  combine_1_factor_over_splits <- function(latent_factor_index){
    test_statistics_before_studentization <- 0
    variance_sample_1 <- variance_sample_2 <- 0
    
    for(split_index in 1:n_folds){
      inner_product_projection_direction <- crossprod(split_data[[1]]$nuisance_collection$estimate_leading_pc[,latent_factor_index], 
                                                      split_data[[split_index]]$nuisance_collection$estimate_leading_pc[,latent_factor_index])
      
      same_sign <- sign(inner_product_projection_direction)
      print(inner_product_projection_direction)
      if(same_sign == 0){
        print('the projection directions are orthogonal')
        same_sign <- 1
      }
      
      test_statistics_before_studentization <- test_statistics_before_studentization + same_sign * split_data[[split_index]]$test_statistic[latent_factor_index]
      
      variance_sample_1 <- variance_sample_1 + split_data[[split_index]]$variance_sample_1[latent_factor_index]
      
      if(!is.null(sample_2)){
        variance_sample_2 <- variance_sample_2 + split_data[[split_index]]$variance_sample_2[latent_factor_index]
      }
    }
    
    test_statistics_before_studentization <- test_statistics_before_studentization/n_folds
    variance_sample_1 <- variance_sample_1/n_folds
    
    if(!is.null(sample_2)) variance_sample_2 <- variance_sample_2/n_folds
    
    if(!is.null(sample_2)){
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1 + 
                               (nrow(sample_2))^(-1) *variance_sample_2)
    }else{
      standard_error <- sqrt((nrow(sample_1))^(-1) * variance_sample_1)
    }
    
    test_statistics <- test_statistics_before_studentization/standard_error
    
    tuple_one_latent_factor <- c(test_statistics_before_studentization, standard_error, test_statistics)
    return(tuple_one_latent_factor)
  }
  
  test_statistics_before_studentization <- standard_error <- test_statistics <- rep(NULL,num_latent_factor)
  for(latent_factor_index in 1:num_latent_factor){
    tuple_one_latent_factor <- combine_1_factor_over_splits(latent_factor_index)
    
    test_statistics_before_studentization[latent_factor_index] <- tuple_one_latent_factor[1]
    standard_error[latent_factor_index] <- tuple_one_latent_factor[2]
    test_statistics[latent_factor_index] <- tuple_one_latent_factor[3]
  }
  
  
  return(list(test_statistics = test_statistics,
              standard_error = standard_error,
              test_statistics_before_studentization = test_statistics_before_studentization,
              split_data = split_data))
}

#' The function for nuisance parameter estimation in anchored_lasso_testing().
#'
#' @param nuisance_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the discriminant direction. Default is logistic Lasso "lasso". 
#' @param num_latent_factor The principle component that lasso coefficient anchors at. The default is PC1 = 1.
#' @export
#' 
estimate_nuisance_parameter_lasso <- function(nuisance_sample_1,
                                              nuisance_sample_2,
                                              pca_method = "sparse_pca",
                                              mean_method = "lasso",
                                              num_latent_factor = 1){
  ###nuisance_sample_1 is required to be a data.frame
  nuisance_sample_1 <<- nuisance_sample_1
  nuisance_sample_2 <<- nuisance_sample_2
  
  feature_number <- ncol(nuisance_sample_1)
  
  nuisance_sample_1 <- as.data.frame(nuisance_sample_1)
  nuisance_sample_size_1 <- nrow(nuisance_sample_1)
  
  
  ###estimate mean vector 
  estimate_mean_1 <- colMeans(nuisance_sample_1)

  centered_sample_1 <- sweep(nuisance_sample_1, 2, estimate_mean_1)
  
  #######
  
  nuisance_sample_2 <- as.data.frame(nuisance_sample_2)
  nuisance_sample_size_2 <- nrow(nuisance_sample_2)
  
  
  ###estimate mean vector
  estimate_mean_2 <- colMeans(nuisance_sample_2)

  centered_sample_2 <- sweep(nuisance_sample_2, 2, estimate_mean_2)
  
  sample_centered <- as.matrix(rbind(centered_sample_1, centered_sample_2))
  
  ########
  
  nuisance_sample_size <- nrow(sample_centered)
  effective_sample_size <- 2*min(nuisance_sample_size_1,nuisance_sample_size_2)
  
  if(pca_method == "sparse_pca"){
    ###estimate sparse principle component
    if(hyperparameter_shared_between_folds < 0){
      cv_result <- PMA::SPC.cv(sample_centered)
      hyperparameter_shared_between_folds <<- max(1, 0.7*cv_result$bestsumabsv)
    }
    print(hyperparameter_shared_between_folds)
    estimate_leading_pc <- PMA::SPC(sample_centered,
                               K = num_latent_factor,
                               sumabsv = hyperparameter_shared_between_folds)$v
    estimate_leading_pc <- estimate_leading_pc[, num_latent_factor]
    # cv_result <- SPC.cv(sample_centered,
    #                     sumabsvs = seq(1, sqrt(floor(feature_number)), length = 5))
    # estimate_leading_pc <- SPC(sample_centered,
    #                            K = 1,
    #                            sumabsv = cv_result$bestsumabsv)$v
    estimate_leading_pc <- estimate_leading_pc/norm(estimate_leading_pc, type = '2')
  }else if(pca_method == "hard"){
    estimate_leading_pc <- dense_leading_pc <- irlba::irlba(sample_centered, nv = num_latent_factor)$v
    pc_threshold <-nuisance_sample_size^(-1/3)
    estimate_leading_pc[abs(estimate_leading_pc) < pc_threshold] <- 0
    
    for(latent_index in 1:num_latent_factor){
      estimate_vector <- estimate_leading_pc[, latent_index]
      if(norm(estimate_vector, type = '2') == 0){
        estimate_leading_pc[, latent_index] <- dense_leading_pc[, latent_index]
      }else{
        estimate_leading_pc[, latent_index] <- estimate_vector/norm(estimate_vector, type = '2')
      }
    }
    estimate_leading_pc <- estimate_leading_pc[, num_latent_factor]
  }else if(pca_method == "dense_pca"){
    estimate_leading_pc <- irlba::irlba(sample_centered, nv = num_latent_factor)$v
    estimate_leading_pc <- estimate_leading_pc[, num_latent_factor]
    # estimate_leading_pc <- estimate_leading_pc/norm(estimate_leading_pc, type = '2')
  }
  
  
  if(mean_method %in% c("lasso", "lasso_no_truncation")){
    ####estimate lasso beta
    # sample_1$class <- 1
    # sample_2$class <- -1
    # nuisance_sample <- rbind(sample_1, sample_2)
    nuisance_sample_1$class <- 1
    nuisance_sample_2$class <- 0
    nuisance_sample <- rbind(nuisance_sample_1, nuisance_sample_2)
    
    non_zero_variance_feature <- (apply(nuisance_sample, 2 ,var) != 0)
    nuisance_sample[,non_zero_variance_feature] <- scale(nuisance_sample[,non_zero_variance_feature])
    
    lasso_cv_result <- glmnet::cv.glmnet(y = nuisance_sample$class,
                                 x = as.matrix(nuisance_sample[, 1:feature_number]),
                                 family = 'binomial',
                                 type.measure = 'mse')
    # results <- glm(class~X16,nuisance_sample, family = 'binomial')
    # summary(results)
    # t.test(nuisance_sample_1[,16], nuisance_sample_2[,16], var.equal = FALSE)
    best_lambda_index <- which(lasso_cv_result$lambda == lasso_cv_result$lambda.1se)
    estimate_lasso_beta <- lasso_cv_result$glmnet.fit$beta[, best_lambda_index]
    # hist((predict(lasso_cv_result$glmnet.fit, newx = as.matrix(nuisance_sample[,1:feature_number])))[,best_lambda_index])
    
    estimate_optimal_direction <- estimate_lasso_beta
    
    if(mean_method == 'lasso'){
      print(norm(estimate_optimal_direction, type = '2'))
      if(norm(estimate_optimal_direction, type = '2') > 0.1*(effective_sample_size)^(-1/3)){ #
        print('bingo')
        # estimate_optimal_direction <- sign(estimate_optimal_direction) * pmax(abs(estimate_optimal_direction) - (effective_sample_size)^(-1/3),0)
        estimate_projection_direction <- (effective_sample_size)^(-1/2)*estimate_leading_pc + (effective_sample_size^(1/2) - 1)* (effective_sample_size)^(-1/2)* estimate_optimal_direction
      }else{
        estimate_projection_direction <- estimate_leading_pc
      }
    }else if(mean_method == 'lasso_no_truncation'){
      if(norm(estimate_optimal_direction, type = '2') > 0){ #
        print('bingo')
        # estimate_optimal_direction <- sign(estimate_optimal_direction) * pmax(abs(estimate_optimal_direction) - (effective_sample_size)^(-1/3),0)
        estimate_projection_direction <- estimate_leading_pc + effective_sample_size^(1/3) * estimate_optimal_direction
      }else{
        estimate_projection_direction <- estimate_leading_pc
      }
    }
  }
  # estimate_projection_direction <- (log(log(effective_sample_size)))^(-1)*estimate_leading_pc + (log(log(effective_sample_size)) - 1)* (log(log(effective_sample_size)))^(-1)* estimate_optimal_direction
  
  estimate_projection_direction <- estimate_projection_direction/norm(estimate_projection_direction, type = '2')
  
  
  return(list(estimate_leading_pc = estimate_leading_pc,
              estimate_mean_1 = estimate_mean_1,
              estimate_mean_2 = estimate_mean_2,
              estimate_lasso_beta = estimate_lasso_beta,
              estimate_projection_direction = estimate_projection_direction,
              estimate_optimal_direction = estimate_optimal_direction))
}

#' Calculate the test statistics on the left-out samples. Called in anchored_lasso_testing().
#'
#' @param cross_fitting_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param cross_fitting_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_collection A collection of nuisance quantities estimated using "nuisance" samples. It is the output of estimate_nuisance_pc().
#' @param mean_method Methods used to estimate the discriminant direction. Default is logistic Lasso "lasso". 
#' @export
#' 
#' 
evaluate_pca_lasso_plug_in <- function(cross_fitting_sample_1,
                                       cross_fitting_sample_2,
                                       nuisance_collection,
                                       mean_method = 'lasso'){
  
  # cross_fitting_sample_1 <<- cross_fitting_sample_1
  # cross_fitting_sample_2 <<- cross_fitting_sample_2
  # nuisance_collection <<- nuisance_collection
  
  estimate_projection_direction <- nuisance_collection$estimate_projection_direction
  cross_fitting_sample_1 <- as.data.frame(cross_fitting_sample_1)
  
  #####inner product statistics
  inner_product_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_projection_direction
  
  influence_each_subject_1 <- inner_product_1
  
  cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)
  
  #####inner product statistics
  inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_projection_direction
  
  influence_each_subject_2 <- inner_product_2 
  
  if(mean_method == 'lasso'){
    for_variance_each_subject_1 <- influence_each_subject_1
    for_variance_each_subject_2 <- influence_each_subject_2
  }else if(mean_method == 'lasso_no_truncation'){
    estimate_leading_pc <- nuisance_collection$estimate_leading_pc
    for_variance_each_subject_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
    for_variance_each_subject_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_leading_pc
  }
  
  
  return(list(influence_each_subject_1 = influence_each_subject_1,
              influence_each_subject_2 = influence_each_subject_2,
              for_variance_each_subject_1 = for_variance_each_subject_1,
              for_variance_each_subject_2 = for_variance_each_subject_2))
} 

#' Simple plug-in test for two-sample mean comparison. 
#' 
#' @param sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the mean vector. Default is sample mean "naive". There is also a hard-thresholding sparse estiamtor "hard".
#' @param num_latent_factor The principle component that lasso coefficient anchors at. The default is PC1 = 1.
#' @param n_folds Number of splits when performing cross-fitting. The default is 5, if computational time allows, you can try to set it to 10.
#' @export
#' @examples 
#' sample_size_1 <- sample_size_2 <- 300
#' true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
#' true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)
#' 
#' sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
#'                                mu = true_mean_1,
#'                                Sigma = diag(1, 100)))
#'  sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
#'                                mu = true_mean_2,
#'                                Sigma = diag(1, 100)))
#'  result <- anchored_lasso_testing(sample_1, sample_2)
#'  result$test_statistics
#'  ##the test statistic. It should follow normal(0,1) when there is no difference between the groups.
#'  summarize_feature_name(result) 
#'  #summarize which features contribute to discriminant vectors (i.e. logistic lasso)
#'  extract_pc(result) # extract the estimated discriminant coefficients
#'  


anchored_lasso_testing <- function(sample_1, sample_2,
                                   pca_method = "sparse_pca",
                                   mean_method = "lasso",
                                   num_latent_factor = 1,
                                   n_folds = 5){
  set.seed(1)
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                        n_folds = n_folds)
  
  hyperparameter_shared_between_folds <<- -1
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
                                                                                       mean_method = mean_method,
                                                                                       num_latent_factor = num_latent_factor)
    split_data[[split_index]]$influence_function_value <- evaluate_pca_lasso_plug_in(sample_1_cross,
                                                                                     sample_2_cross,
                                                                                     split_data[[split_index]]$nuisance_collection,
                                                                                     mean_method = mean_method)
    
    split_data[[split_index]]$test_statistic <- mean(split_data[[split_index]]$influence_function_value$influence_each_subject_1) - mean(split_data[[split_index]]$influence_function_value$influence_each_subject_2)
    split_data[[split_index]]$variance_sample_1 <- var(split_data[[split_index]]$influence_function_value$for_variance_each_subject_1)
    split_data[[split_index]]$variance_sample_2 <- var(split_data[[split_index]]$influence_function_value$for_variance_each_subject_2)
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

# 
# correct_coverage <- function(point_estimate, standard_error, true_parameter){
#   true_parameter <- as.numeric(abs(true_parameter))
#   point_estimate <- abs(point_estimate)
#   lower_bound <- point_estimate - 1.96 * standard_error
#   upper_bound <- point_estimate + 1.96 * standard_error
#   
#   correct_cover <- (lower_bound < true_parameter) & (upper_bound > true_parameter)
#   return(correct_cover)
# }
# 
# 
# naive_pca <- function(sample_1,
#                       sample_2 = NULL,
#                       pca_method = "sparse_pca"){
#   set.seed(1)
#   # sample_1_split_1_index <- sample(1:nrow(sample_1), floor(nrow(sample_1)/ 2))
#   # sample_1_split_1 <- sample_1[sample_1_split_1_index, ]
#   # sample_1_split_2 <- sample_1[-sample_1_split_1_index, ]
#   
#   naive_fitting <- vector("list", 1)
#   # crossffiting_two_fold[[1]]$sample_1_nuisance <- sample_1_split_1
#   # crossffiting_two_fold[[1]]$sample_1_crossfitting <- sample_1_split_2
#   # crossffiting_two_fold[[2]]$sample_1_nuisance <- sample_1_split_2
#   # crossffiting_two_fold[[2]]$sample_1_crossfitting <- sample_1_split_1
#   # 
#   # if(!is.null(sample_2)){
#   #   sample_2_split_1_index <- sample(1:nrow(sample_2), floor(nrow(sample_2)/ 2))
#   #   sample_2_split_1 <- sample_2[sample_2_split_1_index, ]
#   #   sample_2_split_2 <- sample_2[-sample_2_split_1_index, ]
#   #   
#   #   crossffiting_two_fold[[1]]$sample_2_nuisance <- sample_2_split_1
#   #   crossffiting_two_fold[[1]]$sample_2_crossfitting <- sample_2_split_2
#   #   crossffiting_two_fold[[2]]$sample_2_nuisance <- sample_2_split_2
#   #   crossffiting_two_fold[[2]]$sample_2_crossfitting <- sample_2_split_1
#   # }else{
#   #   crossffiting_two_fold[[1]]$sample_2_nuisance <- NULL
#   #   crossffiting_two_fold[[1]]$sample_2_crossfitting <- NULL
#   #   crossffiting_two_fold[[2]]$sample_2_nuisance <- NULL
#   #   crossffiting_two_fold[[2]]$sample_2_crossfitting <- NULL
#   # }
#   
#   # for(split_index in 1:2){
#   #   crossffiting_two_fold[[split_index]]$nuisance_collection <- estimate_nuisance_parameter(crossffiting_two_fold[[split_index]]$sample_1_nuisance, 
#   #                                                                                           crossffiting_two_fold[[split_index]]$sample_2_nuisance,
#   #                                                                                           method = pca_method)
#   #   crossffiting_two_fold[[split_index]]$influence_function_value <- evaluate_pca_plug_in(crossffiting_two_fold[[split_index]]$sample_1_crossfitting,
#   #                                                                                         crossffiting_two_fold[[split_index]]$sample_2_crossfitting,
#   #                                                                                         crossffiting_two_fold[[split_index]]$nuisance_collection)
#   #   if(!is.null(sample_2)){
#   #     crossffiting_two_fold[[split_index]]$test_statistic <- mean(crossffiting_two_fold[[split_index]]$influence_function_value$influence_each_subject_1) - mean(crossffiting_two_fold[[split_index]]$influence_function_value$influence_each_subject_2)
#   #   }else{
#   #     crossffiting_two_fold[[split_index]]$test_statistic <- mean(crossffiting_two_fold[[split_index]]$influence_function_value$influence_each_subject_1)
#   #   }
#   # }
#   
#   naive_fitting[[1]]$nuisance_collection <- estimate_nuisance_parameter(sample_1, 
#                                                                         sample_2,
#                                                                         pca_method = pca_method,
#                                                                         mean_method = mean_method)
#   naive_fitting[[1]]$influence_function_value <- evaluate_pca_plug_in(sample_1,
#                                                                       sample_2,
#                                                                       naive_fitting[[1]]$nuisance_collection)
#   if(!is.null(sample_2)){
#     naive_fitting[[1]]$test_statistic <- mean(naive_fitting[[1]]$influence_function_value$influence_each_subject_1) - mean(naive_fitting[[1]]$influence_function_value$influence_each_subject_2)
#   }else{
#     naive_fitting[[1]]$test_statistic <- mean(naive_fitting[[1]]$influence_function_value$influence_each_subject_1)
#   }
#   
#   ####now combine the two folds
#   
#   # inner_product_leading_pc <- crossprod(crossffiting_two_fold[[1]]$nuisance_collection$estimate_leading_pc, 
#   #                                       crossffiting_two_fold[[2]]$nuisance_collection$estimate_leading_pc)
#   # same_sign <- sign(inner_product_leading_pc)
#   # print(inner_product_leading_pc)
#   # if(same_sign == 0){
#   #   print('the principle componenets are orthogonal')
#   #   same_sign <- 1
#   # }
#   
#   test_statistics_before_studentization <- naive_fitting[[1]]$test_statistic
#   
#   if(!is.null(sample_2)){
#     standard_error <-  sqrt(nrow(sample_1)^(-1)*var(naive_fitting[[1]]$influence_function_value$influence_each_subject_1)+
#                             nrow(sample_2)^(-1)*var(naive_fitting[[1]]$influence_function_value$influence_each_subject_2))
#     
#   }else{
#     standard_error <-  sqrt(nrow(sample_1)^(-1)*var(naive_fitting[[1]]$influence_function_value$influence_each_subject_1))
#   }
#   
#   
#   ####calculate variance of the test statistics
#   # nuisance_collection_not_split <- estimate_nuisance_parameter(sample_1, sample_2, method = pca_method)
#   # influence_function_value_not_split <- evaluate_pca_plug_in(sample_1, sample_2, nuisance_collection_not_split)
#   # 
#   # if(!is.null(sample_2)){
#   #   standard_error <- sqrt(nrow(sample_1)^(-1)*var(influence_function_value_not_split$influence_each_subject_1) + nrow(sample_2)^(-1)*var(influence_function_value_not_split$influence_each_subject_2))
#   # }else{
#   #   standard_error <- sqrt(nrow(sample_1)^(-1)*var(influence_function_value_not_split$influence_each_subject_1))
#   # }
#   
#   test_statistics <- test_statistics_before_studentization/standard_error
#   
#   return(list(test_statistics = test_statistics,
#               standard_error = standard_error,
#               test_statistics_before_studentization = test_statistics_before_studentization,
#               naive_fitting = naive_fitting))
#   
# }
# 
# 
# 
# verify_formula <- function(x = 1){
#   ###verify if the formula is correct
#   #direct calculation
#   # library(MASS)
#   # target_matrix <- 0.6*diag(number_feature) - simulation_covariance
#   # direct_calculation <-   ginv(target_matrix)
#   # summary(array(target_matrix %*% direct_calculation %*% target_matrix - target_matrix))
#   # ###use the spike-model formula
#   # part1 <- (0.6- 0.1)^(-1)*diag(number_feature)
#   # part2 <- (0.6 - 1.1)^(-1)*pc1%*% t(pc1)
#   # part3 <- ((0.6 - 0.1)^(-1)*pc1%*% t(pc1) +(0.6 - 0.1)^(-1)*pc2%*% t(pc2) )
#   # trick_calculation <- part1 + part2 - part3
#   # summary(array(direct_calculation - trick_calculation))
#   # summary(array(target_matrix %*% trick_calculation %*% target_matrix - target_matrix))
#   
#   # verify the performant formula
#   # latent_factor_index <- 1
#   # pc1 <- estimate_leading_pc[,1]
#   # pc2 <- estimate_leading_pc[,2]
#   # 
#   # estimate_covariance <- (estimate_eigenvalue[1] - estimate_noise_variance) * pc1 %*% t(pc1) +(estimate_eigenvalue[2] - estimate_noise_variance) * pc2 %*% t(pc2) + estimate_noise_variance * diag(number_feature)
#   # direct_calculation <- estimate_mean_1 %*%ginv(estimate_eigenvalue[latent_factor_index]*diag(number_feature) - estimate_covariance) %*% t(centered_sample_1)
#   # 
#   # summary(array(direct_calculation-t(mu_inverse_centered)))
# }
# 
