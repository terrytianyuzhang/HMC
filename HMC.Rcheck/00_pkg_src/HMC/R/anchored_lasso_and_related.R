#' The function for nuisance parameter estimation in anchored_lasso_testing().
#'
#' @param nuisance_sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param nuisance_sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the discriminant direction. Default is logistic Lasso "lasso". Can also take value "lasso_no_truncation"
#' @param lasso_tuning_method Method for Lasso penalty hyperparameter tuning. Default is "min", the minimizer of cross-validation error; users can also use "1se" for more sparse solutions.
#' @param num_latent_factor The principle component that lasso coefficient anchors at. The default is PC1 = 1.
#' @param local_environment An environment for hyperparameters shared between folds.
#' @param verbose Print information to the console. Default is TRUE.
#' @export
#' 
#' @return A list of estimated nuisance quantities.
#' \item{estimate_leading_pc}{Leading principle components}
#' \item{estimate_mean_1}{Sample mean for group 1} 
#' \item{estimate_mean_2}{Sample mean for group 1}
#' \item{estimate_lasso_beta}{Logistic Lasso regression coefficients.}
#' \item{estimate_projection_direction}{Anchored projection direction. It is similar to PC1 when signal is weak but similar to estimate_optimal_direction when the signal is moderately large.}
#' \item{estimate_optimal_direction}{Discriminant direction.}

#' 
estimate_nuisance_parameter_lasso <- function(nuisance_sample_1,
                                              nuisance_sample_2,
                                              pca_method = "sparse_pca",
                                              mean_method = "lasso",
                                              lasso_tuning_method = "min",
                                              num_latent_factor = 1,
                                              local_environment = local_environment,
                                              verbose = TRUE){
  ###nuisance_sample_1 is required to be a data.frame
  
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
    if(local_environment$hyperparameter_shared_between_folds < 0){
      cv_result <- PMA::SPC.cv(sample_centered)
      # hyperparameter_shared_between_folds <<- max(1, 0.7*cv_result$bestsumabsv)
      assign("hyperparameter_shared_between_folds", cv_result$bestsumabsv, envir = local_environment)
    }
    # print(local_environment$hyperparameter_shared_between_folds)
    estimate_leading_pc <- PMA::SPC(sample_centered,
                               K = num_latent_factor,
                               sumabsv = local_environment$hyperparameter_shared_between_folds)$v
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
    
    non_zero_variance_feature <- (apply(nuisance_sample, 2 ,stats::var) != 0)
    nuisance_sample[,non_zero_variance_feature] <- scale(nuisance_sample[,non_zero_variance_feature])
    
    lasso_cv_result <- glmnet::cv.glmnet(y = nuisance_sample$class,
                                 x = as.matrix(nuisance_sample[, 1:feature_number]),
                                 family = 'binomial',
                                 type.measure = 'mse')
    
    if(lasso_tuning_method == 'min'){
      best_lambda_index <- which(lasso_cv_result$lambda == lasso_cv_result$lambda.min)
    }else if(lasso_tuning_method == '1se'){
      best_lambda_index <- which(lasso_cv_result$lambda == lasso_cv_result$lambda.1se)
    }else{
     stop("wrong lasso_tuning_method. Values should take either 'min' or '1se'.") 
    }
    
    estimate_lasso_beta <- lasso_cv_result$glmnet.fit$beta[, best_lambda_index]
    # hist((predict(lasso_cv_result$glmnet.fit, newx = as.matrix(nuisance_sample[,1:feature_number])))[,best_lambda_index])
    
    estimate_optimal_direction <- estimate_lasso_beta
    
    if(mean_method == 'lasso'){
      # print(norm(estimate_optimal_direction, type = '2'))
      if(norm(estimate_optimal_direction, type = '2') > 0.1*(effective_sample_size)^(-1/3)){ #
        if(verbose) message('Bingo!Found signal.')
        # estimate_optimal_direction <- sign(estimate_optimal_direction) * pmax(abs(estimate_optimal_direction) - (effective_sample_size)^(-1/3),0)
        estimate_projection_direction <- (effective_sample_size)^(-1/2)*estimate_leading_pc + (effective_sample_size^(1/2) - 1)* (effective_sample_size)^(-1/2)* estimate_optimal_direction
      }else{
        estimate_projection_direction <- estimate_leading_pc
      }
    }else if(mean_method == 'lasso_no_truncation'){
      if(norm(estimate_optimal_direction, type = '2') > 0){ #
        if(verbose) message('Bingo!Found signal.')
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
#' @param mean_method Methods used to estimate the discriminant direction. Default is logistic Lasso "lasso". Can also take value "lasso_no_truncation"
#' @export
#' 
#' @return A list of test statistics.
#' \item{influence_each_subject_1}{Test statistics for sample 1.}
#' \item{influence_each_subject_1}{Test statistics for sample 2.} 
#' \item{for_variance_each_subject_1}{Statistics for variance calculation, sample 1.}
#' \item{for_variance_each_subject_2}{Statistics for variance calculation, sample 2.} 
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
  
  ########
  
  cross_fitting_sample_2 <- as.data.frame(cross_fitting_sample_2)
  
  #####inner product statistics
  inner_product_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_projection_direction
  
  influence_each_subject_2 <- inner_product_2 
  
  
  ####these quantities are for variance estimation in the future
  for_variance_each_subject_1 <- influence_each_subject_1
  for_variance_each_subject_2 <- influence_each_subject_2
  
  # if(mean_method == 'lasso'){
  #   for_variance_each_subject_1 <- influence_each_subject_1
  #   for_variance_each_subject_2 <- influence_each_subject_2
  # }else if(mean_method == 'lasso_no_truncation'){
  #   estimate_leading_pc <- nuisance_collection$estimate_leading_pc
  #   for_variance_each_subject_1 <- as.matrix(cross_fitting_sample_1) %*% estimate_leading_pc
  #   for_variance_each_subject_2 <- as.matrix(cross_fitting_sample_2) %*% estimate_leading_pc
  # }
  
  
  return(list(influence_each_subject_1 = influence_each_subject_1,
              influence_each_subject_2 = influence_each_subject_2,
              for_variance_each_subject_1 = for_variance_each_subject_1,
              for_variance_each_subject_2 = for_variance_each_subject_2))
} 

#' Anchored test for two-sample mean comparison. 
#' 
#' @param sample_1 Group 1 sample. Each row is a subject and each column corresponds to a feature.
#' @param sample_2 Group 2 sample. Each row is a subject and each column corresponds to a feature.
#' @param pca_method Methods used to estimate principle component The default is "sparse_pca", using sparse PCA from package PMA. Other choices are "dense_pca"---the regular PCA; and "hard"--- hard-thresholding PCA, which also induces sparsity.
#' @param mean_method Methods used to estimate the discriminant direction. Default is logistic Lasso "lasso". Can also take value "lasso_no_truncation"
#' @param lasso_tuning_method Method for Lasso penalty hyperparameter tuning. Default is "min", the minimizer of cross-validation error; users can also use "1se" for more sparse solutions.
#' @param num_latent_factor The principle component that lasso coefficient anchors at. The default is PC1 = 1.
#' @param n_folds Number of splits when performing cross-fitting. The default is 5, if computational time allows, you can try to set it to 10.
#' @param verbose Print information to the console. Default is TRUE.
#' 
#' @export
#' @return A list of test statistics.
#' \item{test_statistics}{Test statistics. Each entry corresponds to the test result of one principle component.}
#' \item{standard_error}{Estimated standard error of test_statistics_before_studentization.} 
#' \item{test_statistics_before_studentization}{Similar to test_statistics but does not have variance = 1.}
#' \item{split_data}{Intermediate quantities needed for further assessment and interpretation of the test results.}
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
                                   lasso_tuning_method = "min",
                                   num_latent_factor = 1,
                                   n_folds = 5,
                                   verbose = TRUE){
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  split_data <- vector("list", n_folds)
  sample_1_split_index <- index_spliter(1:nrow(sample_1),
                                        n_folds = n_folds)
  sample_2_split_index <- index_spliter(1:nrow(sample_2),
                                        n_folds = n_folds)
  
  local_environment <- new.env(parent = emptyenv())  
  assign("hyperparameter_shared_between_folds", -1, envir = local_environment)
  # local_environment$hyperparameter_shared_between_folds <- -1
  
  for(split_index in 1:n_folds){
    if(verbose) message(paste0("processiong fold ", split_index))
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
                                                                                       lasso_tuning_method = lasso_tuning_method,
                                                                                       num_latent_factor = num_latent_factor,
                                                                                       local_environment = local_environment,
                                                                                       verbose = verbose)
    split_data[[split_index]]$influence_function_value <- evaluate_pca_lasso_plug_in(sample_1_cross,
                                                                                     sample_2_cross,
                                                                                     split_data[[split_index]]$nuisance_collection,
                                                                                     mean_method = mean_method)
    
    split_data[[split_index]]$test_statistic <- mean(split_data[[split_index]]$influence_function_value$influence_each_subject_1) - mean(split_data[[split_index]]$influence_function_value$influence_each_subject_2)
    split_data[[split_index]]$variance_sample_1 <- stats::var(split_data[[split_index]]$influence_function_value$for_variance_each_subject_1)
    split_data[[split_index]]$variance_sample_2 <- stats::var(split_data[[split_index]]$influence_function_value$for_variance_each_subject_2)
  }
  
  ####now combine the folds
  test_statistics_before_studentization <- 0
  variance_sample_1 <- variance_sample_2 <- 0
  
  for(split_index in 1:n_folds){
    inner_product_projection_direction <- crossprod(split_data[[1]]$nuisance_collection$estimate_projection_direction, 
                                                    split_data[[split_index]]$nuisance_collection$estimate_projection_direction)
    
    same_sign <- sign(inner_product_projection_direction)
    if(same_sign == 0){
      message('the projection directions are orthogonal')
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
