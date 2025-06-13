#' Identify candidate genes based on model coefficients and grouping
#'
#' Filters model-selected features to identify candidate genes for further analysis,
#' optionally aggregating by group if using group lasso.
#'
#' @param beta_1_con A named numeric vector of estimated coefficients from the first-stage model.
#' @param group A vector indicating group membership for each gene (feature). Must be the same length as the number of columns in \code{control_train}.
#' @param classifier_method Character string indicating the classifier method. Must be either \code{"lasso"} or \code{"group_lasso"}.
#' @param control_train A matrix or data frame representing the control group, used to infer gene names if needed.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{candidate_genes}{A character vector of gene names with non-zero coefficients or belonging to selected groups.}
#'   \item{group_subset}{A named vector of group labels corresponding to \code{candidate_genes}, or \code{NULL} if not applicable.}
#' }
#'
#' @export

process_candidate_genes <- function(beta_1_con, group, classifier_method, control_train) {
  # Ensure the group vector has names corresponding to gene identifiers
  if (!is.null(group)) {
    if (is.null(names(group))) {
      if (length(group) != ncol(control_train)) {
        stop("Length of `group` does not match the number of columns in `control_train`.")
      }
      names(group) <- colnames(control_train)
    }
  }
  
  if (classifier_method == "group_lasso") {
    # Count nonzero genes per group
    group_nonzero_counts <- table(group[names(beta_1_con[abs(beta_1_con) > 1e-10])])
    
    selected_groups <- names(group_nonzero_counts[group_nonzero_counts >= 5])  # Keep groups with ≥5 nonzero genes
    
    # Extract genes from selected groups
    candidate_genes <- names(beta_1_con)[group[names(beta_1_con)] %in% selected_groups]
    
    # Subset the group vector for the second call to fit_lasso
    group_subset <- group[candidate_genes]
  } else {
    # Default case: select nonzero coefficients directly
    candidate_genes <- names(beta_1_con[abs(beta_1_con) > 1e-10])
    group_subset <- NULL  # No need to subset
  }
  
  return(list(candidate_genes = candidate_genes, group_subset = group_subset))
}

#' Process a single cross-validation fold for projected mean comparison
#'
#' Performs one fold of a two-stage high-dimensional mean comparison procedure.
#' This includes sparse feature selection via lasso/group-lasso, projection direction
#' adjustment via principal components, and computation of test statistics and scores.
#'
#' @param i Integer. Fold index.
#' @param control A matrix or data frame of control samples.
#' @param treatment1 A matrix or data frame of first treatment group samples.
#' @param treatment2 A matrix or data frame of second treatment group samples.
#' @param control_split_index A list of row indices used to split \code{control} into train/test sets.
#' @param tr1_split_index A list of row indices used to split \code{treatment1}.
#' @param tr2_split_index A list of row indices used to split \code{treatment2}.
#' @param pca_method Character. Method to estimate the principal component. Either \code{"dense_pca"} or \code{"sparse_pca"}.
#' @param classifier_method Character. Classification method. Either \code{"lasso"} or \code{"group_lasso"}.
#' @param lambda_type Character. Type of regularization parameter used in model fitting. Either \code{"lambda.min"} or \code{"lambda.1se"}.
#' @param group Optional grouping vector for group-lasso. Should match the number of columns in \code{control}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{Estimated difference in projected means between control and treatment2.}
#'   \item{variance}{Estimated variance of the test statistic.}
#'   \item{control_score}{Projected scores for control test data.}
#'   \item{tr2_score}{Projected scores for treatment2 test data.}
#'   \item{proj_direction}{Final projection direction vector used for scoring.}
#'   \item{first_beta}{First-stage model coefficients (control vs treatment1).}
#'   \item{final_beta}{Second-stage model coefficients (control vs treatment2).}
#'   \item{second_pc}{Principal component vector from control data.}
#' }
#'
#' @export

process_fold <- function(i, control, treatment1, treatment2, 
                         control_split_index, tr1_split_index, tr2_split_index,
                         pca_method, classifier_method, lambda_type, group, verbose) {
  if (verbose) message(paste0("Processing fold ", i))
  
  # Extract training and testing datasets
  control_test <- control[control_split_index[[i]], ]
  control_train <- control[-control_split_index[[i]], ]
  tr1_train <- treatment1[-tr1_split_index[[i]], ]
  tr2_test <- treatment2[tr2_split_index[[i]], ]
  tr2_train <- treatment2[-tr2_split_index[[i]], ]
  
  # Fit Lasso model for control vs treatment1
  beta_1_con <- tryCatch({
    fit_lasso(control_train, tr1_train, lambda_type, classifier_method, group)
  }, error = function(e) {
    message("LASSO failed for first training set: ", e$message)
    return(NULL)
  })
  
  gene_info <- process_candidate_genes(beta_1_con, group, classifier_method, control_train)
  candidate_genes <- gene_info$candidate_genes
  group_subset <- gene_info$group_subset
  
  # beta_1_con <- fit_lasso(control_train, tr1_train)
  if (length(candidate_genes) == 0) {
    if (verbose) message("Skipping fold, treatment 1 is the same as control.")
    
    return(list(
      statistic = NA, 
      variance = NA,
      control_score = NA,
      tr2_score = NA,
      proj_direction = NA,
      first_beta = NA,
      final_beta = NA,
      second_pc = NA
    ))
  }
  
  control_train_interest <- control_train[, candidate_genes, drop = FALSE]
  tr2_train_interest <- tr2_train[, candidate_genes, drop = FALSE]
  
  # Fit Lasso for control vs treatment2
  beta_2_con <- tryCatch({
    fit_lasso(control_train_interest, tr2_train_interest, lambda_type, classifier_method, group_subset)
  }, error = function(e) {
    message("LASSO failed for second training set: ", e$message)
    return(NULL)
  })
  
  # Estimate Principal Component and Projection Direction
  pc_2_con <- estimate_leading_pc(control_train_interest)
  
  # ===================================================
  # Compute adjusted projection direction
  # ===================================================
  n_effect2 <- 2 * min(nrow(control_train_interest), nrow(tr2_train_interest))  # Effective training size
  a_n <- n_effect2^(1/3)  # Scaling factor based on effective sample size
  
  
  # Adjust projection direction using Lasso coefficients and normalize
  proj_direction <- (pc_2_con + a_n * beta_2_con)  # Adjustment
  proj_direction <- proj_direction / norm(proj_direction, type = "2")  # Normalize
  
  # ===================================================
  # Compute test statistic, scores, and variance
  # ===================================================
  # Extract relevant gene data from control and treatment groups
  control_matrix <- as.matrix(control_test[, candidate_genes], ncol = length(candidate_genes))
  tr2_matrix <- as.matrix(tr2_test[, candidate_genes], ncol = length(candidate_genes))
  
  # Compute scores for control and treatment groups
  control_score <- control_matrix %*% proj_direction
  tr2_score <- tr2_matrix %*% proj_direction
  
  # Compute the test statistic as the difference of means
  T_stat <- mean(control_score) - mean(tr2_score)
  
  # Get sample sizes
  n_x <- nrow(control_test)
  n_z <- nrow(tr2_test)
  
  # Compute the variance of the test statistic
  variance <- (var(control_score) * (n_x - 1) / (n_x^2)) + (var(tr2_score) * (n_z - 1) / (n_z^2))
  
  # Take notes for each split
  return(list(
    statistic = T_stat,  # T_stat is the mean instead of the one with a standard normal distribution
    variance = variance,
    control_score = control_score,
    tr2_score = tr2_score,
    proj_direction = proj_direction,
    first_beta = beta_1_con,
    final_beta = beta_2_con,
    second_pc = pc_2_con
  ))
}

#' Combine results from multiple folds to compute final test statistic
#'
#' Aggregates per-fold statistics from a high-dimensional projection test procedure
#' to compute a single test statistic and two-sided p-value. Projection directions
#' are aligned and padded to ensure consistency across folds.
#'
#' @param fold_data A list of length \code{n_folds}, where each element contains results from a single fold,
#' typically returned by \code{process_fold()}.
#' @param n_folds Integer. Number of folds used in the procedure.
#' @param verbose Logical. Whether to print messages during fold combination (e.g., orthogonal direction warnings).
#'
#' @return A list containing:
#' \describe{
#'   \item{p_value}{Two-sided p-value of the aggregated test.}
#'   \item{test_statistic}{Final test statistic computed as the mean of signed fold statistics normalized by standard error.}
#'   \item{fold_data}{Original \code{fold_data} list, possibly modified to include padded projection directions.}
#' }
#'
#' @export


combine_folds <- function(fold_data, n_folds, verbose = FALSE) {
  numerator_test_statistic<-0
  denominator_variance<-0
  
  missing_proj_directions <- sapply(fold_data[1:n_folds], function(x) any(is.na(x$proj_direction))) 
  
  if (all(missing_proj_directions == TRUE)) {  # All splits degenerate
    return(list(
      p_value = 1,
      test_statistic = 0,
      # standard_error = NA,
      # numerator_test_statistic = 0,
      fold_data = fold_data
    ))
  }
  
  # Identify successful folds (non-degenerate cases)
  valid_folds <- which(!missing_proj_directions)
  first_valid_fold <- valid_folds[1]  # Select the first non-degenerate case as the baseline
  folds_effect_num <- length(valid_folds)

  # ===================================================
  # Padding for each proj_direction so that I can combine them later
  # ===================================================
  
  # Determine the max length among all effective proj_direction and collect unique names
  padding_set <- unique(unlist(lapply(valid_folds, function(i) names(fold_data[[i]]$proj_direction))))
  max_length <- length(padding_set)
  
  # Pad proj_direction for each successful fold
  for (i in valid_folds) {
    new_proj_direction <- numeric(max_length)
    names(new_proj_direction) <- padding_set
    new_proj_direction[names(fold_data[[i]]$proj_direction)] <- fold_data[[i]]$proj_direction
    fold_data[[i]]$proj_direction <- new_proj_direction
  }
  
  # ===================================================
  # compute test statistics
  # ===================================================
  
  for (i in valid_folds) {
    
    denominator_variance <- denominator_variance + fold_data[[i]]$variance
    projection_sign_match <- sign(crossprod(fold_data[[first_valid_fold]]$proj_direction, fold_data[[i]]$proj_direction))
    
    if (projection_sign_match == 0) {
      if (verbose) message("The projection directions are orthogonal")
      projection_sign_match <- 1
    }
    
    numerator_test_statistic <- numerator_test_statistic + projection_sign_match * fold_data[[i]]$statistic
  }
  
  numerator_test_statistic <- numerator_test_statistic / n_folds
  standard_error <- sqrt(denominator_variance / (folds_effect_num^2))
  test_statistic <- numerator_test_statistic / standard_error
  p_value <- 2 * pnorm(-abs(test_statistic))
  
  return(list(
    p_value = p_value,
    test_statistic = test_statistic,
    # standard_error = standard_error,
    # numerator_test_statistic = numerator_test_statistic,
    fold_data = fold_data
  ))
}

#' Mean-difference convergence analysis for high-dimensional data. 
#'
#' The method combines sparse classification (lasso or group-lasso) and principal component analysis
#' to construct interpretable projection directions for hypothesis testing.
#'
#' @param control A matrix or data frame representing the control group. Rows are samples and columns are features.
#' @param treatment1 A matrix or data frame representing the first treatment group.
#' @param treatment2 A matrix or data frame representing the second treatment group.
#' @param pca_method Character. Method for estimating principal components. Must be either \code{"dense_pca"} or \code{"sparse_pca"}.
#' @param classifier_method Character. Method for feature selection. Must be either \code{"lasso"} or \code{"group_lasso"}.
#' @param lambda_type Character. Type of regularization strength to use from cross-validation. Either \code{"lambda.min"} or \code{"lambda.1se"}.
#' @param n_folds Integer. Number of folds for cross-validation. Default is 10.
#' @param group Optional. A grouping vector for use with group-lasso, with length equal to the number of columns in \code{control}.
#' @param standardize_feature Logical. Whether to normalize features using pooled mean and standard deviation before analysis. Default is \code{TRUE}.
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_value}{For the null hypothesis that treatment 1 and treatment 2 do not exhibit similar differences from the control, a smaller test statistic indicates greater similarity between the treatments.}
#'   \item{test_statistic}{Aggregated test statistic across folds.}
#'   \item{fold_data}{List of per-fold results including statistics, scores, and projections.}
#' }
#'
#' @export

convergence_testing <- function(
    control, treatment1, treatment2,
    pca_method = c("dense_pca", "sparse_pca"),
    classifier_method = c("lasso", "group_lasso"),
    lambda_type = 'lambda.1se',
    n_folds = 10,
    group = NULL,
    standardize_feature = TRUE,
    verbose = TRUE
) {
  # ============================================
  # Data Preprocessing: Validation and Conversion
  # ============================================
  control <- validate_and_convert_data(control, "control")
  treatment1 <- validate_and_convert_data(treatment1, "treatment1")
  treatment2 <- validate_and_convert_data(treatment2, "treatment2")
  
  check_non_null_and_identical_colnames(list(control, treatment1, treatment2))
  
  pca_method <- match.arg(pca_method) # match.arg ensures pca_method is one of the allowed values, defaulting to "dense_pca"
  classifier_method <- match.arg(classifier_method)
  
  if (!is.null(group) && classifier_method == 'lasso'){
    message("the grouping vector is not NULL but the method is normal LASSO, set classifier_method as group_lasso in convergence_testing()")
    
  }
  
  if (!is.null(group) && (!is.vector(group) || length(group) != ncol(control))) {
    stop("Error: `group` must be NULL or a vector of the same length as the number of columns in `control`.")
  }
  
  if(standardize_feature){
    # Normalize and split
    normalized_list <- normalize_and_split(control, treatment1, treatment2)
    
    # Access results
    control <- normalized_list$df1
    treatment1 <- normalized_list$df2
    treatment2 <- normalized_list$df3
  }
  # ============================================
  # Split Datasets into Folds
  # ============================================
  split_indices <- lapply(list(control, treatment1, treatment2), function(data) {
    check_data_for_folds(data, n_folds)
    index_spliter(1:nrow(data), n_folds)
  })
  
  control_split_index <- split_indices[[1]]
  tr1_split_index <- split_indices[[2]]
  tr2_split_index <- split_indices[[3]]
  
  fold_data <- vector("list", n_folds)
  
  # ============================================
  # Process data for each fold
  # ============================================  
  for(i in 1:n_folds){
    fold_data[[i]] <- process_fold(i, control, treatment1, treatment2, 
                                   control_split_index, tr1_split_index, tr2_split_index,
                                   pca_method, classifier_method, lambda_type, group, verbose)
  }
  
  # ===================================================
  # Now combine the folds
  # ===================================================
  return(combine_folds(fold_data, n_folds, verbose))
  
}

#' Identify consistently selected features across folds
#'
#' Aggregates non-zero features from cross-validation folds and applies majority voting
#' to identify stable, active features. Optionally returns active feature groups if a
#' grouping structure is provided.
#'
#' @param test_result A list returned by \code{convergence_testing()}, containing fold-wise results in the \code{$fold_data} element.
#' @param voting_method Character. Voting rule used to select features. Currently only \code{"majority_voting"} is supported.
#' @param group Optional. A named vector assigning each feature to a group. Required for identifying active groups.
#' @param group_threshold Integer. Minimum number of selected features in a group for it to be considered active. Default is 1.
#'
#' @return If \code{group} is \code{NULL}, returns a character vector of active feature names.
#' If \code{group} is provided, returns a list with:
#' \describe{
#'   \item{active_features}{Character vector of consistently selected features.}
#'   \item{active_groups}{Character vector of active group labels.}
#' }
#'
#' @export


collect_active_features <- function(test_result, voting_method = c("majority_voting"), 
                                    group = NULL, group_threshold = 1) {
  fold_data <- test_result$fold_data
  voting_method <- match.arg(voting_method)
  n_folds <- length(fold_data)
  active_features_list <- vector("list", n_folds)
  
  # Collect non-zero features for each fold
  for (i in 1:n_folds) {
    if (!is.null(fold_data[[i]]$final_beta)) {
      beta <- fold_data[[i]]$final_beta
    } else {
      beta <- fold_data[[i]]$classifier_coef
    }
    
    non_zero_features <- names(beta[abs(beta) > 1e-10])
    active_features_list[[i]] <- non_zero_features
  }
  
  # Flatten and count
  all_active_features <- unlist(active_features_list)
  feature_counts <- table(all_active_features)
  
  # Apply majority voting
  if (voting_method == 'majority_voting') {
    active_features <- names(feature_counts[feature_counts > n_folds / 2])
  }
  
  # Group handling
  if (!is.null(group)) {
    if (is.null(names(group))) {
      if (!is.null(fold_data[[1]]$final_beta)) {
        names(group) <- names(fold_data[[1]]$final_beta)
      } else {
        names(group) <- names(fold_data[[1]]$classifier_coef)
      }
    }
    
    group_nonzero_counts <- table(group[active_features])
    active_groups <- names(group_nonzero_counts[group_nonzero_counts >= group_threshold])
    
    return(list(
      active_features = active_features,
      active_groups = active_groups
    ))
  }
  
  return(active_features)
}
