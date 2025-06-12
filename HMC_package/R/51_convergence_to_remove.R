# Below is the link to Vivek's code
#https://colab.research.google.com/drive/10ADmaLN83mOHdUQ7hMSgsjUW9Fka4G--?usp=sharing#scrollTo=Zfs4BXQ8iBF9
library(glmnet)
library(caret)
library(pracma)
library(MASS)
library(grpreg)
library(PMA)
# ---------------------------
# Utility Functions
# ---------------------------

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

validate_and_convert_data <- function(data, name) {
  if (!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "should be a matrix or a data frame."))
  }
  
  return(as.matrix(data))
}

check_non_null_and_identical_colnames <- function(data_list) {
  # Check for null or empty column names and ensure all column names are identical
  colnames_data <- lapply(data_list, colnames)
  
  # Check for null or empty column names
  for (i in 1:length(colnames_data)) {
    if (any(is.null(colnames_data[[i]])) || any(colnames_data[[i]] == "")) {
      stop(paste("Dataset", i, "contains null or empty column names. Please give the input data proper column names."))
    }
  }
  
  # Check if all column names are identical
  if (!all(sapply(colnames_data, function(x) all(x == colnames_data[[1]])))) {
    stop("The column names across the datasets are not identical. Please make sure all datasets have the same column names.")
  }
}

normalize_and_split <- function(df1, df2, df3 = NULL) {
  # Combine dataframes
  if (is.null(df3)) {
    combined <- rbind(df1, df2)
  } else {
    combined <- rbind(df1, df2, df3)
  }
  
  # Compute pooled mean and SD
  pooled_mean <- colMeans(combined)
  pooled_sd <- apply(combined, 2, sd)
  
  # Center and scale
  normalized <- scale(combined, center = pooled_mean, scale = pooled_sd)
  
  # Split back
  df1_norm <- normalized[1:nrow(df1), , drop = FALSE]
  df2_norm <- normalized[(nrow(df1) + 1):(nrow(df1) + nrow(df2)), , drop = FALSE]
  
  if (is.null(df3)) {
    return(list(df1 = df1_norm, df2 = df2_norm))
  } else {
    df3_norm <- normalized[(nrow(df1) + nrow(df2) + 1):nrow(combined), , drop = FALSE]
    return(list(df1 = df1_norm, df2 = df2_norm, df3 = df3_norm))
  }
}

check_data_for_folds <- function(data, n_folds) {
  if (nrow(data) < n_folds) stop("Not enough rows to create folds.")
}

# ---------------------------
# Statistical and Modeling Helpers
# ---------------------------

fit_lasso <- function(control_train, treat_train,
                      lambda_type = c("lambda.min", "lambda.1se"),
                      classifier_method = c("lasso", "group_lasso"),
                      group = NULL) {
  
  lambda_type <- match.arg(lambda_type)
  classifier_method <- match.arg(classifier_method)
  
  X_train <- rbind(control_train, treat_train)  
  y_train <- c(rep(0, nrow(control_train)), rep(1, nrow(treat_train)))  
  
  if (classifier_method == "group_lasso") {
    if (is.null(group)) stop("Group vector must be provided for group lasso.")
    
    lasso_result <- cv.grpreg(X = X_train, 
                              y = y_train, 
                              group = group, 
                              penalty = "grLasso", 
                              family = 'binomial')
    
  } else {
    lasso_result <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)
  }
  # Step 3: Extract LASSO coefficients at the optimal lambda
  # Extract coefficients at the optimal lambda
  beta_est <- coef(lasso_result, s = lambda_type)[-1]
  names(beta_est) <- colnames(X_train)
  
  # Calculate threshold for small coefficients
  
  n_effect <- 2 * min(nrow(control_train), nrow(treat_train))
  max_beta_element <- max(abs(beta_est))
  
  # Vectorized operation to threshold small coefficients to zero
  threshold <- max_beta_element * n_effect^(-1/3)
  beta_est[abs(beta_est) < threshold] <- 0
  
  # if(length(beta_est) == 50) beta_est[6:50] <- 0
  
  return(beta_est)
  
}

estimate_leading_pc <- function(control, pca_method = c("dense_pca", "sparse_pca")) {
  # Center the data
  centered_control <- sweep(as.data.frame(control), 2, colMeans(control))
  sample_centered <- as.matrix(centered_control)
  feature_number <- ncol(sample_centered)
  
  # Match PCA method argument
  pca_method <- match.arg(pca_method)
  print(pca_method)
  # Handle edge case: 1-dimensional data
  if (feature_number == 1) {
    warning("There is only one dimension and PCA is requested.")
    pc <- matrix(1, ncol = 1)
    names(pc) <- colnames(control)
    return(pc / sqrt(sum(pc^2)))
  }
  
  # Force dense PCA for low-dimensional data
  if (feature_number <= 30) {
    message("Dimension too small, switching method to dense PCA.")
    pca_method <- "dense_pca"
  }
  
  # Run PCA
  if (pca_method == "dense_pca") {
    print('conducting dense PCA')
    pc <- array(irlba::irlba(sample_centered, nv = 1)$v)
  } else if (pca_method == "sparse_pca") {
    print('conducting sparse PCA')
    cv_result <- PMA::SPC.cv(sample_centered)
    pc <- PMA::SPC(
      sample_centered,
      K = 1,
      sumabsv = cv_result$bestsumabsv
    )$v
  }
  
  # Normalize and return
  names(pc) <- colnames(control)
  return(pc / sqrt(sum(pc^2)))
}

# -------------------------
# all the functions above are also used in the anchored lasso test
# --------------------------

process_candidate_genes <- function(beta_1_con, group, classifier_method, control_train) {
  # Ensure the group vector has names corresponding to gene identifiers
  if (is.null(names(group))) {
    names(group) <- colnames(control_train)  # Assuming genes are columns in control_train
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

# ---------------------------
# Fold Functions
# ---------------------------

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
    
    next
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
  v_hat_tilde_diamond <- (pc_2_con + a_n * beta_2_con)  # Adjustment
  v_hat_tilde_diamond <- v_hat_tilde_diamond / norm(v_hat_tilde_diamond, type = "2")  # Normalize
  
  # ===================================================
  # Compute test statistic, scores, and variance
  # ===================================================
  # Extract relevant gene data from control and treatment groups
  control_matrix <- as.matrix(control_test[, candidate_genes], ncol = length(candidate_genes))
  tr2_matrix <- as.matrix(tr2_test[, candidate_genes], ncol = length(candidate_genes))
  
  # Compute scores for control and treatment groups
  control_score <- control_matrix %*% v_hat_tilde_diamond
  tr2_score <- tr2_matrix %*% v_hat_tilde_diamond
  
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
    proj_direction = v_hat_tilde_diamond,
    first_beta = beta_1_con,
    final_beta = beta_2_con,
    second_pc = pc_2_con
  ))
}

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

# ---------------------------
# Main Function
# ---------------------------

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

# ---------------------------
# Visualization Functions
# ---------------------------

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
