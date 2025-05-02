#' Split indices into folds
#'
#' Randomly splits a given vector of indices into approximately equal-sized folds.
#'
#' @param array A vector of indices (e.g., `1:n`) to be split into folds.
#' @param n_folds Integer. Number of folds. Default is 5.
#'
#' @return A list of length `n_folds`, each containing a subset of the shuffled indices.
#'
#' @examples
#' index_spliter(1:10, n_folds = 3)
#'
#' @export

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

#' Validate and convert input data
#'
#' Checks whether the input is a matrix or data frame, and converts it to a matrix if valid.
#'
#' @param data A matrix or data frame.
#' @param name A string used in error messages to identify the variable name.
#'
#' @return A numeric matrix.
#'
#' @examples
#' validate_and_convert_data(data.frame(x = 1:3, y = 4:6), "example_data")
#'
#' @export


validate_and_convert_data <- function(data, name) {
  if (!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "should be a matrix or a data frame."))
  }
  
  return(as.matrix(data))
}

#' Check non-null and consistent column names across datasets
#'
#' Ensures all input datasets have non-null, non-empty, and identical column names.
#'
#' @param data_list A list of matrices or data frames to be checked.
#'
#' @return NULL (called for side-effect). Throws an error if validation fails.
#'
#' @examples
#' d1 <- data.frame(a = 1:2, b = 3:4)
#' d2 <- data.frame(a = 5:6, b = 7:8)
#' check_non_null_and_identical_colnames(list(d1, d2))
#'
#' @export


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

#' Normalize and split two datasets using pooled mean and standard deviation
#'
#' Combines two datasets, normalizes features using pooled mean and standard deviation,
#' and returns the normalized datasets separately.
#'
#' @param df1 A data frame or matrix. Typically group 1.
#' @param df2 A data frame or matrix. Typically group 2.
#'
#' @return A list with elements:
#' \describe{
#'   \item{df1}{Normalized version of `df1`.}
#'   \item{df2}{Normalized version of `df2`.}
#' }
#'
#' @examples
#' set.seed(123)
#' df1 <- matrix(rnorm(20), nrow = 5)
#' df2 <- matrix(rnorm(20), nrow = 5)
#' normalize_and_split(df1, df2)
#'
#' @export

normalize_and_split <- function(df1, df2) {
  # Combine
  combined <- rbind(df1, df2)
  
  # Compute pooled mean and SD
  pooled_mean <- colMeans(combined)
  pooled_sd <- apply(combined, 2, sd)
  
  # Center and scale
  normalized <- scale(combined, center = pooled_mean, scale = pooled_sd)
  
  # Split back
  df1_norm <- normalized[1:nrow(df1), , drop = FALSE]
  df2_norm <- normalized[(nrow(df1) + 1):nrow(combined), , drop = FALSE]
  
  return(list(df1 = df1_norm, df2 = df2_norm))
}

#' Check that data has enough rows for cross-validation folds
#'
#' Validates that the input data has at least as many rows as the number of desired folds.
#'
#' @param data A data frame or matrix.
#' @param n_folds Integer. The number of folds to check for.
#'
#' @return NULL (called for its side effect). Throws an error if the number of rows is too small.
#'
#' @examples
#' check_data_for_folds(matrix(1:20, nrow = 5), n_folds = 5)
#' \dontrun{
#' check_data_for_folds(matrix(1:4, nrow = 2), n_folds = 5)  # This will throw an error
#' }
#'
#' @export


check_data_for_folds <- function(data, n_folds) {
  if (nrow(data) < n_folds) stop("Not enough rows to create folds.")
}

#' Fit a (group) Lasso logistic regression classifier
#'
#' Performs Lasso or group Lasso logistic regression to distinguish between two groups of samples.
#'
#' @param control_train A matrix or data frame for the control group. Rows are samples, columns are features.
#' @param treat_train A matrix or data frame for the treatment group. Rows are samples, columns are features.
#' @param lambda_type Character. Type of lambda to use from cross-validation. Options are `"lambda.min"` (default) and `"lambda.1se"`.
#' @param classifier_method Character. Choice of classifier. `"lasso"` (default) or `"group_lasso"`.
#' @param group Optional grouping vector for `group_lasso`, same length as the number of columns in the input data.
#'
#' @return A numeric vector of estimated regression coefficients (excluding intercept), thresholded for small values.
#'
#' @details The function fits a logistic regression using either `glmnet` for Lasso or `grpreg` for group Lasso.
#' Coefficients are soft-thresholded by the maximum coefficient times `n^(-1/3)` where `n` is the effective sample size.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom grpreg cv.grpreg
#'
#' @examples
#' \dontrun{
#' X1 <- matrix(rnorm(100), nrow = 10)
#' X2 <- matrix(rnorm(100), nrow = 10)
#' fit_lasso(X1, X2, classifier_method = "lasso")
#' }
#'
#' @export


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

#' Estimate the leading principal component
#'
#' Estimates the leading principal component of the input matrix using dense or sparse PCA.
#'
#' @param control A matrix or data frame. Each row is a sample, and each column is a feature.
#' @param pca_method Character. PCA method to use. Options are `"dense_pca"` (default) or `"sparse_pca"`.
#'
#' @return A normalized numeric vector representing the leading principal component direction.
#'
#' @details For low-dimensional settings (≤ 30 features), the method automatically switches to dense PCA.
#' For sparse PCA, the function uses the `PMA::SPC.cv` cross-validation method.
#'
#' @importFrom irlba irlba
#' @importFrom PMA SPC.cv SPC
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 20)
#' estimate_leading_pc(X, pca_method = "dense_pca")
#' }
#'
#' @export

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
