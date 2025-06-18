rm(list = ls())
library(HMC)
library(MASS)
library(ggplot2)

set.seed(123)
simulate_null_data <- function(n_samples=500, n_features=100) {
  mu_con <- rep(0, n_features)
  mu_tr1 <- c(rep(5, 5), rep(0, n_features - 5))
  mu_tr2 <- c(rep(0, 5),  rep(0, n_features - 10), rep(5, 5))
  
  ratio <- 3
  control <- mvrnorm(n_samples*ratio, mu_con, diag(c(1:5, rep(1, n_features - 5))))
  treatment1 <- mvrnorm(n_samples, mu_tr1, diag(n_features))
  treatment2 <- mvrnorm(n_samples, mu_tr2, diag(n_features))
  
  rownames(control) <- paste("Sample", 1:(n_samples*ratio), sep = "_")
  rownames(treatment1) <- paste("Sample", 1:n_samples, sep = "_")
  rownames(treatment2) <- paste("Sample", 1:n_samples, sep = "_")
  
  colnames(control) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment1) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment2) <- paste("Feature", 1:n_features, sep = "_")
  
  return(list(control = control, treatment1 = treatment1, treatment2 = treatment2))
}



# Simulation settings
n_simulations <- 100
n_samples <- 50
n_features <- 200
n_folds <- 5

test_statistics <- numeric(n_simulations)

# Run simulations
for (i in 1:n_simulations) {
  set.seed(i)
  data <- simulate_null_data(n_samples = n_samples, n_features = n_features)
  result <- convergence_testing(data$control, data$treatment1, data$treatment2,
                                lambda_type = 'lambda.min', n_folds = n_folds)
  test_statistics[i] <- result$test_statistic
}

# Histogram of test statistics
ggplot(data.frame(test_statistics), aes(x = test_statistics)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Test Statistics under Null", x = "Test Statistic", y = "Count") +
  theme_minimal()

# Q-Q Plot using ggplot2
qq_df <- data.frame(
  sample_quantiles = sort(test_statistics),
  theoretical_quantiles = qnorm(ppoints(n_simulations))
)

ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
  geom_point(color = "steelblue", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Test Statistics", x = "Theoretical Quantiles (Normal)", y = "Sample Quantiles") +
  theme_minimal()
