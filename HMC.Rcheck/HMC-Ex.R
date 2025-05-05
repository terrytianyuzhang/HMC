pkgname <- "HMC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "HMC-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('HMC')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("anchored_lasso_testing")
### * anchored_lasso_testing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: anchored_lasso_testing
### Title: Anchored test for two-sample mean comparison.
### Aliases: anchored_lasso_testing

### ** Examples

sample_size_1 <- sample_size_2 <- 300
true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)

sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
                               mu = true_mean_1,
                               Sigma = diag(1, 100)))
 sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
                               mu = true_mean_2,
                               Sigma = diag(1, 100)))
 result <- anchored_lasso_testing(sample_1, sample_2)
 result$test_statistics
 ##the test statistic. It should follow normal(0,1) when there is no difference between the groups.
 summarize_feature_name(result) 
 #summarize which features contribute to discriminant vectors (i.e. logistic lasso)
 extract_pc(result) # extract the estimated discriminant coefficients
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("anchored_lasso_testing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_data_for_folds")
### * check_data_for_folds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_data_for_folds
### Title: Check that data has enough rows for cross-validation folds
### Aliases: check_data_for_folds

### ** Examples

check_data_for_folds(matrix(1:20, nrow = 5), n_folds = 5)
## Not run: 
##D check_data_for_folds(matrix(1:4, nrow = 2), n_folds = 5)  # This will throw an error
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_data_for_folds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_non_null_and_identical_colnames")
### * check_non_null_and_identical_colnames

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_non_null_and_identical_colnames
### Title: Check non-null and consistent column names across datasets
### Aliases: check_non_null_and_identical_colnames

### ** Examples

d1 <- data.frame(a = 1:2, b = 3:4)
d2 <- data.frame(a = 5:6, b = 7:8)
check_non_null_and_identical_colnames(list(d1, d2))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_non_null_and_identical_colnames", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("debiased_pc_testing")
### * debiased_pc_testing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: debiased_pc_testing
### Title: Debiased one-step test for two-sample mean comparison. A small
###   p-value tells us not only there is difference in the mean vectors,
###   but can also indicates which principle component the difference
###   aligns with.
### Aliases: debiased_pc_testing

### ** Examples

sample_size_1 <- sample_size_2 <- 300

true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)
pc1 <- c(rep(1, 10), rep(0, 90))
pc1 <- pc1/norm(pc1, type = '2')

simulation_covariance <- 10 * pc1 %*% t(pc1)
simulation_covariance <- simulation_covariance + diag(1, 100)

sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
                               mu = true_mean_1,
                               Sigma = simulation_covariance))
 sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
                               mu = true_mean_2,
                               Sigma = simulation_covariance))
 result <- debiased_pc_testing(sample_1, sample_2)
 result$test_statistics
 ##these are test statistics. Each one of them corresponds to one PC.
 summarize_pc_name(result, latent_fator_index = 1) #shows which features contribute to PC1
 extract_pc(result) # extract the estimated leading PCs.
 
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("debiased_pc_testing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_leading_pc")
### * estimate_leading_pc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_leading_pc
### Title: Estimate the leading principal component
### Aliases: estimate_leading_pc

### ** Examples

## Not run: 
##D X <- matrix(rnorm(100), nrow = 20)
##D estimate_leading_pc(X, pca_method = "dense_pca")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_leading_pc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_lasso")
### * fit_lasso

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_lasso
### Title: Fit a (group) Lasso logistic regression classifier
### Aliases: fit_lasso

### ** Examples

## Not run: 
##D X1 <- matrix(rnorm(100), nrow = 10)
##D X2 <- matrix(rnorm(100), nrow = 10)
##D fit_lasso(X1, X2, classifier_method = "lasso")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_lasso", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("index_spliter")
### * index_spliter

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: index_spliter
### Title: Split indices into folds
### Aliases: index_spliter

### ** Examples

index_spliter(1:10, n_folds = 3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("index_spliter", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mean_comparison_anchor")
### * mean_comparison_anchor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mean_comparison_anchor
### Title: High-dimensional two-sample mean comparison with anchored
###   projection
### Aliases: mean_comparison_anchor

### ** Examples

## Not run: 
##D X <- matrix(rnorm(200 * 100), nrow = 100)
##D Y <- matrix(rnorm(200 * 100), nrow = 100)
##D result <- mean_comparison_anchor(X, Y, pca_method = "dense_pca", classifier_method = "lasso")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mean_comparison_anchor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalize_and_split")
### * normalize_and_split

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalize_and_split
### Title: Normalize and split two datasets using pooled mean and standard
###   deviation
### Aliases: normalize_and_split

### ** Examples

set.seed(123)
df1 <- matrix(rnorm(20), nrow = 5)
df2 <- matrix(rnorm(20), nrow = 5)
normalize_and_split(df1, df2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalize_and_split", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simple_pc_testing")
### * simple_pc_testing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simple_pc_testing
### Title: Simple plug-in test for two-sample mean comparison.
### Aliases: simple_pc_testing

### ** Examples

sample_size_1 <- sample_size_2 <- 300
true_mean_1 <- matrix(c(rep(1, 10), rep(0, 90)), ncol = 1)
true_mean_2 <- matrix(c(rep(1.5, 10), rep(0, 90)), ncol = 1)
pc1 <- c(rep(1, 10), rep(0, 90))
pc1 <- pc1/norm(pc1, type = '2')

simulation_covariance <- 10 * pc1 %*% t(pc1)
simulation_covariance <- simulation_covariance + diag(1, 100)

sample_1 <- data.frame(MASS::mvrnorm(sample_size_1,
                                     mu = true_mean_1,
                                     Sigma = simulation_covariance))
sample_2 <- data.frame(MASS::mvrnorm(sample_size_2,
                                     mu = true_mean_2,
                                     Sigma = simulation_covariance))
result <- simple_pc_testing(sample_1, sample_2)
result$test_statistics
##these are test statistics. Each one of them corresponds to one PC.
summarize_pc_name(result, latent_fator_index = 1) #shows which features contribute to PC1
extract_pc(result) # extract the estimated leading PCs.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simple_pc_testing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("validate_and_convert_data")
### * validate_and_convert_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: validate_and_convert_data
### Title: Validate and convert input data
### Aliases: validate_and_convert_data

### ** Examples

validate_and_convert_data(data.frame(x = 1:3, y = 4:6), "example_data")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("validate_and_convert_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
