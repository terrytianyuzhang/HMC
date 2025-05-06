# Clean up environment
rm(list = ls())

# Load development tools
library(devtools)

# Define package path
pkg_path <- "/Users/tianyuzhang/Documents/HMC/HMC_package"

# Step 1: Document package (generate NAMESPACE and man/)
document(pkg_path)

# Step 2: Build source tar.gz (optional if you're just installing locally)
build(pkg_path)

# Step 3: Install the package from source
# 
check(pkg = pkg_path, args = "--as-cran")
install(pkg_path)
# Step 4: Load the package
library(HMC)

# Step 5: Check documentation for key functions
help("index_spliter", package = "HMC")
help("check_data_for_folds", package = "HMC")
help("estimate_leading_pc", package = "HMC")
help("collect_active_features_proj", package = "HMC")
?HMC

