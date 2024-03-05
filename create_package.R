library(devtools)
library(roxygen2)
# setwd("/Users/tianyuzhang/Documents/")
# create("HMC")
setwd("/Users/tianyuzhang/Documents/HMC/HMC_package/")
document()
setwd("..")
build("HMC_package")
install("HMC")

library(HMC)
?index_spliter
?anchored_lasso_testing

