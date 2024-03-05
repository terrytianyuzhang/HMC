library(devtools)
library(roxygen2)
# setwd("/Users/tianyuzhang/Documents/")
# create("HMC")
setwd("/Users/tianyuzhang/Documents/HMC/")
document()
setwd("..")
build("HMC")
install("HMC")

library(HMC)
?index_spliter
?anchored_lasso_testing

