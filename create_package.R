library(devtools)
library(roxygen2)
# setwd("/Users/tianyuzhang/Documents/")
# create("HMC")
setwd("/Users/tianyuzhang/Documents/HMC/")
document()
setwd("..")
install("HMC")

library(HMC)
?index_spliter
