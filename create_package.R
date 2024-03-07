library(devtools)
library(roxygen2)
# setwd("/Users/tianyuzhang/Documents/")
# create("HMC")
setwd("/Users/tianyuzhang/Documents/HMC/HMC_package/")
document()
setwd("..")
build("HMC_package")
install("HMC_package")

library(HMC)
?index_spliter
?anchored_lasso_testing

# 
# 
# foo <- function(x){
#   pkgEnv <- new.env(parent=emptyenv())  
#   pkgEnv$hyperparameter <- -1
#   
#   print(pkgEnv$hyperparameter)
#   foo1(1, pkgEnv)
#   foo1(1, pkgEnv)
# }
# 
# foo1 <- function(x, env){
#   if(env$hyperparameter == -1){
#     print('ok')
#     assign("hyperparameter", 10, envir = env)
#     print(env$hyperparameter)
#   }else{
#     print('no')
#     print(env$hyperparameter)
#   }
# }
# 
# foo(1)
