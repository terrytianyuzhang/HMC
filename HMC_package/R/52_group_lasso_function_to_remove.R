fit_grlasso <- function(control, treatment, clustering){
  expression <- rbind(control, treatment)
  group_label <- c(rep(0, nrow(control)), rep(1, nrow(treatment)))
  
  cv_fit <- cv.grpreg(X = expression, 
                      y = group_label, 
                      group = clustering$gene_name, 
                      penalty = "grLasso", 
                      family = 'binomial')
  selected_beta <- coef(cv_fit, s = "lambda.min")
  plot(selected_beta)
  
  return(selected_beta)
}

group_truncation_num <- 5
control_subsetting_sample_size <- 1000