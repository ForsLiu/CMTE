setwd("C:/Users/yh95l/Desktop/CMTE")


source("./parameters.R")
source("./DataGen.R")
source("./cmte.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

set.seed(123)





results_log <- foreach(i = 1:nrow(param_grid), .packages = c("rTensor", "MASS", "TRES"), .combine = rbind) %dopar% {
  r_vec <- r_vec_list[[param_grid$r_index[i]]]
  p     <- param_grid$p[i]
  eps   <- param_grid$eps[i]
  n     <- param_grid$n[i]
  Omega <- param_grid$Omega[i]
  f_num <- param_grid$f_num[i]
  
  r_str <- paste(r_vec, collapse = "x")
  input_file <- sprintf("data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_rep)
  
  if (!file.exists(input_file)) {
    return(data.frame(n = n, f_num = f_num, message = "File missing"))
  }
  
  load(input_file)  # loads Y_list, X_list, B_list_all
  
  cmte_acc_list   <- numeric(n_rep)
  trr_acc_list    <- numeric(n_rep)
  tmddm_acc_list  <- numeric(n_rep)
  
  for (rep in 1:n_rep) {
    Y <- Y_list[[rep]]
    X <- X_list[[rep]]
    beta_list <- B_list_all[[rep]]$beta_list
    
    # CMTE
    M_xy <- TMDDM(X@data, Y)
    cmte_est <- CMTE(X@data, Y, M_xy)
    cmte_acc_list[rep] <- beta_acc(cmte_est, beta_list)
    
    # TRR
    TReg <- TRR.fit(X@data, Y, u = rep(1, length(r_vec)), method = "1D")
    trr_est <- TReg$Gamma
    trr_acc_list[rep] <- beta_acc(trr_est, beta_list)
    
    # TMDDM
    M_xy <- TMDDM(X@data, Y)
    tmddm_est <- lapply(M_xy, function(Mk) eigen(Mk)$vectors[, 1])
    tmddm_acc_list[rep] <- beta_acc(tmddm_est, beta_list)
  }
  
  # Save result
  output_file <- sprintf("results/coef_est_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                         n, p, r_str, eps, f_num, n_rep)
  save(cmte_acc_list, trr_acc_list, tmddm_acc_list, file = output_file)
  
  # Return log message
  msg <- sprintf("n=%d, f_num=%d, r=%s | CMTE=%.4f, TRR=%.4f, TMDDM=%.4f",
                 n, f_num, r_str,
                 mean(cmte_acc_list, na.rm = TRUE),
                 mean(trr_acc_list, na.rm = TRUE),
                 mean(tmddm_acc_list, na.rm = TRUE))
  data.frame(n = n, f_num = f_num, r = r_str, message = msg)
}

print(results_log$message)



#generate table
source("./table_output.R")

#r_str <- paste(r_vec, collapse = "x")
#filename <- sprintf("data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",n, p, r_str, eps, f_num, n_rep)
#
#load(filename)
#
#
## CMTE
#time_CMTE <- system.time({
#  M_xy <- TMDDM(X@data, Y)
#  cmte_est <- CMTE(X@data, Y, M_xy)
#})[3] 
#cmte_accuracy <- beta_acc(cmte_est, beta_list)
#
## Tensor response regression
#time_TRR <- system.time({
#  TReg <- TRR.fit(X@data, Y, u = c(1,1), method = "1D")
#  trr_est <- TReg$Gamma
#})[3]
#trr_accuracy <- beta_acc(trr_est, beta_list)
#
##TMDDM
#time_tmddm <- system.time({
#  M_xy <- TMDDM(X@data, Y)
#  tmddm_est <- lapply(M_xy, function(Mk) {
#    eigen(Mk)$vectors[, 1]
#  })
#})[3] 
#tmddm_accuracy <- beta_acc(tmddm_est, beta_list)
#
#
## Save results
#filename_all <- sprintf("results/coef_est_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData", n, p, r_str, eps, f_num, n_rep)
#save(cmte_est, trr_est, tmddm_est,
#     cmte_accuracy, trr_accuracy, tmddm_accuracy,
#     time_CMTE, time_TRR, time_tmddm,
#     file = filename_all)
#message("Saved to: ", filename_all)



