setwd("C:/Users/yh95l/Desktop/CMTE")

file.remove(list.files("Data", full.names = TRUE))
file.remove(list.files("results", full.names = TRUE))


source("./parameters.R")
source("./DataGen_1.3.R")
source("./cmte.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

set.seed(123)


results_log <- foreach(i = 1:nrow(param_grid), .packages = c("rTensor", "MASS", "TRES"), .combine = rbind) %dopar% {
  
  tryCatch({
  
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
      return(data.frame(n = n, f_num = f_num, r = r_str, message = "File missing"))
    }
    
    load(input_file)  # loads Y_list, X_list, B_list_all
    
    cmte_1d_acc_list  <- numeric(n_rep)
    cmte_ecd_acc_list <- numeric(n_rep)
    trr_1d_acc_list   <- numeric(n_rep)
    trr_ecd_acc_list  <- numeric(n_rep)
    tmddm_acc_list    <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      Y <- Y_list[[rep]]
      X <- X_list[[rep]]
      beta_list <- B_list_all[[rep]]$beta_list
      M_list <- TMDDM(X@data, Y)
      
      # CMTE - 1D
      cmte_1d_est <- CMTE(X@data, Y, M_list, eps = 1e-6, method = "1D")
      cmte_1d_acc_list[rep] <- beta_acc(cmte_1d_est, beta_list)
      
      # CMTE - ECD
      cmte_ecd_est <- CMTE(X@data, Y, M_list, eps = 1e-6, method = "ECD")
      cmte_ecd_acc_list[rep] <- beta_acc(cmte_ecd_est, beta_list)
      
      # TRR - 1D
      TReg_1D <- TRR.fit(X@data, Y, u = rep(1, length(r_vec)), method = "1D")
      trr_1d_est <- TReg_1D$Gamma
      trr_1d_acc_list[rep] <- beta_acc(trr_1d_est, beta_list)
      
      # TRR - ECD
      TReg_ECD <- TRR.fit(X@data, Y, u = rep(1, length(r_vec)), method = "ECD")
      trr_ecd_est <- TReg_ECD$Gamma
      trr_ecd_acc_list[rep] <- beta_acc(trr_ecd_est, beta_list)
      
      # TMDDM
      M_list <- TMDDM(X@data, Y)
      tmddm_est <- lapply(M_list, function(Mk) eigen(Mk)$vectors[, 1])
      tmddm_acc_list[rep] <- beta_acc(tmddm_est, beta_list)
    }
    
    # Save result
    output_file <- sprintf("results/coef_est_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                           n, p, r_str, eps, f_num, n_rep)
    save(
      cmte_1d_acc_list,
      cmte_ecd_acc_list,
      trr_1d_acc_list,
      trr_ecd_acc_list,
      tmddm_acc_list,
      file = output_file
    )
    
    # Return log message
    data.frame(
      n        = rep(n, n_rep),
      f_num    = rep(f_num, n_rep),
      r        = rep(r_str, n_rep),
      rep      = seq_len(n_rep),
      CMTE_1D  = cmte_1d_acc_list,
      CMTE_ECD = cmte_ecd_acc_list,
      TRR_1D   = trr_1d_acc_list,
      TRR_ECD  = trr_ecd_acc_list,
      TMDDM    = tmddm_acc_list
    )
  
  }, error = function(e) {
    data.frame(n = n, f_num = f_num, r = r_str, rep = NA, CMTE_1D = NA, CMTE_ECD = NA,
               TRR_1D = NA, TRR_ECD = NA, TMDDM = NA, message = paste("ERROR:", e$message))
  })
}

print(results_log)



#generate table
source("./table_output.R")
