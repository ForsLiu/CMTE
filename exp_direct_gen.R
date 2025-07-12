setwd("C:/Users/yh95l/Desktop/CMTE")

library(foreach)
library(doParallel)
library("TRES")

file.remove(list.files("results", full.names = TRUE))

source("./parameters.R")
source("./data_gen_functions.R")
source("./Evaluation.R")

source("./alg_CMTE/cmte.R")

n_cores <- max(1, min(parallel::detectCores() - 2, nrow(param_grid)))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

foreach(i = 1:nrow(param_grid), .packages = c("rTensor", "MASS", "TRES")) %dopar% {
  tryCatch({
    r_vec   <- r_vec_list[[param_grid$r_index[i]]]
    p       <- param_grid$p[i]
    eps     <- param_grid$eps[i]
    n       <- param_grid$n[i]
    Omega   <- param_grid$Omega[i]
    f_num   <- param_grid$f_num[i]
    n_dir   <- param_grid$n_dir[i]
    
    r_str <- paste(r_vec, collapse = "x")
    start_time <- Sys.time()
    
    cmte_ecd_acc_list  <- vector("list", n_rep)
    cmte_ecd_time_list <- numeric(n_rep)
    
    tmddm_acc_list     <- vector("list", n_rep)
    tmddm_time_list    <- numeric(n_rep)
    
    trr_ecd_acc_list   <- vector("list", n_rep)
    trr_ecd_time_list  <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      ## --- Data Generation ---
      B_list    <- DGen_B(r_vec, p)
      B         <- B_list$B
      beta_list <- B_list$beta_list
      Sigma     <- DGen_Sigma(beta_list, exp(Omega), n_dir)
      L         <- chol(Sigma)
      X         <- DGen_X(n, p)
      Y         <- DGen_Y(B, X, eps, L, function_num = f_num)
      
      ## --- TMDDM ---
      t_tmddm <- system.time({
        M_list <- TMDDM(X@data, Y)
        tmddm_est <- lapply(M_list, function(Mk) eigen(Mk)$vectors[, 1])
      })["elapsed"]
      tmddm_acc_list[[rep]] <- beta_acc(tmddm_est, beta_list)
      tmddm_time_list[rep]  <- t_tmddm
      
      ## --- CMTE ---
      t_cmte <- system.time({
        cmte_ecd_est <- CMTE(X@data, Y, M_list, eps = 1e-6, method = "ECD")
      })["elapsed"]
      cmte_ecd_acc_list[[rep]] <- beta_acc(cmte_ecd_est, beta_list)
      cmte_ecd_time_list[rep]  <- t_cmte
      
      ## --- TRR ---
      t_trr <- system.time({
        TReg_ECD <- TRR.fit(X@data, Y, u = rep(1, length(r_vec)), method = "ECD")
      })["elapsed"]
      trr_ecd_acc_list[[rep]] <- beta_acc(TReg_ECD$Gamma, beta_list)
      trr_ecd_time_list[rep]  <- t_trr
      
    }
    
    ## --- Save Results ---
    save(cmte_ecd_acc_list, cmte_ecd_time_list,
         file = sprintf("results/coef_est_cmte_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep))
    save(tmddm_acc_list, tmddm_time_list,
         file = sprintf("results/coef_est_tmddm_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep))
    save(trr_ecd_acc_list, trr_ecd_time_list,
         file = sprintf("results/coef_est_trr_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep))
    
  }, error = function(e) {
    message(sprintf("ERROR: %s | n=%d p=%d r=%s fn=%d dir=%d",
                    e$message, param_grid$n[i], param_grid$p[i],
                    paste(r_vec_list[[param_grid$r_index[i]]], collapse = "x"),
                    param_grid$f_num[i], param_grid$n_dir[i]))
  })
}


try(stopCluster(cl), silent = TRUE)  
closeAllConnections()                
gc()      


#generate table
source("./table_output.R")

