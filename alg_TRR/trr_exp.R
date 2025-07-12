#source("./parameters.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

unlink(list.files("results", pattern = "^coef_est_trr", full.names = TRUE), force = TRUE)

n_cores <- max(1, min(parallel::detectCores() - 2, nrow(param_grid)))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

set.seed(123)

results_log <- foreach(i = 1:nrow(param_grid), .packages = c("rTensor", "MASS", "TRES"), .combine = rbind) %dopar% {
  
  tryCatch({
    r_vec   <- r_vec_list[[param_grid$r_index[i]]]
    p       <- param_grid$p[i]
    eps     <- param_grid$eps[i]
    n       <- param_grid$n[i]
    Omega   <- param_grid$Omega[i]
    f_num   <- param_grid$f_num[i]
    n_dir   <- param_grid$n_dir[i]
    
    r_str <- paste(r_vec, collapse = "x")
    
    trr_ecd_acc_list  <- vector("list", n_rep)
    trr_ecd_time_list <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      input_file <- sprintf("Data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                            n, p, r_str, eps, f_num, n_dir, n_rep,
                            n, p, r_str, eps, f_num, n_dir, rep)
      
      if (!file.exists(input_file)) {
        message(sprintf("Missing file: %s", input_file))
        trr_ecd_acc_list[[rep]] <- NA
        trr_ecd_time_list[rep] <- NA
        next
      }
      
      load(input_file)  # loads Y, X, B_list
      beta_list <- B_list$beta_list
      
      # TRR - ECD only
      t2 <- system.time({
        TReg_ECD <- TRR.fit(X@data, Y, u = rep(1, length(r_vec)), method = "ECD")
      })["elapsed"]
      trr_ecd_acc_list[[rep]] <- beta_acc(TReg_ECD$Gamma, beta_list)
      trr_ecd_time_list[rep] <- t2
    }
    
    # Save result
    output_file <- sprintf("results/coef_est_trr_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                           n, p, r_str, eps, f_num, n_dir, n_rep)
    save(
      trr_ecd_acc_list,
      trr_ecd_time_list,
      file = output_file
    )
    
    # Return result log
    data.frame(
      n             = rep(n, n_rep),
      f_num         = rep(f_num, n_rep),
      r             = rep(r_str, n_rep),
      n_dir         = rep(n_dir, n_rep),
      rep           = seq_len(n_rep),
      TRR_ECD       = I(trr_ecd_acc_list),
      TRR_ECD_time  = trr_ecd_time_list,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    message(sprintf("ERROR: %s | n=%d p=%d r=%s fn=%d dir=%d", e$message, n, p, r_str, f_num, n_dir))
    data.frame(
      n             = rep(n, n_rep),
      f_num         = rep(f_num, n_rep),
      r             = rep(r_str, n_rep),
      n_dir         = rep(n_dir, n_rep),
      rep           = seq_len(n_rep),
      TRR_ECD       = I(as.list(rep(NA, n_rep))),
      TRR_ECD_time  = rep(NA, n_rep),
      stringsAsFactors = FALSE
    )
  })
}

try(stopCluster(cl), silent = TRUE)
closeAllConnections()
gc()
