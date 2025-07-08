source("./parameters.R")
source("./Evaluation.R")
source("./alg_CMTE/cmte.R")  # for TMDDM and beta_acc

library(foreach)
library(doParallel)
library("TRES")

# Remove old TMDDM result files
unlink(list.files("results", pattern = "^coef_est_tmddm", full.names = TRUE), force = TRUE)

# Set up parallel backend
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
    
    tmddm_acc_list  <- vector("list", n_rep)
    tmddm_time_list <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      input_file <- sprintf("Data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                            n, p, r_str, eps, f_num, n_dir, n_rep,
                            n, p, r_str, eps, f_num, n_dir, rep)
      
      if (!file.exists(input_file)) {
        message(sprintf("Missing file: %s", input_file))
        tmddm_acc_list[[rep]] <- NA
        tmddm_time_list[rep] <- NA
        next
      }
      
      load(input_file)  # loads Y, X, B_list
      beta_list <- B_list$beta_list
      
      t1 <- system.time({
        M_list <- TMDDM(X@data, Y)
        tmddm_est <- lapply(M_list, function(Mk) eigen(Mk)$vectors[, 1])
      })["elapsed"]
      
      tmddm_acc_list[[rep]] <- beta_acc(tmddm_est, beta_list)
      tmddm_time_list[rep] <- t1
    }
    
    # Save result
    output_file <- sprintf("results/coef_est_tmddm_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                           n, p, r_str, eps, f_num, n_dir, n_rep)
    save(
      tmddm_acc_list,
      tmddm_time_list,
      file = output_file
    )
    
    # Return result log
    data.frame(
      n            = rep(n, n_rep),
      f_num        = rep(f_num, n_rep),
      r            = rep(r_str, n_rep),
      n_dir        = rep(n_dir, n_rep),
      rep          = seq_len(n_rep),
      TMDDM        = I(tmddm_acc_list),
      TMDDM_time   = tmddm_time_list,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    message(sprintf("ERROR: %s | n=%d p=%d r=%s fn=%d dir=%d", e$message, n, p, r_str, f_num, n_dir))
    data.frame(
      n            = rep(n, n_rep),
      f_num        = rep(f_num, n_rep),
      r            = rep(r_str, n_rep),
      n_dir        = rep(n_dir, n_rep),
      rep          = seq_len(n_rep),
      TMDDM        = I(as.list(rep(NA, n_rep))),
      TMDDM_time   = rep(NA, n_rep),
      stringsAsFactors = FALSE
    )
  })
}

# Clean up
try(stopCluster(cl), silent = TRUE)
closeAllConnections()
gc()
