source("./parameters_1.3.R")
source("./alg_CMTE/cmte.R")
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
    input_file <- sprintf("Data_1.3/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                          n, p, r_str, eps, f_num, n_rep)
    
    if (!file.exists(input_file)) {
      return(data.frame(n = n, f_num = f_num, r = r_str, message = "File missing"))
    }
    
    load(input_file)  # loads Y_list, X_list, B_list_all
    
    cmte_1d_acc_list  <- numeric(n_rep)
    cmte_ecd_acc_list <- numeric(n_rep)
    cmte_1d_time_list <- numeric(n_rep)
    cmte_ecd_time_list <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      Y <- Y_list[[rep]]
      X <- X_list[[rep]]
      beta_list <- B_list_all[[rep]]$beta_list
      M_list <- TMDDM(X@data, Y)
      
      # CMTE - 1D
      t1 <- system.time({
        cmte_1d_est <- CMTE(X@data, Y, M_list, eps = 1e-6, method = "1D")
      })["elapsed"]
      cmte_1d_acc_list[rep] <- beta_acc(cmte_1d_est, beta_list)
      cmte_1d_time_list[rep] <- t1
      
      # CMTE - ECD
      t2 <- system.time({
        cmte_ecd_est <- CMTE(X@data, Y, M_list, eps = 1e-6, method = "ECD")
      })["elapsed"]
      cmte_ecd_acc_list[rep] <- beta_acc(cmte_ecd_est, beta_list)
      cmte_ecd_time_list[rep] <- t2
    }
    
    # Save result
    output_file <- sprintf("results_1.3/coef_est_cmte_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                           n, p, r_str, eps, f_num, n_rep)
    save(
      cmte_1d_acc_list,
      cmte_ecd_acc_list,
      cmte_1d_time_list,
      cmte_ecd_time_list,
      file = output_file
    )
    
    # Return log message
    data.frame(
      n        = rep(n, n_rep),
      f_num    = rep(f_num, n_rep),
      r        = rep(r_str, n_rep),
      rep      = seq_len(n_rep),
      CMTE_1D  = cmte_1d_acc_list,
      CMTE_1D_time = cmte_1d_time_list,
      CMTE_ECD = cmte_ecd_acc_list,
      CMTE_ECD_time = cmte_ecd_time_list
    )
    
  }, error = function(e) {
    data.frame(n = n, f_num = f_num, r = r_str, rep = NA, 
               CMTE_1D = NA, CMTE_1D_time = NA,
               CMTE_ECD = NA, CMTE_ECD_time = NA,
               message = paste("ERROR:", e$message))
  })
}

numeric_cols <- sapply(results_log, is.numeric)

numeric_cols["rep"]   <- FALSE
numeric_cols["n"]     <- FALSE
numeric_cols["f_num"] <- FALSE

group_cols <- setdiff(names(results_log), c(names(results_log)[numeric_cols], "rep"))

mean_results <- aggregate(
  results_log[, numeric_cols],
  by = results_log[, group_cols, drop = FALSE],
  FUN = mean
)
rep_counts <- aggregate(rep ~ ., data = results_log[, c(group_cols, "rep")], FUN = length)
names(rep_counts)[names(rep_counts) == "rep"] <- "total_rep"

mean_results <- merge(mean_results, rep_counts, by = group_cols)
mean_results$rep <- NULL

print(mean_results)

try(stopCluster(cl), silent = TRUE)  # stop old cluster if it exists
closeAllConnections()                # force-close socket and file connections
gc()                                 # trigger garbage collection