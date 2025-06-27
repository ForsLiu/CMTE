source("./parameters_1.3.R")
source("./Evaluation.R")

library(foreach)
library(doParallel)
library("TRES")

unlink(list.files("results_1.3", pattern = "^coef_est_cmte", full.names = TRUE), force = TRUE)

n_cores <- max(1, min(parallel::detectCores() - 2, nrow(param_grid)))
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
    
    tmddm_acc_list  <- numeric(n_rep)
    tmddm_time_list <- numeric(n_rep)
    
    for (rep in 1:n_rep) {
      input_file <- sprintf("Data_1.3/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                            n, p, r_str, eps, f_num, n_rep,
                            n, p, r_str, eps, f_num, rep)
      
      if (!file.exists(input_file)) {
        message(sprintf("Missing file: %s", input_file))
        next
      }
      
      load(input_file)  # loads Y, X, B_list
      beta_list <- B_list$beta_list
      
      t1 <- system.time({
        M_list <- TMDDM(X@data, Y)
        tmddm_est <- lapply(M_list, function(Mk) eigen(Mk)$vectors[, 1])
      })["elapsed"]
      
      tmddm_acc_list[rep] <- beta_acc(tmddm_est, beta_list)
      tmddm_time_list[rep] <- t1
    }
    
    # Save result
    output_file <- sprintf("results_1.3/coef_est_tmddm_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                           n, p, r_str, eps, f_num, n_rep)
    save(
      tmddm_acc_list,
      tmddm_time_list,
      file = output_file
    )
    
    # Return log message
    data.frame(
      n          = rep(n, n_rep),
      f_num      = rep(f_num, n_rep),
      r          = rep(r_str, n_rep),
      rep        = seq_len(n_rep),
      TMDDM      = tmddm_acc_list,
      TMDDM_time = tmddm_time_list
    )
    
  }, error = function(e) {
    data.frame(n = n, f_num = f_num, r = r_str, rep = NA,
               TMDDM = NA, TMDDM_time = NA,
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