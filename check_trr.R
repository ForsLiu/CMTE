results <- list()

for (i in 1:nrow(param_grid)) {
  r_vec <- r_vec_list[[param_grid$r_index[i]]]
  p     <- param_grid$p[i]
  n     <- param_grid$n[i]
  eps   <- param_grid$eps[i]
  f_num <- 2
  n_dir <- param_grid$n_dir[i]
  r_str <- paste(r_vec, collapse = "x")
  
  if (f_num != 2) next  # only check exp(2|x|)
  
  for (rep in 1:n_rep) {
    file_path <- sprintf("Data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                         n, p, r_str, eps, f_num, n_dir, rep,
                         n, p, r_str, eps, f_num, n_dir, rep)
    
    if (!file.exists(file_path)) next
    load(file_path)  # loads X (already transformed)
    
    X_exp <- t(X@data)                 # n x p matrix of exp(2|x|)
    X_orig_abs <- log(X_exp) / 2       # |x| approximately
    X_corrs <- diag(cor(X_exp, X_orig_abs))
    mean_corr <- mean(X_corrs, na.rm = TRUE)
    results[[length(results)+1]] <- data.frame(
      n = n, r = r_str, rep = rep,
      mean_corr = mean_corr
    )
  }
}

# Combine results
df_corr <- do.call(rbind, results)
print(df_corr)
