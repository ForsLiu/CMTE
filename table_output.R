rows <- list()

for (i in 1:nrow(param_grid)) {
  r_vec <- r_vec_list[[param_grid$r_index[i]]]
  p     <- param_grid$p[i]
  eps   <- param_grid$eps[i]
  n     <- param_grid$n[i]
  Omega <- param_grid$Omega[i]
  f_num <- param_grid$f_num[i]
  n_dir <- param_grid$n_dir[i]
  
  r_str <- paste(r_vec, collapse = "x")
  r_latex <- paste0("\\(", paste(r_vec, collapse = " \\times "), "\\)")
  Omega_latex <- paste0("\\(", n_dir, " \\times ", n_dir, "\\)")
  
  f_expr <- switch(as.character(f_num),
                   "1" = "\\(X\\)",
                   "2" = "\\(\\exp(2|X|)\\)",
                   "3" = "\\(X^2\\)",
                   "f(X) = \\text{Unknown}")
  
  file_cmte  <- sprintf("results/coef_est_cmte_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep)
  file_trr   <- sprintf("results/coef_est_trr_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep)
  file_tmddm <- sprintf("results/coef_est_tmddm_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData",
                        n, p, r_str, eps, f_num, n_dir, n_rep)
  
  acc_lists <- list()
  time_lists <- list()
  
  if (file.exists(file_cmte)) {
    env <- new.env(); load(file_cmte, envir = env)
    acc_lists[["CMTE-ECD"]] <- env$cmte_ecd_acc_list
    time_lists[["CMTE-ECD"]] <- env$cmte_ecd_time_list
  }
  
  if (file.exists(file_trr)) {
    env <- new.env(); load(file_trr, envir = env)
    acc_lists[["TRR-ECD"]] <- env$trr_ecd_acc_list
    time_lists[["TRR-ECD"]] <- env$trr_ecd_time_list
  }
  
  if (file.exists(file_tmddm)) {
    env <- new.env(); load(file_tmddm, envir = env)
    acc_lists[["TMDDM"]] <- env$tmddm_acc_list
    time_lists[["TMDDM"]] <- env$tmddm_time_list
  }
  
  if (length(acc_lists) == 0) next
  
  method_acc_map <- lapply(acc_lists, function(acc) {
    acc <- round(unlist(acc), 4)
    while (length(acc) < 2) acc <- c(acc, NA)
    acc
  })
  
  method_time_map <- lapply(time_lists, function(t) {
    avg <- if (length(t) > 0) mean(t, na.rm = TRUE) else NA
    round(avg, 2)
  })
  
  acc_matrix <- do.call(rbind, method_acc_map)
  min_indices <- apply(acc_matrix, 2, function(col) {
    if (all(is.na(col))) return(rep(FALSE, length(col)))
    col == min(col, na.rm = TRUE)
  })
  
  methods <- names(acc_lists)
  n_methods <- length(methods)
  
  for (j in seq_along(methods)) {
    method <- methods[j]
    acc_vec <- method_acc_map[[method]]
    
    avg_time <- method_time_map[[method]]
    min_time <- suppressWarnings(min(unlist(method_time_map), na.rm = TRUE))
    
    if (!is.null(avg_time) && !is.na(avg_time)) {
      formatted_time <- formatC(avg_time, format = "f", digits = 2)
      if (!is.na(min_time) && avg_time == min_time) {
        time_fmt <- paste0("\\(\\textbf{", formatted_time, "}\\)")
      } else {
        time_fmt <- paste0("\\(", formatted_time, "\\)")
      }
    } else {
      time_fmt <- "\\(\\text{NA}\\)"
    }
    
    acc_fmt <- sapply(seq_along(acc_vec), function(k) {
      val <- acc_vec[k]
      if (is.na(val)) return("\\(\\text{NA}\\)")
      formatted <- formatC(val, format = "f", digits = 4)
      if (min_indices[j, k]) paste0("\\(\\textbf{", formatted, "}\\)") else paste0("\\(", formatted, "\\)")
    })
    
    if (j == 1) {
      row <- c(
        paste0("\\multirow{", n_methods, "}{*}{", r_latex, "}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", p, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", eps, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", n, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{", Omega_latex, "}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", n_rep, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{", f_expr, "}"),
        paste0("\\(", method, "\\)"),
        acc_fmt[1], acc_fmt[2], time_fmt
      )
    } else {
      row <- c(rep("", 7), paste0("\\(", method, "\\)"), acc_fmt[1], acc_fmt[2], time_fmt)
    }
    
    rows[[length(rows) + 1]] <- row
  }
}

headers <- c(
  "\\(\\text{Dimension}\\)",
  "\\(p\\)",
  "\\(\\epsilon\\)",
  "\\(n\\)",
  "\\(\\Omega\\)",
  "\\(\\text{Replications}\\)",
  "\\(f\\)",
  "\\(\\text{Method}\\)",
  "\\(\\mathcal{D}(\\beta_1, \\hat{\\beta}_1)\\)",
  "\\(\\mathcal{D}(\\beta_2, \\hat{\\beta}_2)\\)",
  "\\(\\text{Avg Time (s)}\\)"
)

cat("\\begin{longtable}{", paste(rep("l", length(headers)), collapse = ""), "}\n", sep = "")
cat("\\hline\n")
cat(paste(headers, collapse = " & "), " \\\\\n")
cat("\\hline\n")
cat("\\endfirsthead\n")
cat("\\hline\n")
cat(paste(headers, collapse = " & "), " \\\\\n")
cat("\\hline\n")
cat("\\endhead\n")
cat("\\hline\n")
cat("\\endfoot\n")
cat("\\endlastfoot\n")

block_size <- 3
for (i in seq_along(rows)) {
  cat(paste(rows[[i]], collapse = " & "), " \\\\\n")
  if (i %% block_size == 0) {
    cat("\\hline\n")
  }
}

cat("\\caption{Estimation Accuracy and Runtime for \\(\\beta_i\\)'s from CMTE, TRR, and TMDDM.}\n")
cat("\\end{longtable}\n")