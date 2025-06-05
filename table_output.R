# Load results
rows <- list()

for (i in 1:nrow(param_grid)) {
  r_vec <- c(param_grid$r1[i], param_grid$r2[i])
  p     <- param_grid$p[i]
  eps   <- param_grid$eps[i]
  n     <- param_grid$n[i]
  Omega <- param_grid$Omega[i]
  f_num <- param_grid$f_num[i]
  
  r_str <- paste(r_vec, collapse = "x")
  Omega_latex <- "\\(\\exp(3)\\)"
  
  f_expr <- switch(as.character(f_num),
                   "1" = "\\(X\\)",
                   "2" = "\\(\\exp(2|X|)\\)",
                   "f(X) = \\text{Unknown}")
  
  file <- sprintf("results/coef_est_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",
                  n, p, r_str, eps, f_num, n_rep)
  
  if (!file.exists(file)) {
    warning("Missing file: ", file)
    next
  }
  
  env <- new.env()
  load(file, envir = env)
  
  acc_lists <- list(
    "CMTE"  = env$cmte_acc_list,
    "TRR"   = env$trr_acc_list,
    "TMDDM" = env$tmddm_acc_list
  )
  
  methods <- names(acc_lists)
  n_methods <- length(methods)
  
  # Precompute accuracy map
  method_acc_map <- lapply(acc_lists, function(acc) {
    acc <- round(acc, 4)
    while (length(acc) < 2) acc <- c(acc, NA)
    acc
  })
  
  # Determine bold minimum for each beta_i
  acc_matrix <- do.call(rbind, method_acc_map)
  min_indices <- apply(acc_matrix, 2, function(col) {
    if (all(is.na(col))) return(rep(FALSE, length(col)))
    col == min(col, na.rm = TRUE)
  })
  
  for (j in seq_along(methods)) {
    method <- methods[j]
    acc_vec <- method_acc_map[[method]]
    
    acc_fmt <- sapply(seq_along(acc_vec), function(k) {
      val <- acc_vec[k]
      if (is.na(val)) return("\\(\\text{NA}\\)")
      formatted <- formatC(val, format = "f", digits = 4)
      if (min_indices[j, k]) paste0("\\(\\textbf{", formatted, "}\\)") else paste0("\\(", formatted, "\\)")
    })
    
    if (j == 1) {
      row <- c(
        paste0("\\multirow{", n_methods, "}{*}{\\(", gsub("x", " \\\\times ", r_str), "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", p, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", eps, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", n, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{", Omega_latex, "}"),
        paste0("\\multirow{", n_methods, "}{*}{\\(", n_rep, "\\)}"),
        paste0("\\multirow{", n_methods, "}{*}{", f_expr, "}"),
        paste0("\\(", method, "\\)"),
        acc_fmt[1], acc_fmt[2]
      )
    } else {
      row <- c(
        rep("", 7),
        paste0("\\(", method, "\\)"),
        acc_fmt[1], acc_fmt[2]
      )
    }
    
    rows[[length(rows) + 1]] <- row
  }
}

# Table headers (no time column)
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
  "\\(\\mathcal{D}(\\beta_2, \\hat{\\beta}_2)\\)"
)

# Output LaTeX table
header_row <- paste(headers, collapse = " & ")
latex_rows <- lapply(rows, function(row) paste(row, collapse = " & "))

cat("\\begin{table}[ht]\n")
cat("\\centering\n")

cat("\\label{tab:beta_results}\n")
cat("\\begin{tabular}{", paste(rep("l", length(headers)), collapse = ""), "}\n")
cat("\\hline\n")
cat(header_row, " \\\\\n")
cat("\\hline\n")

for (i in seq_along(latex_rows)) {
  cat(latex_rows[[i]], " \\\\\n")
  if (i %% 3 == 0) cat("\\hline\n")  # Add line after every 3 rows (i.e., each file)
}

cat("\\end{tabular}\n")
cat("\\caption{Estimation Accuracy for \\(\\beta_i\\)'s from CMTE, TRR, and TMDDM}\n")
cat("\\end{table}\n")