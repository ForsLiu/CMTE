results_dir <- "results"
files <- list.files(results_dir, pattern = "^coef_est_n\\d+_p\\d+_r.*_eps.*_fn\\d+_rep\\d+\\.RData$", full.names = TRUE)

rows <- list()
max_beta_len <- 0

for (file in files) {
  matches <- regexec("coef_est_n(\\d+)_p(\\d+)_r([[:alnum:]x]+)_eps([0-9.]+)_fn(\\d+)_rep(\\d+).RData", basename(file))
  parts <- regmatches(basename(file), matches)[[1]]
  
  if (length(parts) != 7) next
  
  n <- as.integer(parts[2])
  p <- as.integer(parts[3])
  r_str <- parts[4]
  eps <- as.numeric(parts[5])
  f_num <- as.integer(parts[6])
  rep <- as.integer(parts[7])
  Omega <- "\\(\\exp(3)\\)"
  
  env <- new.env()
  load(file, envir = env)
  
  if (!all(c("beta_accuracy", "B_accuracy", "time_CMTE", "time_TRR") %in% ls(env))) next
  
  beta_acc <- round(env$beta_accuracy, 4)
  B_acc <- round(env$B_accuracy, 4)
  max_beta_len <- max(max_beta_len, length(beta_acc))
  
  row <- c(
    paste0("\\(", gsub("x", " \\\\times ", r_str), "\\)"),
    paste0("\\(", p, "\\)"),
    paste0("\\(", eps, "\\)"),
    paste0("\\(", n, "\\)"),
    Omega,
    paste0("\\(", f_num, "\\)"),
    paste0("\\(", rep, "\\)"),
    paste0("\\(", beta_acc, "\\)"),
    paste0("\\(", B_acc, "\\)"),
    paste0("\\(", round(env$time_CMTE, 2), "\\)"),
    paste0("\\(", round(env$time_TRR, 2), "\\)")
  )
  
  rows[[length(rows) + 1]] <- row
}

# Construct LaTeX table header
beta_headers <- paste0("\\(\\mathcal{D}(\\beta_{", seq_len(max_beta_len), "}, \\hat{\\beta}_{", seq_len(max_beta_len), "})\\)")
headers <- c(
  "\\(\\text{Dimension}\\)",
  "\\(p\\)",
  "\\(\\epsilon\\)",
  "\\(n\\)",
  "\\(\\Omega\\)",
  "\\(f\\)",
  "\\(\\text{rep}\\)",
  beta_headers,
  "\\(\\mathcal{D}(B, \\hat{B})\\)",
  "\\(\\text{Time}_{\\text{CMTE}}\\)",
  "\\(\\text{Time}_{\\text{TRR}}\\)"
)
header_row <- paste(headers, collapse = " & ")
hline <- "\\hline"

# Format rows with proper padding
latex_rows <- lapply(rows, function(row) {
  pad <- length(headers) - length(row)
  if (pad > 0) row <- c(row[1:(length(row) - 1)], rep("\\(\\text{NA}\\)", pad), row[length(row)])
  paste(row, collapse = " & ")
})

# Print full LaTeX table
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Estimation Accuracy and Runtime}\n")
cat("\\label{tab:all_results}\n")
cat(paste0("\\begin{tabular}{", paste(rep("l", length(headers)), collapse = ""), "}\n"))
cat(hline, "\n")
cat(header_row, " \\\\\n")
cat(hline, "\n")
cat(paste0(unlist(latex_rows), " \\\\\n", collapse = ""))
cat(hline, "\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
