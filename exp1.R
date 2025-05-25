setwd("C:/Users/yh95l/Desktop/CMTE")

#change parameter and run all below
source("./parameters.R")
source("./DataGen.R")
source("./cmte.R")
source("./Evaluation.R")
library("TRES")

set.seed(123)

r_str <- paste(r_vec, collapse = "x")
filename <- sprintf("data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData",n, p, r_str, eps, f_num, n_rep)

load(filename)


# CMTE
time_CMTE <- system.time({
  M_xy <- TMDDM(X@data, Y)
  beta_est <- CMTE(X@data, Y, M_xy)
})[3]  # elapsed time
beta_accuracy <- beta_acc(beta_est, beta_list)

# Tensor regression
time_TRR <- system.time({
  TReg <- TRR.fit(X@data, Y)
  B_est <- TReg$coefficients@data
})[3]
B_accuracy <- B_acc(B_est, B@data)

# Save results
filename_all <- sprintf("results/coef_est_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData", n, p, r_str, eps, f_num, n_rep)
save(beta_est, B_est, beta_accuracy, B_accuracy, time_CMTE, time_TRR, file = filename_all)
message("Saved to: ", filename_all)

#generate table
source("./table_output.R")
