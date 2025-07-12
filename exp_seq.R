setwd("C:/Users/yh95l/Desktop/CMTE")

library(foreach)
library(doParallel)
library("TRES")

unlink(list.files("Data", full.names = TRUE), recursive = TRUE, force = TRUE)
file.remove(list.files("results", full.names = TRUE))

# Define parameter lists
r_vec_list <- list(c(10,10), c(50,50), c(100,100))
r_vec_list <- list(c(10,10), c(50,50))
p_list     <- c(3)
eps_list   <- c(0.1)
n_list     <- c(50, 100, 200, 500)
n_list     <- c(50, 100)
Omega_list <- c(3)
f_num_list <- c(1, 2, 3)
f_num_list <- c(2, 3)
n_rep      <- 2

for (n_dir in 1:3) {
  cat(sprintf("=== [%s] Running for n_dir = %d ===\n", Sys.time(), n_dir))
  
  param_grid <- expand.grid(
    r_index = seq_along(r_vec_list),
    p       = p_list,
    eps     = eps_list,
    n       = n_list,
    Omega_c = Omega_list,
    f_num   = f_num_list,
    n_dir   = n_dir
  )
  
  source("./DataGen.R")
  source("./Evaluation.R")
  
  cat(sprintf("=== [%s] Running CMTE ===\n", Sys.time()))
  source("./alg_CMTE/cmte_exp.R")
  
  cat(sprintf("=== [%s] Running TRR ===\n", Sys.time()))
  source("./alg_TRR/trr_exp.R")
  
  cat(sprintf("=== [%s] Running TMDDM ===\n", Sys.time()))
  source("./alg_TMDDM/tmddm_exp.R")
  
  # Cleanup Data
  unlink(list.files("Data", full.names = TRUE), recursive = TRUE, force = TRUE)
}


#generate table
source("./table_output.R")