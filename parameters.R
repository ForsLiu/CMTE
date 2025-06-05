#exp1.1

#r_vec <- c(10,10)
#p     <- 3
#eps   <- .1
#n     <- 50
#Omega <- 3
#f_num <- 1
#n_rep <- 50


# parameter grids
r_vec_list <- list(c(10,10)) 
p_list     <- c(3)
eps_list   <- c(0.1)
n_list     <- c(10, 50, 100)
Omega_list <- c(3)
f_num_list <- c(1, 2)
n_rep      <- 2

#combinations
param_grid <- expand.grid(
  r1 = sapply(r_vec_list, `[`, 1),
  r2 = sapply(r_vec_list, `[`, 2),
  p = p_list,
  eps = eps_list,
  n = n_list,
  Omega = Omega_list,
  f_num = f_num_list
)