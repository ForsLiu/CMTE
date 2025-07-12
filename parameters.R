# parameter grids
r_vec_list <- list(c(10,10),c(50,50),c(100,100))
#r_vec_list <- list(c(10,10),c(50,50))
p_list     <- c(3)
eps_list   <- c(0.1)
n_list     <- c(50,100,200,500)
Omega_list <- c(3)
f_num_list <- c(1,2,3)
n_rep      <- 2
n_dir_list <- c(1,2,3,4,5)
#n_dir_list <- c(1,2,3)

# combinations
param_grid <- expand.grid(
  r_index = seq_along(r_vec_list),
  p = p_list,
  eps = eps_list,
  n = n_list,
  Omega_c = Omega_list,
  f_num = f_num_list,
  n_dir = n_dir_list
)