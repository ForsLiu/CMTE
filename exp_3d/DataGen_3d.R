library(MASS)
library(rTensor)
library(doParallel)

set.seed(123)
DGen_B <- function(r_vec, p) {
  m <- length(r_vec)
  
  Theta <- array(runif(p), dim = c(rep(1, m), p))
  Theta <- as.tensor(Theta)
  
  beta_list <- lapply(r_vec, function(r) {
    beta <- runif(r, -1, 1)
    qr.Q(qr(matrix(beta, ncol = 1)))
  })
  
  I_p <- diag(p)
  
  B <- Theta
  for (i in seq_len(m)) {
    B <- ttl(B, list(beta_list[[i]]), ms = i)
  }
  B <- ttl(B, list(I_p), ms = m + 1)
  
  return(list(Theta = Theta,B = B, beta_list = beta_list))
}

DGen_X <- function(n, p) {
  p <- prod(p)
  X <- matrix(rnorm(n * p), nrow = n)
  X_tensor <- array(t(X), dim = c(p, n))
  return(as.tensor(X_tensor))
}

DGen_Omega <- function(m,c) {
  Omega <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      Omega[i, j] <- 0.5^abs(i - j) * c
    }
  }
  return(Omega)
}

DGen_Omega0 <- function(r,n_dir) {
  if (r - n_dir <= 1) {
    return(matrix(exp(2), nrow = 1, ncol = 1))
  }
  
  Or <- qr.Q(qr(matrix(rnorm((r - n_dir)^2), nrow = r - n_dir)))
  diag_vals <- exp(seq(2, 2 - 0.5 * (r - n_dir - 1), by = -0.5))
  Or %*% diag(diag_vals) %*% t(Or)
}

DGen_Sigma <- function(beta_list, Omega_c, n_dir) {
  m <- length(beta_list)
  Sigma_list <- vector("list", m)
  
  for (i in seq_len(m)) {
    beta <- beta_list[[i]]  # r x d
    r <- nrow(beta)
    d <- ncol(beta)
    
    stopifnot(d + n_dir <= r)
    
    # alpha: r x n_dir
    alpha <- matrix(runif(r * n_dir, -1, 1), nrow = r)
    
    # Combine beta and alpha, then orthonormalize
    tilde_beta <- qr.Q(qr(cbind(beta, alpha)))  # r x (d + n_dir)
    
    # Extend to full orthonormal basis
    tilde_beta_full <- qr.Q(qr(cbind(tilde_beta, matrix(runif(r * (r - (d + n_dir))), nrow = r))))
    
    # Omega for informative subspace
    Omega <- DGen_Omega(n_dir, Omega_c)  # n_dir x n_dir
    
    # Omega0 for complement subspace
    Omega0 <- DGen_Omega0(r,n_dir)  # (r - n_dir) x (r - n_dir)
    
    # Select complement directions
    tilde_beta0 <- tilde_beta_full[, (d + n_dir):r, drop = FALSE]  # r x (r - d - n_dir + 1)
    
    # Sigma = informative + complement parts
    Sigma_i <- tilde_beta[, (d + 1):(d + n_dir), drop = FALSE] %*% Omega %*% 
               t(tilde_beta[, (d + 1):(d + n_dir), drop = FALSE]) +
               tilde_beta0 %*% Omega0 %*% t(tilde_beta0)
    
    Sigma_list[[i]] <- Sigma_i
  }

  Sigma <- Sigma_list[[1]]
  if (m > 1) {
    for (i in 2:m) {
      Sigma <- kronecker(Sigma_list[[i]], Sigma)
    }
  }

  Sigma <- Sigma + diag(1e-6, nrow(Sigma))
  
  return(Sigma)
}

DGen_E <- function(n, L) {
  Z <- matrix(rnorm(n * ncol(L)), nrow = n)
  Z %*% t(L)
}

DGen_Y <- function(B, X, eps, L, function_num) {
  m <- length(dim(B)) - 1
  dims_B <- dim(B)
  r_vec <- dims_B[-(m + 1)]
  p <- dims_B[m + 1]
  n <- dim(X@data)[2]
  
  # Get X matrix of shape n x p
  X_mat <- t(X@data)
  
  # Apply transformation
  if (function_num == 1) {
    X_tilde <- X_mat
  } else if (function_num == 2) {
    X_tilde <- exp(2 * abs(X_mat))
  } else if (function_num == 3) {
    X_tilde <- X_mat ^ 2
  } else {
    stop("Invalid function_num.")
  }
  
  # Unfold B along last mode
  Bmat <- t(k_unfold(B, m = m + 1)@data)
  
  # Generate noise
  E <- DGen_E(n, L)
  
  # Compute Y
  Y_vec <- X_tilde %*% t(Bmat) + eps * E
  Y <- array(t(Y_vec), dim = c(r_vec, n))
  
  return(Y)
}

rp_grid <- unique(param_grid[c("r_index", "p")])


n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

for (i in 1:nrow(param_grid)) {
  
  r_vec    <- r_vec_list[[param_grid$r_index[i]]]
  p        <- param_grid$p[i]
  eps      <- param_grid$eps[i]
  n        <- param_grid$n[i]
  Omega_c  <- param_grid$Omega_c[i]
  f_num    <- param_grid$f_num[i]
  n_dir    <- param_grid$n_dir[i]
  
  r_str <- paste(r_vec, collapse = "x")
  dir_name <- sprintf("Data/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d", 
                      n, p, r_str, eps, f_num, n_dir, n_rep)
  dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
  
  message(sprintf("=== Start: n = %d, p = %d, r = %s, eps = %.2f, f_num = %d, n_dir = %d ===", 
                  n, p, r_str, eps, f_num, n_dir))
  
  t_start <- Sys.time()
  
  B_list    <- DGen_B(r_vec, p)
  B         <- B_list$B
  beta_list <- B_list$beta_list
  Sigma     <- DGen_Sigma(beta_list, exp(Omega_c), n_dir)
  L         <- chol(Sigma)
  
  foreach(rep = 1:n_rep, .packages = c("rTensor", "MASS")) %dopar% {
    set.seed(rep)
    
    X <- DGen_X(n, p)
    Y <- DGen_Y(B, X, eps, L, function_num = f_num)
    
    file_name <- sprintf("%s/SimData_n%d_p%d_r%s_eps%.2f_fn%d_dir%d_rep%d.RData", 
                         dir_name, n, p, r_str, eps, f_num, n_dir, rep)
    save(Y, X, B_list, file = file_name)
  }
  
  t_end <- Sys.time()
  message(sprintf("=== Completed: n = %d, p = %d, r = %s, eps = %.2f, f_num = %d, n_dir = %d | Total Time: %.2f sec ===\n",
                  n, p, r_str, eps, f_num, n_dir, as.numeric(difftime(t_end, t_start, units = "secs"))))
}

# Stop parallel backend
stopCluster(cl)