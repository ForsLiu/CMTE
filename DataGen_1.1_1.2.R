library(MASS)
library(rTensor)

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

DGen_Omega0 <- function(r) {
  Or <- qr.Q(qr(matrix(rnorm((r - 1)^2), nrow = r - 1)))
  diag_vals <- exp(seq(2, -2, length.out = r - 1))
  Or %*% diag(diag_vals) %*% t(Or)
}

DGen_Sigma <- function(beta_list, Omega_c) {
  m <- length(beta_list)
  Sigma_list <- vector("list", m)
  
  for (i in seq_len(m)) {
    beta <- beta_list[[i]]
    r <- nrow(beta)
    d <- ncol(beta)
    d0 <- r - d                             
    
    Q_full <- qr.Q(qr(cbind(beta, matrix(runif(r * d0), nrow = r))))
    beta0 <- Q_full[, (d + 1):r, drop = FALSE]
    
    Or <- qr.Q(qr(matrix(rnorm(d0 * d0), nrow = d0)))
    diag_vals <- exp(seq(2, -2, length.out = d0))
    Omega0 <- Or %*% diag(diag_vals) %*% t(Or)
    
    Sigma_list[[i]] <- beta %*% Omega_c %*% t(beta) + beta0 %*% Omega0 %*% t(beta0)
  }
  
  Reduce(kronecker, rev(Sigma_list))
}



DGen_E <- function(n, Sigma) {
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * ncol(Sigma)), nrow = n)
  Z %*% t(L)
}


DGen_Y <- function(B, X, eps, Sigma, function_num) {
  m <- length(dim(B)) - 1
  r_vec <- dim(B)[1:m]
  n <- dim(X@data)[2]
  
  X_mat <- t(X@data)
  
  X_tilde <- switch(
    as.character(function_num),
    "1" = X_mat,
    "2" = exp(2 * abs(X_mat)),
    "3" = sin(X_mat),
    stop("Invalid function_num")
  )
  
  Bmat <- t(k_unfold(B, m = m + 1)@data)
  E <- DGen_E(n, Sigma)
  
  Y_vec <- X_tilde %*% t(Bmat) + eps * E
  array(t(Y_vec), dim = c(r_vec, n))
}


rp_grid <- unique(param_grid[c("r_index", "p")])

#generate parameter once


for (i in 1:nrow(param_grid)) {

  r_vec <- r_vec_list[[param_grid$r_index[i]]]
  p     <- param_grid$p[i]
  eps   <- param_grid$eps[i]
  n     <- param_grid$n[i]
  Omega_c <- param_grid$Omega_c[i]
  f_num <- param_grid$f_num[i]
  
  B_list    <- DGen_B(r_vec, p)
  B         <- B_list$B
  beta_list <- B_list$beta_list
  Sigma     <- DGen_Sigma(beta_list, exp(Omega_c))
  
  # Lists to collect replications
  Y_list <- vector("list", n_rep)
  X_list <- vector("list", n_rep)
  B_list_all <- vector("list", n_rep)

  start_time <- Sys.time()
  
  for (rep in 1:n_rep) {
    
    if (rep %% ceiling(n_rep / 5) == 0 || rep == 1 || rep == n_rep) {
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      message(sprintf("Progress: rep %d / %d | elapsed = %.1f sec | n = %d, p = %d, r = %s, eps = %.2f, f_num = %d",
                      rep, n_rep, as.numeric(elapsed),
                      n, p, paste(r_vec, collapse = "x"), eps, f_num))
    }
    
    X         <- DGen_X(n, p)
    Y         <- DGen_Y(B, X, eps, Sigma, function_num = f_num)
    
    Y_list[[rep]]     <- Y
    X_list[[rep]]     <- X
    B_list_all[[rep]] <- B_list
  }
  
  r_str <- paste(r_vec, collapse = "x")
  filename <- sprintf("Data_1.1_1.2/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData", n, p, r_str, eps, f_num, n_rep)
  save(Y_list, X_list, B_list_all, file = filename)
  message("Saved to: ", filename)
}
