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

DGen_Omega <- function(m,c) {
  Omega <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      Omega[i, j] <- 0.5^abs(i - j) * c
    }
  }
  return(Omega)
}

DGen_Omega0 <- function(r) {
  Or <- qr.Q(qr(matrix(rnorm((r - 2)^2), nrow = r - 2)))
  diag_vals <- exp(seq(2, -1.5, length.out = r - 2))
  Or %*% diag(diag_vals) %*% t(Or)
}

DGen_Sigma <- function(beta_list, Omega_c) {
  m <- length(beta_list)
  Sigma_list <- vector("list", m)
  #alpha_list <- vector("list", m)
  
  for (i in seq_len(m)) {
    beta <- beta_list[[i]]
    r <- nrow(beta)
    d <- ncol(beta)
    d0 <- r - d
    
    alpha <- matrix(runif(r, -1, 1), nrow = r) #10x1
    tilde_beta <- qr.Q(qr(cbind(beta, alpha)))  # 10 x 2 orthogonalized, t(tilde_beta) %*% tilde_beta in I_2x2
    tilde_beta_full <- qr.Q(qr(cbind(tilde_beta, matrix(runif(r * d0), nrow = r))))
    
    Omega <- DGen_Omega(2,Omega_c) #2x2
    Omega0 <- DGen_Omega0(r)  #(ri-2)x(ri-2) = 8x8
    
    tilde_beta0 <- tilde_beta_full[, 3:r, drop = FALSE] #10x8
    
    Sigma_i <- tilde_beta %*% Omega %*% t(tilde_beta) + tilde_beta0 %*% Omega0 %*% t(tilde_beta0)
    #10x2 * 2x2 2x10 + 10x8 * 8x8 *8x10  
    Sigma_list[[i]] <- Sigma_i
    #alpha_list[[i]] <- alpha
  }
  
  Sigma <- Sigma_list[[1]]
  if (m > 1) {
    for (i in 2:m) {
      Sigma <- kronecker(Sigma_list[[i]], Sigma)
    }
  }
  
  return(Sigma)
}


DGen_E <- function(n, Sigma) {
  d <- ncol(Sigma)
  E <- mvrnorm(n, mu = rep(0, d), Sigma = Sigma)
  return(E)
}


DGen_Y <- function(B, X, eps, Sigma, function_num) {
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
    X_tilde <- sin(X_mat)
  } else {
    stop("Invalid function_num.")
  }
  
  # Unfold B along last mode
  Bmat <- t(k_unfold(B, m = m + 1)@data)
  
  # Generate noise
  E <- DGen_E(n, Sigma)
  
  # Compute Y
  Y_vec <- X_tilde %*% t(Bmat) + eps * E
  Y <- array(t(Y_vec), dim = c(r_vec, n))
  
  return(Y)
}


rp_grid <- unique(param_grid[c("r_index", "p")])



for (i in 1:nrow(param_grid)) {
  r_vec <- r_vec_list[[param_grid$r_index[i]]]
  p     <- param_grid$p[i]
  eps   <- param_grid$eps[i]
  n     <- param_grid$n[i]
  Omega_c <- param_grid$Omega_c[i]
  f_num <- param_grid$f_num[i]
  
  #generate parameter once
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
  filename <- sprintf("Data_1.3/SimData_n%d_p%d_r%s_eps%.2f_fn%d_rep%d.RData", n, p, r_str, eps, f_num, n_rep)
  save(Y_list, X_list, B_list_all, file = filename)
  message("Saved to: ", filename)
}