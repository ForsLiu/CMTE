library(rTensor)
library(pracma)
library(parallel)
library(doParallel)
library(foreach)

TMDDM <- function(X, Y) {
  n <- nrow(X)
  Y_dims <- dim(Y)
  M <- length(Y_dims) - 1
  M_list <- vector("list", M)
  
  for (k in 1:M) {
    Yk <- k_unfold(as.tensor(Y), m = k)@data
    Yk_samples <- matrix(Yk, nrow = Y_dims[k])
    d_k <- Y_dims[k]
    Yk_mat <- matrix(0, nrow = n, ncol = d_k)
    prod_other <- prod(Y_dims[-c(k, M + 1)])
    
    for (i in 1:n) {
      idx_start <- (i - 1) * prod_other + 1
      idx_end <- i * prod_other
      Yk_mat[i, ] <- rowMeans(Yk[, idx_start:idx_end])
    }
    
    Yk_bar <- colMeans(Yk_mat)
    M_k <- matrix(0, d_k, d_k)
    
    for (i in 1:n) {
      for (j in 1:n) {
        Yi <- Yk_mat[i, ] - Yk_bar
        Yj <- Yk_mat[j, ] - Yk_bar
        diff <- X[i, ] - X[j, ]
        M_k <- M_k - (Yi %*% t(Yj)) * sqrt(sum(diff^2)) / (n^2)
      }
    }
    
    M_list[[k]] <- M_k
  }
  
  return(M_list)
}


CMTE <- function(X, Y, M_xy, eps = 1e-6) {

  Y_dims <- dim(Y)
  n <- Y_dims[length(Y_dims)]
  M <- length(Y_dims) - 1
  
  beta_list <- vector("list", M)
  
  for (k in 1:M) {
    d_k <- Y_dims[k]
    prod_other <- prod(Y_dims[-c(k, M + 1)])
    
    Yk <- k_unfold(as.tensor(Y), m = k)@data
    Yk_mat <- matrix(0, nrow = n, ncol = d_k)
    for (i in 1:n) {
      idx <- ((i - 1) * prod_other + 1):(i * prod_other)
      Yk_mat[i, ] <- rowMeans(Yk[, idx])
    }
    
    Sigma_Y <- cov(Yk_mat)
    M_k <- M_xy[[k]]
    B_k <- Sigma_Y + M_k + eps * diag(d_k)
    A_k <- Sigma_Y
    
    eig_out <- eigen(solve(B_k, A_k), symmetric = FALSE)
    beta_k <- Re(eig_out$vectors[, 1])
    beta_k <- beta_k / sqrt(sum(beta_k^2))
    beta_list[[k]] <- beta_k
  }
  
  return(beta_list)
}


