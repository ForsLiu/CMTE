library(rTensor)
library(pracma)
library(parallel)
library(doParallel)
library(foreach)
library("TRES")

TMDDM <- function(X, Y) {
  n <- ncol(X)
  Y_dims <- dim(Y)
  M <- length(Y_dims) - 1 #number of mode of response
  M_list <- vector("list", M) #TMDDM output
  
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
        diff <- X[, i] - X[, j] 
        M_k <- M_k - (Yi %*% t(Yj)) * sqrt(sum(diff^2)) / (n^2)
      }
    }
    
    M_list[[k]] <- M_k
  }
  
  return(M_list)
}


CMTE <- function(X, Y, M_xy, eps = 1e-6, method) {

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
    Sigma_Y <- Sigma_Y + diag(1e-6, nrow(Sigma_Y)) #computationally singular 
    
    if (method == "1D") {
      beta_list[[k]] <- OptM1D(Sigma_Y, M_xy[[k]], 1)
    }
    else if (method == "ECD") {
      beta_list[[k]] <- ECD(Sigma_Y, M_xy[[k]], 1)
    }
  }
  
  return(beta_list)
}


