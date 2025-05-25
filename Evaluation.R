beta_acc <- function(beta_est, beta_true) {
  M <- length(beta_true)
  D_vec <- numeric(M)
  
  for (i in 1:M) {
    B_true <- as.matrix(beta_true[[i]])
    B_est  <- as.matrix(beta_est[[i]])
    
    B_true <- qr.Q(qr(B_true))
    B_est  <- qr.Q(qr(B_est))
    
    P_true <- B_true %*% t(B_true)
    P_est  <- B_est %*% t(B_est)
    
    # Frobenius norm
    D_vec[i] <- norm(P_true - P_est, type = "F")
  }
  
  return(D_vec)
}


B_acc <- function(B_est, B_true) {
  B_est <- as.array(B_est)
  B_true <- as.array(B_true)
  
  diff_norm <- sqrt(sum((B_est - B_true)^2))
  true_norm <- sqrt(sum(B_true^2))
  
  rel_error <- diff_norm / true_norm
  return(rel_error)
}
