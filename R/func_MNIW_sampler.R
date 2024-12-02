require(matrixsampling)
require(invgamma)
check_pds <- function(C, eps = 10^(-5), per = 0.05) {
  # print(C)
  norm_f <- norm(C - t(C), type="F")
  condition <- ((norm_f / mean(C)) < per) || (norm_f < eps)
  if (condition) {
    C <- (C + t(C))/2
    if (!is.positive.definite(C)) {
      C <- C + 3 * 10^(-4) * diag(dim(C)[1])
      warning("Matrix is not positive definite.")
      # print(tail(eigen(C)$values))
    }
    return(C)
  } else {
    print(norm_f)
    stop("Matrix symmetry is outside of tolerance: ",eps)
  }
}


make_pds <- function(C, eps = 10^(-4)) {
  C <- (C + t(C))/2
  if (!is.positive.definite(C)) {
    C <- C + eps * diag(dim(C)[1])
    warning("Matrix is not positive definite.")
  }
  return(C)
}
# if(F){
#   nsam = nsam; X = Ft; H = Vt; v = para_nT$nt;
#   S = para_nT$Dt; C = para_nT$mt; Vb = para_nT$Mt
# }
# TODO: vectorize sampling?
#' Sample from matrix-normal inverse wishart distribution
#' Y = XB + E, E~MN(O, H, Sigma)
#' Sigma~IW(v, S)
#' B | Sigma~MN(C, Vb, Sigma)
#' @param nsam Number of samples
#' @param X Covariate matrix
#' @param H Row covariance matrix of Y
#' @param v Degrees of freedom of IW
#' @param S Scale matrix of IW
#' @param C Mean matrix of B | Sigma
#' @param Vb Row covariance matrix of  B | Sigma
#' @returns List of samples of B and Sigma
#' @export
MNIW_sampler <- function(nsam, X, v, S, C, Vb){
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(S)
  
  # adjust symmetric
  # H <- check_pds(H)
  # S <- check_pds(S)
  # Vb <- check_pds(Vb)
  
  sigma <- matrixsampling::rinvwishart(nsam, nu = v, Omega = S, epsilon = 0, checkSymmetry = F)
  B <- array(dim = c(p, q, nsam))
  
  # zero <- matrix(0, n, q)
  # E <- array(dim = c(n, q, nsam))
  # Y <- array(dim = c(n, q, nsam))
  
  for (i in 1:nsam) {
    B[,,i] <- rmn_cpp(m = C, U = Vb, V = sigma[,,i])
    # B[,,i] <- rmatrixnormal(1, M = C, U = Vb, V = sigma[,,i], checkSymmetry = F, keep = F)
    # E[,,i] <- rmatrixnormal(1, zero, H, sigma[,,i], checkSymmetry = F, keep = F)
    # Y[,,i] <- X %*% B[,,i] + E[,,i]
  }
  
  # out <- list(sigma = sigma, B = B, Y = Y)
  out <- list(sigma = sigma, B = B)
  return(out)
}

#' @param nsam Number of samples
#' @param X Covariate matrix
#' @param C Mean matrix of B | Sigma
#' @param Vb Row covariance matrix of  B | Sigma
#' @param v shape of inverse gamma
#' @param S rate of inverse gamma
#' @param R covariance kernel
#' @returns List of samples of B and Sigma
#' @export
MNIG_sampler <- function(nsam, X, v, S, C, Vb, R){
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(C)
  
  # adjust symmetric
  sigma2 <- invgamma::rinvgamma(n = nsam, shape = v, rate = S)
  B <- array(dim = c(p, q, nsam))
  
  for (i in 1:nsam) {
    B[,,i] <- rmn_cpp(m = C, U = Vb, V = sigma2[i] * R)
    # B[,,i] <- rmatrixnormal(1, M = C, U = Vb, V = sigma2[i] * R, checkSymmetry = F, keep = F)
  }
  
  out <- list(sigma2 = sigma2, B = B)
  return(out)
}


# MNIW_exact <- function(X, Y, B, sigma, H, H_til, J, seed = 123){
#   set.seed(seed)
#   Hinv = solve(H)
#   mean = X %*% B + t(J) %*% Hinv %*% (Y - X %*% B) 
#   Hnew = H_til - t(J) %*% Hinv %*% J
#   
#   out <- rmatrixnormal(1, mean, Hnew, sigma)
#   return(out)
# }

