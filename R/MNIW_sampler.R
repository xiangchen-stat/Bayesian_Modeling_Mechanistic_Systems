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
MNIW_sampler <- function(nsam, X, H, v, S, C, Vb){
  n <- nrow(X)
  q <- ncol(S)
  p <- ncol(X)
  
  sigma <- rinvwishart(nsam, nu = v, Omega = S, epsilon = 0, checkSymmetry = F)
  zero <- matrix(0, n, q)
  B <- array(dim = c(p, q, nsam))
  E <- array(dim = c(n, q, nsam))
  # Y <- array(dim = c(n, q, nsam))
  for (i in 1:nsam) {
    Bpre <- rmatrixnormal(1, M = C, U = Vb, V = sigma[,,i], checkSymmetry = F, keep = TRUE)
    B[,,i] <- matrix(Bpre, nrow = p, ncol = q)
    Epre <- rmatrixnormal(1, zero, H, sigma[,,i])
    E[,,i] <- matrix(Epre, nrow = n, ncol = q)
    # Y[,,i] <- X %*% B[,,i] + E[,,i]
  }
  
  # out <- list(sigma = sigma, B = B, Y = Y)
  out <- list(sigma = sigma, B = B)
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
