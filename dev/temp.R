dmatrixt <- function(q, m, M, nu, D, log=TRUE){
  S <- nrow(m)
  p <- ncol(m)
  
  uM <- chol(M)
  uD <- chol(D)
  Qbig <- q - m
  Qbig <- backsolve(uM, Qbig, transpose=TRUE)
  Qbig <- crossprod(Qbig)
  Qbig <- backsolve(uD, backsolve(uD, Qbig, transpose=TRUE))
  Qbig <- Qbig + diag(p)
  
  logp <- -p*S/2 * log(pi) + lmvgamma((nu + p)/2, S) - lmvgamma(nu/2, S) -
    S/2 * logdet(M) - p/2 * logdet(D) - (nu + p)/2 * logdet(Qbig)
  
  if (!log) return(exp(logp))
  return(logp)
}

lppd_IW_1t_dan <- function(Yt, Ft, st, St, Vt, nt, Dt) {
  lppd <- dmatrixt(Yt, Ft %*% st, Ft %*% St %*% t(Ft) + Vt,
                   nt, Dt, log=TRUE)
  return(lppd)
}



Yt = Y[[t]]; Ft = F_ls[[t]]; Vt = V_ls; st = para_ffbs_mniw$bs$st[,,t]; St = para_ffbs_mniw$bs$St[,,t];
nt = para_ffbs_mniw$ff[[t]]$nt; Dt = para_ffbs_mniw$ff[[t]]$Dt

dmatrixt(Yt, Ft %*% st, Ft %*% St %*% t(Ft) + Vt,
                 nt, Dt, log=TRUE)
mniw::dMT(X = Yt, Lambda = Ft %*% st, SigmaR = Ft %*% St %*% t(Ft) + Vt, 
                  nu = nt, SigmaC = Dt, log = TRUE)
MixMatrix::dmatrixt(x = Yt, mean = Ft %*% st, L = Ft %*% St %*% t(Ft) + Vt,
                            df = nt, R = Dt, log = TRUE)
q = Yt; m = Ft %*% st; M = Ft %*% St %*% t(Ft) + Vt; 
nu = nt; d = Dt; R = R; log = TRUE
dmatrixtscS <- function(q, m, M, nu, d, R, log=TRUE) {
  S <- nrow(q)
  p <- ncol(q)
  
  uM <- chol(M)
  uR <- chol(R)
  Qbig <- q - m
  Qbig <- backsolve(uM, Qbig, transpose=TRUE)
  Qbig <- backsolve(uR, t(Qbig), transpose=TRUE)
  
  logp <- -(p * S/2 + nu) * log(1 + sum(Qbig * Qbig)/(2 * d)) + CholWishart::lmvgamma(p*S/2 + nu, 1) - CholWishart::lmvgamma(nu, 1) - 
    p*S/2 * log(2 * pi * d) - S/2 * logdet(M) - p/2 * logdet(R)
  
  if (!log) return(exp(logp))
  return(logp)
}