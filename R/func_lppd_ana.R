## lppd analytical functions
dMTig <- function(Y, m, M, nu, d, R, log=TRUE) {
  p <- nrow(Y)
  S <- ncol(Y)
  Rinv <- solve(R)
  Minv <- solve(M)
  q <- Y - m
  Q <- Rinv %*% t(q) %*% Minv %*% q
  temp <- 1 + 1/(2*d) * LaplacesDemon::tr(Q)
  out <- lgamma(p*S/2 + nu) - lgamma(nu) - p*S/2 * log(2 * pi * d) - 
    S/2 * LaplacesDemon::logdet(M) - p/2 * LaplacesDemon::logdet(R) - 
    (p*S/2 + nu) * log(temp)
  
  if (log == TRUE) {
    return(out)
  } else{
    return(exp(out))
  }
}


lppd_IW_1t <- function(Yt, Ft, Vt, st, St, nt, Dt) {
  # note nu = nt - S + 1; hyper T is not matrix t
  S <- dim(Yt)[2]
  lppd <- mniw::dMT(X = Yt, Lambda = Ft %*% st, SigmaR = Ft %*% St %*% t(Ft) + Vt,
                    nu = nt - S + 1, SigmaC = Dt, log = TRUE)
  # lppd <- mniw::dMT(X = Yt, Lambda = Ft %*% st, SigmaR = Ft %*% St %*% t(Ft) + Vt,
  #                   nu = nt, SigmaC = Dt, log = TRUE)
  return(lppd)
}

# lppd_IG_1t_dan <- function(Yt, Ft, Vt, st, St, nt, Dt, R) {
#   lppd <- dmatrixtscS(q = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt, 
#                       nu = nt, d = Dt, R = R, log = TRUE)
#   return(lppd)
# }


lppd_IG_1t <- function(Yt, Ft, Vt, st, St, nt, Dt, R) {
  # note nu = nt - S + 1; hyper T is not matrix t
  S <- dim(Yt)[2]
  # lppd <- dMTig(Y = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt, 
  #               nu = nt - S + 1, d = Dt, R = R, log = TRUE)
  lppd <- dMTig(Y = Yt, m = Ft %*% st, M = Ft %*% St %*% t(Ft) + Vt, 
                nu = nt, d = Dt, R = R, log = TRUE)
  return(lppd)
}

lppd_id_1t <- function (Yt, Ft, Vt, st, St) {
  lppd <- mniw::dMNorm(X = Yt, Lambda = Ft %*% st,
                       SigmaR = make_pds(Ft %*% St %*% t(Ft) + Vt), SigmaC = diag(ncol(Yt)), log=TRUE)
  return(lppd)
}