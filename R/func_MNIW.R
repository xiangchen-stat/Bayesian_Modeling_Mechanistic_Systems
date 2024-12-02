MNIW_R_naiive <- function(X, Y, newX = none, newY = none, H, v, S, C, Vb){
  out <- list()
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  Vbinv <- solve(Vb)
  Hinv <- solve(H)
  Xt <- t(X)
  Yt <- t(Y)
  Ct <- t(C)
  Oinv <- Vbinv + Xt %*% Hinv %*% X
  O <- solve(Oinv)
  M <- Vbinv %*% C + Xt %*% Hinv %*% Y
  Mt <- t(M)
  Snew <- S + Ct %*% Vbinv %*% C + Yt %*% Hinv %*% Y - Mt %*% O %*% M
  vnew <- v + n
  Cnew <- O %*% M
  Vbnew <- O

  out[["vnew"]] <- vnew  
  out[["Snew"]] <- Snew
  out[["Cnew"]] <- Cnew
  out[["Vbnew"]] <- Vbnew
  return(out)
}

MNIW_R <- function(X, Y, newX = none, newY = none, H, v, S, C, Vb){
  out <- list()
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  
  Vbinv <- solve(Vb)
  Hinv <- solve(H)
  XtHinv <- crossprod(X, Hinv)
  Oinv <- Vbinv + XtHinv %*% X
  O <- solve(Oinv)
  VbinvC <- Vbinv %*% C
  M <- VbinvC + XtHinv %*% Y
  Snew <- S + crossprod(C, VbinvC) + crossprod(Y, Hinv) %*% Y - crossprod(M, O) %*% M
  vnew <- v + n
  Cnew <- O %*% M
  Vbnew <- O
  
  out[["vnew"]] <- vnew  
  out[["Snew"]] <- Snew
  out[["Cnew"]] <- Cnew
  out[["Vbnew"]] <- Vbnew
  return(out)
}