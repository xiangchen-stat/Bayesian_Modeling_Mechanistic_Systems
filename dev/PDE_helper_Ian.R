{
library(lhs)
library(mbr)
library(deSolve)
library(ReacTran)
library(latex2exp)

#NOTE: Ian's code here.
SI <- function(t, y, parms) {
  S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  dS <- -eta1*S*I/N
  dI <- eta1*S*I/N - eta2*I +
    tran.2D(C = I, D.x = alpha, D.y = alpha,
            dx = Gridx, dy = Gridy,
            C.x.up = 0,
            C.y.up = 0,
            C.x.down = 0,
            C.y.down = 0)$dC
  list(c(dS, dI))
}

SIR <- function(t, y, parms) {
  S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  R <- matrix(nrow = Nx, ncol = Ny, data = y[(2*Nx*Ny+1) : (3*Nx*Ny)])
  #NOTE: tran.2D as called here is equivalent to computing the Laplacian of S, times alpha1;
  #	with dC being the result.
  dS <- -eta1*S*I/N +
	  tran.2D(C = S, D.x = alpha1, D.y = alpha1,
	  dx = Gridx, dy= Gridy,
	  C.x.up = 0,
	  C.y.up = 0,
	  C.x.down = 0,
	  C.y.down = 0)$dC
  dI <- eta1*S*I/N - eta2*I +
    tran.2D(C = I, D.x = alpha2, D.y = alpha2,
            dx = Gridx, dy = Gridy,
            C.x.up = 0,
            C.y.up = 0,
            C.x.down = 0,
            C.y.down = 0)$dC
  dR <- eta2*I + 
	  tran.2D(C = R, D.x = alpha3, D.y = alpha3,
            dx = Gridx, dy = Gridy,
            C.x.up = 0,
            C.y.up = 0,
            C.x.down = 0,
            C.y.down = 0)$dC
  list(c(dS, dI, dR))
}

# output parameters
out.dir <- ".."

# NOTE: Set initial conditions for ODE here.
Nx <- 10
Ny <- 10
Gridx <- setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
Gridy <- setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)
N = 100000
X1ini <- matrix(nrow = Nx, ncol = Ny, data = N)
X2ini <- matrix(nrow = Nx, ncol = Ny, data = 0)
X2ini[5,5] <- 50
X3ini <- matrix(nrow = Nx, ncol = Ny, data = 0)
yini <- c(X1ini, X2ini, X3ini)
times <- seq(0,100)#,by=10)
# print(system.time(
#   out <- ode.2D(y = yini, parms = NULL, func = SI,
#                 nspec = 2, dimens = c(Nx, Ny), times = times,
#                 lrw = 2000000, names=c("X1", "X2"))
# ))
# par(oma = c(0,0,1,0))
# image(out, which = "X2", xlab = "x", ylab = "y",
#       mfrow = c(3, 3), ask = FALSE,
#       main = paste("t = ", times),
#       grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
# mtext(side = 3, outer = TRUE, cex = 1.25, line = -1,
#       "2-D SI")
temp = expand.grid(x = 1:Nx, y = 1:Nx)
spPts = 15
locs = temp[round(qunif(maximinLHS(spPts,k = 1),1,Nx^2)),]
locs = as.numeric(rownames(locs))
# TODO: sims=2 just to keep the code from breaking. We'll have to get rid of it later.
sims = 50 #25
cube = maximinLHS(sims,k = 5)
set.seed(1)
e1 = qunif(cube[,1],.5,1)
#set.seed(1)
e2 = qunif(cube[,2],0.01,.5)
set.seed(1)
a1 = qunif(cube[,3],0,.05)
a2 = qunif(cube[,4],0,.05)
a3 = qunif(cube[,5],0,.05)

# check the parameter estimates
print(sprintf("eta = (%.3f, %.3f), alpha = (%.3f, %.3f, %.3f)", e1[sims], e2[sims], a1[sims], a2[sims], a3[sims]))

S = length(locs)
ysp_old = array(0, dim=c(sims,length(times),length(locs)))
ysp = array(0, dim=c(sims,S,length(times)))
ysp2 = array(0, dim=c(sims,length(times),S))
FF1 = array(0, dim=c(sims, S, length(times)-1))
FFc = array(0, dim=c(sims, 1, length(times)-1, S))
FF = array(0, dim=c(sims, 1, length(times)-1, length(locs)))
#theta_design = cbind(b, g, a)
theta_design = cbind(e1, e2, a1, a2, a3)

#write out infected
write.csv(theta_design, 
          file=file.path(out.dir, "data", paste0("pde_para.csv", sep = "")), 
          row.names=FALSE)
#TODO: generate output per simulation
lrw <- Nx * Ny * length(times) * sims
i <- 1
for(i in 1:sims){
  if(i %% 5 == 0){
    print(Sys.time())
    print(paste(i, "/", sims))
  }
  eta1 = e1[i]
  eta2 = e2[i]
  alpha1 = a1[i]
  alpha2 = a2[i]
  alpha3 = a3[i]
  out <- ode.2D(y = yini, parms = c(eta1=eta1,eta2=eta2,
				    alpha1=alpha1, alpha2=alpha2, alpha3=alpha3), 
                func = SIR,
                nspec = 3, dimens = c(Nx, Ny), times = times,
                lrw = lrw, names=c("X1", "X2", "X3"))[,-1]

  #write out infected
  write.csv(out[,(1 + Nx * Ny):(2 * Nx * Ny)], 
            file=file.path(out.dir, "data", paste0("pde_solution_", i, ".csv", sep = "")), 
	      row.names=FALSE)

  # save the matrix of infected people in one row
  ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
  ysp_ = log(ysp_+1)
  for(sp in 1:length(locs)){
    ysp_old[i,,sp] = (ysp_[,locs[sp]])
    FF[i,1,,sp] = (ysp_old[i,-length(times),sp])
    ysp[i,sp,] = (ysp_[,locs[sp]])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    ysp2[i,,sp] = (ysp_[,locs[sp]])
    FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
  }
}

}
#TODO: what goes on past this point?

tempArr = array(0,dim=c(dim(FF)[1], dim(FF)[2], dim(FF)[3]))
FFsp = list()
# for each location
for(i in 1:dim(ysp_old)[3]){

  for(j in 1:dim(tempArr)[3]){
    tempArr[,,j] = FF[,1,j,i]
  }
  FFsp[[i]] = tempArr
}
#library(LICORS)
unscaled = theta_design
x = theta_design
for(i in 1:ncol(theta_design)){
  x[,i] = (theta_design[,i])
  temp = (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
  att.temp = attributes(temp)
  theta_design[,i] = temp
}
set.seed(1)
(test = c(runif(ncol(theta_design),.3,.8)))
goback = c()
for(i in 1:ncol(theta_design)){
  goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
}
RcppParallel::setThreadOptions(numThreads = "auto", stackSize = "auto")
RcppParallel::defaultNumThreads()
eta1 = goback[1]
eta2= goback[2]
alpha1= goback[3]
alpha2= goback[4]
alpha3= goback[5]

out <- ode.2D(y = yini, parms = list(alpha1=alpha1,
                                     alpha2=alpha2,
                                     alpha3=alpha3,
                                     eta1=eta1,
                                     eta2=eta2), func = SIR,
              nspec = 3, dimens = c(Nx, Ny), times = times,
              lrw = 2000000, names=c("X1", "X2", "X3"))[,-1]
#NOTE: record infected population I_t(s)
zz = log(out[,(Nx*Ny+1) : (2*Nx*Ny)]+1)
sigma_z = 1
set.seed(1)
z = zz[,locs] + matrix(rnorm(length(dim(ysp)[2]*dim(ysp)[3]), 0, sigma_z), nrow(zz), length(locs))
# matplot(z,type="l")
st_time = Sys.time();

FF = array(0, dim=c(S*sims, S, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(Matrix::bdiag(l))
}
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp_old)[2])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp_old[,i,])
}
calcSigma <- function(X,l) {
  Sigma <- matrix(rep(0, nrow(X)*nrow(X)), nrow=nrow(X))
  Sigma = diag(nrow(X))
  for(i in 1:(nrow(Sigma)-1)){
    for(j in (i+1):nrow(Sigma)){
      Sigma[i, j] = exp(-l * t(X[i, ] - X[j, ])%*%(X[i, ] - X[j, ]));
      Sigma[j, i] = Sigma[i, j];
    }
  }
  return(Sigma)
}
X = as.matrix(expand.grid(Gridx$x.mid,Gridx$x.mid)[locs,], ncol=2)
zsp = t(z)
X = as.matrix(expand.grid(Gridx$x.mid,Gridx$x.mid)[locs,], ncol=2)
# NOTE: spatial correlation matrix?
C0=calcSigma(X,.09)
heatmap(C0,Rowv = NA,Colv = NA)
# library(lhs)
# library(mbr)
# library(deSolve)
# library(ReacTran)
# library(latex2exp)
# # https://arxiv.org/pdf/1310.0505.pdf
# SI <- function(t, y, parms) {
#   S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
#   I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
#   dS <- -eta*S*I/N 
#   dI <- eta*S*I/N - gamma*I +
#     tran.2D(C = I, D.x = alpha1, D.y = alpha1,
#             dx = Gridx, dy = Gridy,
#             C.x.up = 0,
#             C.y.up = 0,
#             C.x.down = 0,
#             C.y.down = 0)$dC
#   list(c(dS, dI))
# }
# Nx <-5
# Ny <-5
# Gridx <- setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
# Gridy <- setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)
# eta=.6; gamma=.2; alpha1=1; alpha2 = 0
# 
# N = 100000
# X1ini <- matrix(nrow = Nx, ncol = Ny, data = N)
# X2ini <- matrix(nrow = Nx, ncol = Ny, data = 0) 
# X2ini[5,5] = 10
# 
# 
# yini <- c(X1ini, X2ini)
# times <- seq(0,100,by=10)
# print(system.time(
#   out <- ode.2D(y = yini, parms = NULL, func = SI,
#                 nspec = 2, dimens = c(Nx, Ny), times = times,
#                 lrw = 2000000, names=c("X1", "X2"))
# ))
# par(oma = c(0,0,1,0))
# image(out, which = "X2", xlab = "x", ylab = "y",
#       mfrow = c(3, 3), ask = FALSE,
#       main = paste("t = ", times),
#       grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
# mtext(side = 3, outer = TRUE, cex = 1.25, line = -1,
#       "2-D SI")
# temp = expand.grid(x = 1:Nx, y = 1:Nx)
# spPts = 15
# locs = temp[round(qunif(maximinLHS(spPts,k = 1),1,Nx^2)),]
# plot(locs, pch=19, col="red")
# locs = as.numeric(rownames(locs))
# sims=25
# cube = maximinLHS(sims,k = 4)
# set.seed(1)
# b = qunif(cube[,1],.5,1)
# set.seed(1)
# g = qunif(cube[,2],0.01,.5)
# set.seed(1)
# a = qunif(cube[,3],0,.01)
# 
# S = length(locs)
# times = seq(1,nrow(out), by=1)
# # ysp = array(0, dim=c(sims,length(times),length(locs)))
# ysp_old = array(0, dim=c(sims,length(times),length(locs)))
# ysp = array(0, dim=c(sims,S,length(times)))
# ysp2 = array(0, dim=c(sims,length(times),S))
# FF1 = array(0, dim=c(sims, S, length(times)-1))
# FFc = array(0, dim=c(sims, 1, length(times)-1, S))
# FF = array(0, dim=c(sims, 1, length(times)-1, length(locs)))
# theta_design = cbind(b, g, a1)
# plot(theta_design)
# for(i in 1:sims){
#   eta = b[i]
#   gamma = g[i]
#   alpha = a[i]
#   out <- ode.2D(y = yini, parms = NULL, func = SI,
#                 nspec = 2, dimens = c(Nx, Ny), times = times,
#                 lrw = 2000000, names=c("X1", "X2"))[,-1]
#   ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
#   # ysp_ = s[,2]
#   ysp_ = log(ysp_+1)
#   for(sp in 1:length(locs)){
#     ysp_old[i,,sp] = (ysp_[,locs[sp]])
#     FF[i,1,,sp] = (ysp_old[i,-length(times),sp])
#     
#     ysp[i,sp,] = (ysp_[,locs[sp]])
#     FF1[i,sp,] = (ysp[i,sp,-length(times)])
#     
#     ysp2[i,,sp] = (ysp_[,locs[sp]])
#     FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
#   }
# }
# tempArr = array(0,dim=c(dim(FF)[1], dim(FF)[2], dim(FF)[3]))
# FFsp = list()
# for(i in 1:dim(ysp_old)[3]){
#   for(j in 1:dim(tempArr)[3]){
#     tempArr[,,j] = FF[,1,j,i]
#   }
#   FFsp[[i]] = tempArr
# }
# library(LICORS)
# unscaled = theta_design
# x = theta_design
# for(i in 1:ncol(theta_design)){
#   x[,i] = (theta_design[,i])
#   temp = (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
#   att.temp = attributes(temp)
#   theta_design[,i] = temp
# }
# matplot(t(ysp_old[,,1]),add=0,type="l")
# for(i in 2:S){
#   matplot(t(ysp_old[,,i]),add=T, type="l")
# }
# set.seed(1)
# (test = c(runif(ncol(theta_design),.3,.7)))
# pairs(rbind(theta_design,test),
#       cex=c(rep(1, nrow(theta_design)), 2),
#       pch=c(rep(1, nrow(theta_design)), 19), 
#       col=c(rep("black",nrow(theta_design)),"red"),
#       labels = c(TeX('$\\eta$'),TeX('$\\gamma$'), TeX('$\\alpha$')))
# goback = c()
# for(i in 1:ncol(theta_design)){
#   goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
# }
# goback
# # 
# RcppParallel::setThreadOptions(numThreads = "auto", stackSize = "auto")
# RcppParallel::defaultNumThreads()
# eta = goback[1]
# gamma= goback[2]
# alpha1= goback[3]
# alpha2= goback[4]
# out <- ode.2D(y = yini, parms = list(alpha1=alpha1,
#                                      alpha2=alpha2,
#                                      eta=eta,
#                                      gamma=gamma), func = SI,
#               nspec = 2, dimens = c(Nx, Ny), times = times,
#               lrw = 2000000, names=c("X1", "X2"))[,-1]
# zz = log(out[,(Nx*Ny+1) : (2*Nx*Ny)]+1)
# sigma_z = 1
# set.seed(1)
# z = zz[,locs] + matrix(rnorm(length(dim(ysp)[2]*dim(ysp)[3]), 0, sigma_z), nrow(zz), length(locs))
# st_time = Sys.time();
# p2 = 1
