library(Rcpp)
sourceCpp("../cpp models/parallel_mvcalibratorDiscrepancy.cpp")

# initialize graph
N = 20
network_orig=network <- make_tree(N, children = 2,mode = "undirected")
plot(network)
l = layout_nicely(network)
mat=as_adjacency_matrix(network)
act = sample(c(1:10),size = N,replace = T)
act = rep(0,N)
act[1] = 10
initial_df <- data.frame(node=c(1,2:(N)), activation=act, stringsAsFactors=FALSE)
a = .01; b = .05; c = 0
results <- spreadr(mat, initial_df, include_t0 = T,
                   decay = a, retention = b, suppress = c)
results$activation=(results$activation+0)
dat = (reshape(results, idvar = "node", timevar = "time", direction = "wide")[,-1])

S = nrow(dat)
t_max = 10
times = (0:t_max)

FF = array(0, dim=c(S*sims, S, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(Matrix::bdiag(l))
}
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp)[3])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp[,,i])
}
st_time = Sys.time();
p2=1
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
st_time = Sys.time();
X = as.matrix(expand.grid(Gridx$x.mid,Gridx$x.mid)[locs,], ncol=2)
W=(.4*as.matrix(mat)+diag(S))
heatmap(W,Rowv = NA,Colv = NA)
fullySpatial = spDLMGP(z=t(z)[,-1], y=yin[,-1],
                 yy = array(ysp2[,-1,],
                            dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
                 FF=FF, FFc=FFsp,
                 m0=matrix(rep(0,S)), 
                 C0=W, W=W,
                 yinit = ysp2[1,1,], alpha0=.1, beta0=.1,
                 delta_v=.9, delta_w=.9,
                 params=theta_design, spparams = theta_design,
                 pred_params = test,
                 tune=rep(.001,2), tune2=rep(.1,2), 
                 tune3=1, tune4=rep(.001,2),
                 niter=500, burnin=500/2,
                 nugget=0, sigma_z=rep(1,1),
                 a0z = 1, b0z=1, S=S)

heatmap(apply(fullySpatial$mc_C,1:2,mean), Rowv = NA,Colv = NA)
matplot(t(fullySpatial$mc_beta),type="l")
(Sys.time()-st_time)
for(i in 1:ncol(theta_design)){
  matplot((fullySpatial$mc_beta[i,]),type="l",
          main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
}
matplot(fullySpatial$calibrate,type="l")
# hist((spgasp$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(spgasp$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=.2, col="red", lwd=3)
# apply(spgasp$calibrate, 2, acf, lag=200)
plot(sqrt(t(spgasp$mc_sigma2)),type="l")
# test
# apply(spgasp$calibrate,2,quantile)
test
save = fullySpatial$calibrate
round(colMeans(fullySpatial$calibrate),3)
{
  library(latex2exp)
  par(mfrow=c(2,2))
  for(i in 1:length(test)){
    drw = matrix(rbeta(length(fullySpatial$calibrate),1,1), ncol=ncol(fullySpatial$calibrate))
    if(i==1){
      hist(fullySpatial$calibrate[,i], breaks=10, main="", ylab="",
           col=rgb(0/255,0/255,.502,.7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(spgasp$calibrate[,i],, breaks=10, main="", ylab="",
           col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=1)
      hist(drw[,i], breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")),
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }else{
      hist(fullySpatial$calibrate[,i], breaks=10, main="",ylab="",
           col=rgb(0/255, 0/255, .502, .7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(spgasp$calibrate[,i],, breaks=10, main="", ylab="",
           col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=1)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }
    points(y=0,x=test[i], col="red", lwd=5,pch=19)
  }
  hist(sqrt(fullySpatial$mc_sigma2), xlim=c(0,.25), cex.lab=2,ylab="",
       col=rgb(0/255,0/255,.502,.7), breaks=10, main="",
       freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
  hist(sqrt(spgasp$mc_sigma2),, breaks=10, main="", ylab="",
       col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
       freq = 0, xlab=TeX(
         paste0('$\\eta_', i, "$")), add=1)
  hist(sqrt(rgamma(nrow(fullySpatial$calibrate),1,1)), xlim=c(0,2), breaks=10, freq = 0, xlab=TeX(
    paste0('$\\eta_', i, "$")),
    col=rgb(.678, .847, .902, alpha=.5),add=T)
  points(y=0,x=sigma_z, col="red", lwd=4,pch=19)
  plot(c(0,0), xlim=c(-10,10),ylim=c(-10,10),
       col="white",xlab="",ylab="",xaxt="n",yaxt="n",frame.plot = FALSE)
  legend(-13, 8, c("Prior", "Spatial Posterior", "Non-spatial Posterior"),
         fill=c(rgb(.678, .847, .902),
                rgb(0/255,0/255,.502, .7),
                rgb(1, 0, 0, alpha=.7)),
         bty = 'n',text.font=2,cex=1.5)
}
# points(rbind(test[1:2]),pch=19,cex=2,col="red")
# theta_loc = 7
# (test=theta_design[theta_loc,] + runif(4,0,0.05))
goback = c()
for(i in 1:ncol(theta_design)){
  goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
}
goback
# 
RcppParallel::setThreadOptions(numThreads = 4)
a = goback[1]
b = goback[2]
c = goback[3]
library(ggplot2)
result <- spreadr(pnet, start_run, time = t_max,
                  decay = a,
                  retention = b,
                  suppress = c,
                  include_t0=F)
# ggplot(result, aes(x=time, y=activation, color=node)) +
#   geom_point() +
#   geom_line()
