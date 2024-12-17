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

FF = array(0, dim=c(S*sims, 1, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(unlist(l))
}
nonspatial = DLMGP(z=z[-1,], y=yin[,-1],
               yy = array(ysp2[,-1,],
                          dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
               FF=FF, FFc=FFsp,
               m0=matrix(rep(0,1)), 
               C0=matrix(1), W=matrix(1),
               yinit = ysp[1,,1], alpha0=.1, beta0=.1,
               delta_v=.95, delta_w=.95,
               params=theta_design,
               pred_params = test,
               tune=rep(.001,2), tune2=rep(.001,2), 
               tune3=.25, tune4=rep(.01,2),
               niter=500, burnin=500/2,
               nugget=0, sigma_z=rep(1,1),
               a0z = 1, b0z=1, S=1, S2=ncol(z))

(Sys.time()-st_time)
for(i in 1:ncol(theta_design)){
  matplot((spgasp$mc_beta[i,]),type="l",
          main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
}
matplot(spgasp$calibrate,type="l")
# hist((spgasp$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(spgasp$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=.2, col="red", lwd=3)
# apply(spgasp$calibrate, 2, acf, lag=200)
# plot(t(spgasp$mc_sigma2))
# test
# apply(spgasp$calibrate,2,quantile)
test
round(colMeans(spgasp$calibrate),3)
{
  library(latex2exp)
  par(mfrow=c(2,2))
  for(i in 1:length(test[])){
    drw = matrix(rbeta(length(spgasp$calibrate),1,1), ncol=ncol(spgasp$calibrate))
    if(i==1){
      hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=15, 
           col=rgb(0/255,0/255,.502,.7),
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), main="Calibration Posterior", add=0)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="Calibration Posterior",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
      legend(.1, 6, c("Prior", "Posterior", "Truth"), 
             fill=c(rgb(.678, .847, .902),
                    rgb(0/255,0/255,.502, .7),
                    "red"))
    }else{
      hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=15, 
           col=rgb(0/255, 0/255, .502, .7),
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), main="", add=0)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }
    abline(v=test[i], col="red", lwd=3)
  }
  hist(sqrt(spgasp$mc_sigma2),
       col=rgb(0/255,0/255,.502,.7), breaks=30, main="",
       freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
  hist(rgamma(nrow(spgasp$calibrate),1,1), xlim=c(0,2), breaks=50, freq = 0, xlab=TeX(
    paste0('$\\eta_', i, "$")),
    col=rgb(.678, .847, .902, alpha=.5),add=T)
  abline(v=sigma_z, col="red", lwd=3)
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
