source("../R/gen_data.R")
source("../R/FF.R")
source("../R/pds.R")
source("naiive_FF.R")
source("test_utils.R")
library(MNIW)
library(matrixsampling)
library(microbenchmark)
library(Rcpp)
sourceCpp("../src/FF.cpp")


seed <- 1234
set.seed(seed)
N <- 100
S <- 100
p <- 2
nT <- 20
dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)


## FF_1step----
# compare mine with Daniel
res1 <- FF_1step_R_naiive(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                          Wt = dt$para$W0, Vt = dt$para$V0, 
                          mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                          nt_1 = dt$para$n0, Dt_1 = dt$para$D0)

res2 <- FF.1step(Y = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, m = dt$para$m0, 
                 C = dt$para$M0, W = dt$para$W0, V = dt$para$V0, a = dt$para$n0,
                 B = dt$para$D0, n = dim(dt$para$V0)[1], delta = 1)

res1$nt - res2$a
res1$delta - res2$delta
norm(res1$Dt - res2$B, "F")
norm(res1$Qt - crossprod(res2$uQt,res2$uQt), "F")
norm(res1$qt - res2$f, "F")
norm(res1$at - res2$mstar, "F")
norm(res1$At - crossprod(res2$uRt,res2$uRt), "F")
norm(res1$mt - res2$m, "F")
norm(res1$Mt - res2$C, "F")

# Compare naiive R with optimized R
res3 <- FF_1step_R_naiive(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                          Wt = dt$para$W0, Vt = dt$para$V0, 
                          mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                          nt_1 = dt$para$n0, Dt_1 = dt$para$D0)

res4 <- FF_1step_R(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                   Wt = dt$para$W0, Vt = dt$para$V0, 
                   mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                   nt_1 = dt$para$n0, Dt_1 = dt$para$D0)

compare_R_cpp(res3, res4)


# Compare R with RcppEigen
res5 <- FF_1step_R(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                   Wt = dt$para$W0, Vt = dt$para$V0, 
                   mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                   nt_1 = dt$para$n0, Dt_1 = dt$para$D0)

res6 <- FF_1step_cpp(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                     Wt = dt$para$W0, Vt = dt$para$V0, 
                     mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                     nt_1 = dt$para$n0, Dt_1 = dt$para$D0)

compare_R_cpp(res5, res6)

microbenchmark(FF_1step_R(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                          Wt = dt$para$W0, Vt = dt$para$V0, 
                          mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                          nt_1 = dt$para$n0, Dt_1 = dt$para$D0),
               FF_1step_cpp(Yt = dt$Y$Y1, Ft = dt$para$F0, Gt = dt$para$G0, 
                            Wt = dt$para$W0, Vt = dt$para$V0, 
                            mt_1 = dt$para$m0, Mt_1 = dt$para$M0, 
                            nt_1 = dt$para$n0, Dt_1 = dt$para$D0))

## FF----
# Y = dt$Y
# F_ls = dt$para$F0; 
# G_ls = dt$para$G0; 
# W_ls = dt$para$W0; 
# V_ls = dt$para$V0; 
# m0 = dt$para$m0; 
# M0 = dt$para$M0; 
# n0 = dt$para$n0;
# D0 = dt$para$D0
# delta = 1.0

res_ff1 <- FF_R_naiive(Y = dt$Y, F_ls = dt$para$F0, G_ls = dt$para$G0, 
                       W_ls = dt$para$W0, V_ls = dt$para$V0, 
                       m0 = dt$para$m0, M0 = dt$para$M0, 
                       n0 = dt$para$n0, D0 = dt$para$D0, 
                       nT = nT, delta = 1.0)

mean_iw <- res_ff1$T20$Dt / (res_ff1$T20$nt - dim(res_ff1$T20$Dt)[1] - 1)
norm(mean_iw - dt$para$Sigma, "F")

ff.file.paths <- FF.bigdata(~ V1 + V2, test_data_set, n0, D0, 
                            m0, M0, Gtmat, Wtmat, Vt = diag,
                            out.head = ff.out.head) 
