source("../R/gen_data.R")
source("../R/BS.R")
source("../R/pds.R")
source("../R/MNIW_sampler.R")
source("naiive_BS.R")
source("test_utils.R")
source("naiive_FF.R")
source("naiive_FFBS.R")
library(MNIW)
library(matrixsampling)
library(microbenchmark)
library(Rcpp)
sourceCpp("../src/FF.cpp")
sourceCpp("../src/BS.cpp")


seed <- 1234
set.seed(seed)
N <- 100
S <- 1000
p <- 20
nT <- 20
dt <- gen_ffbs_data(N = N, S = S, p = p, nT = nT)

# Generate FF data
res_ff <- FF_R_naiive(Y = dt$Y, F_ls = dt$para$F0, G_ls = dt$para$G0, 
                       W_ls = dt$para$W0, V_ls = dt$para$V0, 
                       m0 = dt$para$m0, M0 = dt$para$M0, 
                       n0 = dt$para$n0, D0 = dt$para$D0, 
                       nT = nT, delta = 1.0)

## BS_1step----
# compare mine with Daniel
res1 <- BS_1step_R_naiive(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                          st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                          at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                          Gt1 = dt$para$G0, delta = 1)

res2 <- smooth.1step(s = res_ff$T2$mt, S = res_ff$T2$Mt, 
                     m = res_ff$T1$mt, C = res_ff$T1$Mt, 
                     G = dt$para$G0, uRt = chol(res_ff$T2$At), delta = 1)

norm(res1$st - res2$s, "F")
norm(res1$St - res2$S, "F")


# Compare naiive R with optimized R
res3 <- BS_1step_R_naiive(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                          st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                          at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                          Gt1 = dt$para$G0, delta = 1)

res4 <- BS_1step_R(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                          st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                          at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                          Gt1 = dt$para$G0, delta = 1)

compare_R_cpp(res3, res4)
microbenchmark(BS_1step_R_naiive(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                                 st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                                 at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                                 Gt1 = dt$para$G0, delta = 1),
               BS_1step_R(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                                  st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                                  at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                                  Gt1 = dt$para$G0, delta = 1),
               smooth.1step(s = res_ff$T2$mt, S = res_ff$T2$Mt, 
                            m = res_ff$T1$mt, C = res_ff$T1$Mt, 
                            G = dt$para$G0, uRt = chol(res_ff$T2$At), delta = 1)
)

# Compare R with RcppEigen
res5 <- BS_1step_R(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                   st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                   at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                   Gt1 = dt$para$G0, delta = 1)

res6 <- BS_1step_cpp(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                   st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                   at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                   Gt1 = dt$para$G0, delta = 1)

compare_R_cpp(res5, res6)

microbenchmark(BS_1step_R(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                                  st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                                  at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                                  Gt1 = dt$para$G0, delta = 1),
               BS_1step_cpp(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
                                    st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
                                    at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
                                    Gt1 = dt$para$G0, delta = 1))

# Compare cpp with cpp_opt
# res7 <- BS_1step_cpp(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
#                    st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
#                    at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
#                    Gt1 = dt$para$G0, delta = 1)
# 
# res8 <- BS_1step_cpp_opt(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
#                      st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
#                      at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
#                      Gt1 = dt$para$G0, delta = 1)
# 
# compare_R_cpp(res7, res8)
# 
# microbenchmark(BS_1step_cpp(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
#                           st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
#                           at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
#                           Gt1 = dt$para$G0, delta = 1),
#                BS_1step_cpp_opt(mt = res_ff$T1$mt, Mt = res_ff$T1$Mt, 
#                             st1 = res_ff$T2$mt, St1 = res_ff$T2$Mt,
#                             at1 = res_ff$T2$at, At1 = res_ff$T2$At, 
#                             Gt1 = dt$para$G0, delta = 1))

## BS----
nsam <- 30
F_ls = dt$para$F0;
G_ls = dt$para$G0;
V_ls = dt$para$V0;
delta = 1.0
res_bs1 <- BS_R_naiive(nsam = nsam, res_ff, F_ls, V_ls, G_ls, nT = nT, delta = 1)
norm(res_bs1$Sigma[,,1] - dt$para$Sigma, "F")



## FFBS----
res_ffbs <- FFBS_R_naiive(nsam = nsam, Y = dt$Y, F_ls = dt$para$F0, G_ls = dt$para$G0, 
                      W_ls = dt$para$W0, V_ls = dt$para$V0, 
                      m0 = dt$para$m0, M0 = dt$para$M0, 
                      n0 = dt$para$n0, D0 = dt$para$D0, 
                      nT = nT, delta = 1.0)

