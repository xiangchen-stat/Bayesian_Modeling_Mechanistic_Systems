{
source("../R/gen_data.R")
source("../R/FF.R")
source("../R/MNIW_sampler.R")
source("naiive_FF.R")
source("naiive_BS.R")
source("naiive_FFBS.R")
source("test_utils.R")
source("../R/generate_grid.R")
source("../R/big_data_set.R")
library(MNIW)
library(matrixsampling)
library(microbenchmark)
library(Rcpp)
sourceCpp("../src/FF.cpp")
sourceCpp("../src/BS.cpp")


seed <- 1234
set.seed(seed)
N <- 100
S <- 100
p <- 2
nT <- 5
nsam <- 5
fnrow <- N
fncol <- S
bnrow <- 25 # block nrow
bncol <- 25

path_out <- "D:/Documents/UCLA/0-Administrative/GSR/MNIW/out"
data.out.head <- file.path(path_out, "data_")
gen_ffbs_csv(N = N, S = S, p = p, nT = nT, path = data.out.head)

Y_filename_list <- paste0(data.out.head, "Y", seq(1, nT), ".csv")
F_filename <- paste0(data.out.head, "F0.csv")
G_filename <- paste0(data.out.head, "G0.csv")
W_filename <- paste0(data.out.head, "W0.csv")
V_filename <- paste0(data.out.head, "V0.csv")

m0_filename <- paste0(data.out.head, "m0.csv")
M0_filename <- paste0(data.out.head, "MM0.csv")
# n0_filename <- paste0(data.out.head, "n0.csv")
D0_filename <- paste0(data.out.head, "Sigma.csv")

# dt <- gen_ffbs_data(N = bnrow, S = bncol, p = p, nT = nT)


## FF----
Y_ls = Y_filename_list
F_ls = F_filename
V_ls = V_filename 
G_ls = G_filename
W_ls = W_filename
m0 = as.matrix(read_big_csv_quick(m0_filename, cols = c(1, bncol)))
M0 = as.matrix(read_big_csv_quick(M0_filename))
D0 = as.matrix(read_big_csv_quick(D0_filename, rows = c(1, bncol), cols = c(1, bncol)))
n0 = dim(D0)[1] + 1

# F_ls = dt$para$F0;
# G_ls = dt$para$G0;
# W_ls = dt$para$W0;
# V_ls = dt$para$V0;
# m0 = dt$para$m0;
# M0 = dt$para$M0;
# n0 = dt$para$n0;
# D0 = dt$para$D0
nT = nT; fnrow = fnrow; fncol = fncol;
bnrow = bnrow; bncol = bncol;
delta = 1.0
i <- 1
j <- 1
f <- 1
}
Sys.time()
## FF big data----
FF_bigdata_R_naiive(Y_ls = Y_filename_list, F_ls = F_ls, V_ls = V_ls, 
                    G_ls = G_ls, W_ls = W_ls, 
                     m0 = m0, M0 = M0, 
                     n0 = n0, D0 = D0,  
                     fnrow = fnrow, fncol = fncol,
                     bnrow = bnrow, bncol = bncol, path_out = path_out,
                     nT = nT, delta = 1.0, verbose = TRUE)


## BS big data----
{
FF.out.head <- file.path(path_out, "output_FF_")
nt_ls <- paste0(FF.out.head, "nt_F", seq(1, nT), ".csv")
Dt_ls <- paste0(FF.out.head, "Dt_F", seq(1, nT), ".csv")
mt_ls <- paste0(FF.out.head, "mt_F", seq(1, nT), ".csv")
Mt_ls <- paste0(FF.out.head, "MMt_F", seq(1, nT), ".csv")
at_ls <- paste0(FF.out.head, "at_F", seq(1, nT), ".csv")
At_ls <- paste0(FF.out.head, "AAt_F", seq(1, nT), ".csv")
FT <- paste0(FF.out.head, "Ft_F", nT, ".csv")
VT <- paste0(FF.out.head, "vt_F", nT, ".csv")

f <- nT
# i <- n_b
}
Sys.time()
BS_bigdata_R_naiive(nt_ls = nt_ls, Dt_ls = Dt_ls, 
                    mt_ls = mt_ls, Mt_ls = Mt_ls, at_ls = at_ls, At_ls = At_ls,
                    FT = FT, VT = VT, G_ls = G_ls, 
                    m0 = m0, M0 = M0, 
                    fnrow = fnrow, fncol = fncol,
                    bnrow = bnrow, bncol = bncol,
                    path_out = path_out, 
                    nT = nT, delta = 1, verbose = T)

Sys.time()
FFBS_bigdata_R_naiive(Y_ls = Y_filename_list, F_ls = F_ls, V_ls = V_ls, 
                      G_ls = G_ls, W_ls = W_ls, 
                      m0 = m0, M0 = M0, 
                      n0 = n0, D0 = D0,  
                      fnrow = fnrow, fncol = fncol,
                      bnrow = bnrow, bncol = bncol, path_out = path_out,
                      nT = nT, delta = 1.0, verbose = TRUE)

Sys.time()

# Predict----
{
n_new <- 5
Ft_new <- gen_ran_matrix(nrow = n_new, ncol = p)
var_aug <- gen_pd_matrix(dim = bnrow + n_new)
Vt <- var_aug[1:bnrow, 1:bnrow]
Vt_new <- var_aug[(bnrow+1):(bnrow+n_new), (bnrow+1):(bnrow+n_new)]
Jt <- var_aug[(1:bnrow), (bnrow+1):(bnrow+n_new)]
write.csv(x = Ft_new, file = paste0(path_out, "/output_fix_Ft_new.csv", sep = ""), row.names = F)
write.csv(x = Vt, file = paste0(path_out, "/output_fix_Vt.csv", sep = ""), row.names = F)
write.csv(x = Vt_new, file = paste0(path_out, "/output_fix_Vt_new.csv", sep = ""), row.names = F)
write.csv(x = Jt, file = paste0(path_out, "/output_fix_Jt.csv", sep = ""), row.names = F)

Sigma_ls <- paste0(path_out, "/output_BS_Sigma.csv", sep = "")
Theta_ls <- paste0(path_out, "/output_BS_theta_F", seq(1, nT), ".csv", sep = "")
F_new_ls <- paste0(path_out, "/output_fix_Ft_new.csv", sep = "")
V_ls <- paste0(path_out, "/output_fix_Vt.csv", sep = "")
V_new_ls <- paste0(path_out, "/output_fix_Vt_new.csv", sep = "")
J_ls <- paste0(path_out, "/output_fix_Jt.csv", sep = "")

f <- 1
}
FFBS_predict_bigdata_R_naiive(Sigma_ls = Sigma_ls, Theta_ls = Theta_ls, Y_ls = Y_ls,
                              F_ls, F_new_ls = F_new_ls, J_ls = J_ls,  
                              V_ls = V_ls, V_new_ls = V_new_ls, 
                              fnrow = fnrow, fncol = fncol,
                              bnrow = bnrow, bncol = bncol, path_out = path_out,
                              nT = nT, delta = 1.0, verbose = TRUE)




# # compare
# # Daniel's para
# dt_Dan <- list()
# dt_Dan$n0 <- 2
# dt_Dan$D0 <- diag(c)
# dt_Dan$m0 <- matrix(rep(0, p*c), nrow=p)
# dt_Dan$M0 <- diag(p)
# dt_Dan$Gtmat <- matrix(diag(p), nrow=1)
# dt_Dan$Wtmat <- matrix(diag(p), nrow=1)
# dt_Dan$Vtmat <- diag(bnrow)
# dt_Dan$F0 <- F_filename
# dt_Dan$V0 <- V_filename
