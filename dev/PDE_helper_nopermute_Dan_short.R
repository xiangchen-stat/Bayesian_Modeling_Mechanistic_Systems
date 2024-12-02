PC <- "mac"
switch (PC,
        "x1c" = setwd("D:/Documents/UCLA/0-Administrative/GSR/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "mac" = setwd("/Users/xiangchen/Documents/Bayesian-Modeling-Mechanistic-Systems/dev"),
        "xps" = setwd("C:/Users/wilson/Desktop/X1C/MNIW/dev")
)
library(deSolve)
library(ReacTran)

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

Nx <- 6
Ny <- 6
N = 10000
times <- seq(0, 5)

eta1 = runif(1, 2, 4)
eta2 = runif(1, .2, .4)
alpha1 = runif(1, .2, .4)
alpha2 = runif(1, .2, .4)
alpha3 = runif(1, .2, .4)

# debug
# eta1 = input_cal[1]
# eta2 = input_cal[2]
# alpha1 = input_cal[3]
# alpha2 = input_cal[4]
# alpha3 = input_cal[5]

parms <- c(eta1 = eta1, eta2 = eta2, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3)


# generate grid
locs <- expand.grid(y=1:Ny, x=1:Nx)
locs <- locs[,c("x", "y")]
locs = as.numeric(rownames(locs))
S = length(locs)
Gridx <- ReacTran::setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
Gridy <- ReacTran::setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)

# set up initial values for SIR
S0 <- matrix(nrow = Nx, ncol = Ny, data = N)
I0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
I0[Nx/2,Ny/2] <- 1 # set midpoint to have 50 infected people
I0[round(Nx/10), round(Ny/10)] <- N * 0.05
R0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
yini <- c(S0, I0, R0)


# Generate SIR data from the SIR ODE.
out <- ode.2D(y = yini, parms = parms,
              func = SIR,
              nspec = 3, dimens = c(Nx, Ny), times = times,
              lrw = 200000000, names=c("X1", "X2", "X3"))[,-1]
lgout <- log(out + 1)
res <- lgout[,(1 + Nx * Ny):(2 * Nx * Ny)]
    
print(dim(res))
