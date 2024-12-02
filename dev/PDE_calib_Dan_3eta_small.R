{
  PC <- "mac"
  switch (PC,
          "x1c" = setwd("D:/Documents/UCLA/0-Administrative/GSR/Bayesian-Modeling-Mechanistic-Systems/dev"),
          "mac" = setwd("/Users/xiangchen/Documents/Bayesian-Modeling-Mechanistic-Systems/dev"),
          "xps" = setwd("C:/Users/wilson/Desktop/X1C/MNIW/dev")
  )
  library(readr)
  library(lhs)
  library(deSolve)
  library(ReacTran)
  
  source("../dev-Dan/generate_grid.R")
  
  p = 2
  # NOTE: Nx is length of x-dimension, or rows, Ny is the length of the y-dimension, or columns.
  Nx <- 6
  Ny <- 6
  N = 10000
  nsims <- 50
  times <- seq(0, 25)
  
  
  bnrow <- nsims
  bncol <- 25
  
  locs <- expand.grid(y=1:Ny, x=1:Nx)
  locs <- locs[,c("x", "y")]
  locs = as.numeric(rownames(locs))
  
  S = length(locs)
  
  
  top.dir <- ".."
  out.dir <- ".."
  # train.out.dir <- file.path(top.dir, sprintf("data/SIR_nopermute_train_%dx%d", Nx, Ny))
  train.out.dir <- file.path(top.dir, sprintf("data/Dan_pde_nopermute_train_%dx%d_%d_%d_%d", 
                                              Nx, Ny, length(times), nsims, 1))
  dir.create(train.out.dir, showWarnings=FALSE)
  # train.out.head <- "SIR_nopermute_train"
  train.out.head <- "pde_reform"
  print(train.out.dir)
  # parameter_fname_train <- file.path(train.out.dir, "parameters_train.csv")
  parameter_fname_train <- file.path(train.out.dir, "pde_para.csv")
  
  e1 <- NULL
  e2 <- NULL
  a1 <- NULL
  a2 <- NULL
  a3 <- NULL
  
  
  print("Initializing parameters.")
  
  set.seed(123)
  
  cube = maximinLHS(nsims,k = 5)
  # write this cube
  cube_df <- as.data.frame(cube)
  names(cube_df) <- c("eta1", "eta2", "alpha1", "alpha2", "alpha3")
  cube_fname <- file.path(train.out.dir, "hypercube_samples_train.csv")
  # write.csv(cube_df, file=cube_fname, row.names=FALSE)
  
  e1 = qunif(cube[,1], 2, 4) # infection rate 
  e2 = qunif(cube[,2], 0.2, 0.4) # recovery rate
  
  # set.seed(1)
  # a1 = qunif(cube[,3],0,.2)
  a2 = qunif(cube[,4], 0, 0.2) # spatial infection rate
  # a3 = qunif(cube[,5],0,.2) # recovery rate
  a1 = rep(0, nsims)
  # a2 = rep(0, nsims)
  a3 = rep(0, nsims)
  
  # set some specific values in test data
  e1[nsims] <- 3
  e2[nsims] <- 0.3
  a2[nsims] <- 0.1
  
  write.csv(data.frame(eta1 = e1, eta2 = e2, alpha1 = a1, alpha2 = a2, alpha3 = a3),
            file=parameter_fname_train, row.names=FALSE)
  
  
  
  Gridx <- setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
  Gridy <- setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)
  
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
  
  S0 <- matrix(nrow = Nx, ncol = Ny, data = N)
  I0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  I0[Nx/2,Ny/2] <- 1 # set midpoint to have 50 infected people
  I0[round(Nx/10), round(Ny/10)] <- N * 0.05
  R0 <- matrix(nrow = Nx, ncol = Ny, data = 0)
  yini <- c(S0, I0, R0)
  
  
  eta1_list <- e1
  eta2_list <- e2
  alpha1_list <- a1
  alpha2_list <- a2
  alpha3_list <- a3
  
  # permutes data so that every block is rearranged so that every bncol slice 
  #   of the block captures a sqrt(bncol) x sqrt(bncol) spatial block of the data.
  # The number of rows is assumed to be trivial.
  permute_data <- function(data, Nx, Ny, bncolx, bncoly, grid.traversal="rowsnake") {
    fnrow <- bnrow <- nrow(data)
    out <- matrix(nrow=bnrow, ncol=ncol(data), data=0)
    
    # NOTE: bncol is the total number of squares in the spatial block. bncolx and bncoly is simply how it is partitioned.
    bncol <- bncolx * bncoly
    
    # Generate a grid for the data grid irrespective of partitioning.
    data.grid <- generate_grid(Ny, Nx, bncoly, bncolx,
                               traversal.mode = grid.traversal,
                               is.flexible=FALSE)
    
    K <- Nx * Ny / bncol
    
    # number of slices
    for (k in 1:K) {
      block_ixs <- data.grid[k,]
      for (cix in 1:bncolx) {
        # Get the corresponding column indices in the data.
        # out's indices proceed linearly.
        # data's pseudocode: (left_ix - 1) * Ny + up_ix:down_ix.
        out[,((k-1)*bncol + (cix-1)*bncoly +1):((k-1)*bncol + cix*bncoly)] <- 
          data[,((block_ixs$L + cix - 2) * Ny + block_ixs$U):((block_ixs$L + cix - 2) * Ny + block_ixs$D)]
      }
    }
    
    return(out)
  }
  
  
  # Generate SIR data from the SIR ODE.
  # If last Y file and last F file do not exist, generate the data.
  if (!file.exists(file.path(train.out.dir, sprintf("%s_%d.csv", train.out.head, max(times)))) | 
      !file.exists(file.path(train.out.dir, sprintf("%s_F%d.csv", train.out.head, max(times))))) {
    
    #gen_SIR_data <- function(sims, eta1_list, eta2_list,
    #			 alpha1_list, alpha2_list, alpha3_list, yini,
    #			 Nx, Ny, times, locs, out.head = "SIR") {
    
    ysp = array(0, dim=c(nsims,S,length(times)))
    
    for(i in 1:nsims){
      if((i %% round(nsims*0.1) == 0) || nsims <= 10){
        print(paste("Iteration", i, ":", Sys.time()))
      }
      eta1 = eta1_list[i]
      eta2 = eta2_list[i]
      alpha1 = alpha1_list[i]
      alpha2 = alpha2_list[i]
      alpha3 = alpha3_list[i]
      
      #TODO: The parms aren't making it in if the stuff is warpped in a function.
      parms <- c(eta1 = eta1, eta2 = eta2, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3)
      
      out <- ode.2D(y = yini, parms = parms,
                    func = SIR,
                    nspec = 3, dimens = c(Nx, Ny), times = times,
                    lrw = 200000000, names=c("X1", "X2", "X3"))[,-1]
      lgout <- log(out + 1)
      # write out infected
      write.csv(lgout[,(1 + Nx * Ny):(2 * Nx * Ny)], 
                file=file.path(train.out.dir, paste0("/pde_solution_", i, ".csv", sep = "")), 
                row.names=FALSE)
      
      #  if (i == nsims) write.table(out, file=file.path(train.out.dir, "SIR_solution.txt"), 
      #	      sep=",", row.names=FALSE, col.names=FALSE)
      
      # save the matrix of infected people in one row
      ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
      ysp_ = log(ysp_+1)
      
      for(sp in 1:length(locs)){
        # NOTE: Get I_t(s) at s = locs[sp]
        #ysp_old[i,,sp] = (ysp_[,locs[sp]])
        # NOTE: Ian sets the FF at each spatial location to be the log (I_t(s) + 1) at that point in space across time?
        #FF[i,1,,sp] = (ysp_old[i,,sp])
        ysp[i,sp,] = (ysp_[,locs[sp]])
        # NOTE: set FF1 at each simulation and spatial coordinate the same way he set FF? 
        #FF1[i,sp,] = (ysp[i,sp,-length(times)])
        #ysp2[i,,sp] = (ysp_[,locs[sp]])
        #FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
      }
    }
  }
}
