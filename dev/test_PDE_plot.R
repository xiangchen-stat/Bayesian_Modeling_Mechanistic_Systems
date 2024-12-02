{
# need change: n_input, Nx, Ny, nsam
setwd("D:/Documents/UCLA/0-Administrative/GSR/MNIW/dev")
source("init_lib.r")
seed <- 1234
set.seed(seed)
path <- "Dan_pde_nopermute_train_100x100_51_20_1"
path_data <- file.path("..", "data", path)
path_fig <- file.path("..", "figures", path)

# read in pde data----
ind_old_data <- F
n_input <- 20
nT <- 51
Nx <- 100
Ny <- 100
N_people <- 100000
N_sp <- Nx * Ny
prop_train <- 0.5
n_train <- round(n_input * prop_train)
n_test <- n_input - n_train
pde_para <- as.matrix(read.csv(file = paste0(path_data, "/pde_para.csv", sep = "")))[1:n_input,]
pde_para_train <- pde_para[1:n_train,]
pde_para_test <- pde_para[(n_train + 1):n_input,]
ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))

if(ind_old_data == TRUE){
  # read old data
  load(file.path(paste(path_data, "/dat_pde.RData", sep = "")))
  load(file.path(paste(path_data, "/dt_pde_train.RData", sep = "")))
  load(file.path(paste(path_data, "/dt_pde_test.RData", sep = "")))
} else{
  # manipulate new data
  dat_pde <- list()
  for (i in 1:n_input) {
    dat_pde[[i]] <- as.matrix(read.csv(file = paste0(path_data, "/pde_solution_", i, ".csv", sep = "")))
    if(i %% 10 == 0){
      print(paste("read in pde_solution", i, "/", n_input))
      print(Sys.time())
    }
  }
  save(dat_pde, file = file.path(paste(path_data, "/dat_pde.RData", sep = "")))
  
  # transform the pde train data
  # read in transformed train data
  dt_pde_train <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in 1:n_train) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_train[[i]] <- temp
    write.csv(temp,
              file=file.path(path_data, paste0("pde_train_", i, ".csv", sep = "")),
              row.names=FALSE)
  }
  names(dt_pde_train) <- paste0("T", seq(1:nT))
  save(dt_pde_train, file = file.path(paste(path_data, "/dt_pde_train.RData", sep = "")))
  
  # transform the pde test data
  # read in transformed test data
  dt_pde_test <- list()
  for (i in 1:nT) {
    temp <- c()
    for (j in (n_train + 1):n_input) {
      temp <- rbind(temp, dat_pde[[j]][i,])
    }
    dt_pde_test[[i]] <- temp
    write.csv(temp,
              file=file.path(path_data, paste0("pde_test_", i, ".csv", sep = "")),
              row.names=FALSE)
  }
  names(dt_pde_test) <- paste0("T", seq(1:nT))
  save(dt_pde_test, file = file.path(paste(path_data, "/dt_pde_test.RData", sep = "")))
}

## Plot PDE results ----
# ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
# plot_ls <- list() 
# tstamp <- as.integer(seq(1, nT, length.out = 9))

for (i in 1:n_train) {
  dat <- dt_pde_train
  input_num <- i
  tstamp <- as.integer(seq(1, nT, length.out = 9))
  max_y <- max(as.vector(unlist(dat))) # set max limit for all plots
  
  pde_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
}

for (i in 1:n_test) {
  dat <- dt_pde_test
  input_num <- i
  tstamp <- as.integer(seq(1, nT, length.out = 9))
  max_y <- max(as.vector(unlist(dat))) # set max limit for all plots
  
  pde_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
}

}

