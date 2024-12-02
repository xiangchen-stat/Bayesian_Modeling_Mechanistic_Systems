
gen_F_ls_AR1_stdrow <- function(N, nT, Y){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft_temp <- 1/2 * (Y[[1]] + Y[[2]])
      Ft <- t(apply(Ft_temp, MARGIN = 1, FUN = scale))
      F_ls[[i]] <- Ft
    } else{
      Ft_temp <- Y[[i - 1]]
      Ft <- t(apply(Ft_temp, MARGIN = 1, FUN = scale))
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}

gen_F_ls_AR1_stdall <- function(N, nT, Y){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft_temp <- 1/2 * (Y[[1]] + Y[[2]])
      mu_Ft <- mean(Ft_temp)
      sd_Ft <- sd(as.vector(Ft_temp))
      Ft <- (Ft_temp - mu_Ft) / sd_Ft
      F_ls[[i]] <- Ft
    } else{
      Ft_temp <- Y[[i - 1]]
      mu_Ft <- mean(Ft_temp)
      sd_Ft <- sd(as.vector(Ft_temp))
      Ft <- (Ft_temp - mu_Ft) / sd_Ft
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}

gen_F_ls_AR1_stdall_var <- function(N, nT, Y){
  F_ls <- list()
  for (i in 1:nT) {
    if(i == 1){
      Ft_temp <- 1/2 * (Y[[1]] + Y[[2]])
      sd_Ft <- sd(as.vector(Ft_temp))
      Ft <- Ft_temp / sd_Ft
      F_ls[[i]] <- Ft
    } else{
      Ft_temp <- Y[[i - 1]]
      sd_Ft <- sd(as.vector(Ft_temp))
      Ft <- Ft_temp / sd_Ft
      F_ls[[i]] <- Ft
    }
  }
  return(F_ls)
}