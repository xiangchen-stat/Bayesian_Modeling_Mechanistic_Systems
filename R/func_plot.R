library(ggpubr)
# p_load(MBA, reshape2, ggmap, sf) # for spatial analysis

# Set ggplot theme
theme_set(theme_minimal(base_size = 22))
col_epa <- c("#00e400", "#ffff00", "#ff7e00", "#ff0000", "#99004c", "#7e0023")
col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
             "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429",
             "#a62021", "#871b1c")

quick_heat <- function(dt, max_y){
  p <- ggplot(dt, aes(x = x, y = y, fill = fill)) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(0, max_y),
                         oob = scales::squish) +
    labs(x = "x", y = "y", fill = "Value")+
    # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
    theme(text = element_text(size=28),
          legend.text = element_text(size = 28),
          legend.key.size = unit(1.5, "cm"))
  return(p)
}
quick_save <- function(filename, path_fig, plot){
  ggsave(filename = paste(filename, "_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = plot,
         device = "png",
         width = 42,
         height = 35,
         units = "cm",
         dpi = 100
  ) 
}

# Emulation ----
# plot heatmap as a 3x3 panel
plot_panel_heatmap_9 <- function(dat, input_num, tstamp, max_y, Nx, Ny, nT){
  plot_ls <- list()
  # tstamp <- as.integer(seq(1, nT, length.out = 10))
  ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  ind_plot <- 1
  
  for (i in tstamp) {
    if(i == 1){
      ind_plot <- 1 # start counting
    }
    temp <- dat[[i]][input_num,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol = temp)%>% 
      as.data.frame()
    
    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value")+
      # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
      # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
      theme(text = element_text(size=28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(1.5, "cm"))
    
    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }
  ggsave(filename = paste("plot_panel_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                          plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                          plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                          ncol = 3, nrow = 3,
                          labels = c(paste("PDE: t =", tstamp[1]-1),
                                     paste("PDE: t =", tstamp[2]-1), 
                                     paste("PDE: t =", tstamp[3]-1), 
                                     paste("PDE: t =", tstamp[4]-1),
                                     paste("PDE: t =", tstamp[5]-1), 
                                     paste("PDE: t =", tstamp[6]-1), 
                                     paste("PDE: t =", tstamp[7]-1),
                                     paste("PDE: t =", tstamp[8]-1), 
                                     paste("PDE: t =", tstamp[9]-1)),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          # hjust = -1,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 50,
         units = "cm",
         dpi = 100
  )
  return(plot_ls)
}


cal_errorbar <- function(X){
  out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = median),
                    lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
                    upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
  return(out)
}


cal_errorbar_mean <- function(X){
  out <- data.frame(med = apply(X = X, MARGIN = 1, FUN = mean),
                    lower = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.025),
                    upper = apply(X = X, MARGIN = 1, FUN = quantile, prob = 0.975))
  return(out)
}

# calibration ----
# plot heatmap as a 3x3 panel
plot_panel_heatmap_9_cal <- function(dat, tstamp, max_y, loc_cal, Nx, Ny){
  plot_ls <- list()
  ind_plot <- 1 # start counting
  
  for (i in tstamp) {
    temp <- dat[i,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = loc_cal$row, col = loc_cal$col, sol = temp)%>% 
      as.data.frame()
    
    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value")+
      # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
      # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
      theme(text = element_text(size=28),
            legend.text = element_text(size = 28),
            legend.key.size = unit(1.5, "cm"))
    
    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }
  ggsave(filename = paste("plot_panel_cal_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                          plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                          plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                          ncol = 3, nrow = 3,
                          labels = c(paste("PDE: t =", tstamp[1]-1),
                                     paste("PDE: t =", tstamp[2]-1), 
                                     paste("PDE: t =", tstamp[3]-1), 
                                     paste("PDE: t =", tstamp[4]-1),
                                     paste("PDE: t =", tstamp[5]-1), 
                                     paste("PDE: t =", tstamp[6]-1), 
                                     paste("PDE: t =", tstamp[7]-1),
                                     paste("PDE: t =", tstamp[8]-1), 
                                     paste("PDE: t =", tstamp[9]-1)),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          # hjust = -1,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 50,
         units = "cm",
         dpi = 100
  )
  return(plot_ls)
}
