## Set up for plotting ----
source("../R/func_plot.R")

## Plot PDE results ----
### Heat Map----
{
  # ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  # plot_ls <- list() 
  # tstamp <- as.integer(seq(1, nT, length.out = 9))
  input_num <- 2
  tstamp <- as.integer(seq(1, nT_ori, length.out = 9))
  dat <- dt_pde_test
  max_y <- max(as.vector(unlist(dat))) # set max limit for all plots
  # max_y <- max(as.vector(unlist(res_pre_exact))) # set max limit for all plots
  
  pde_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  # dat <- postm_ls # this is monte carlo
  dat <- res_pre_exact # this is exact
  ffbs_heat <- plot_panel_heatmap_9(dat = dat, tstamp = tstamp,
                                    input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  # ffbs_heat19 <- plot_panel_heatmap_9(dat = dat, tstamp = 1:9,
  #                                   input_num = input_num, max_y = max_y, Nx = Nx, Ny = Ny, nT = nT)
  
  ggsave(filename = paste("plot_heat_compare_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = ggarrange(pde_heat[[4]], pde_heat[[6]], pde_heat[[8]],
                          ffbs_heat[[4]], ffbs_heat[[6]], ffbs_heat[[8]],
                          ncol = 3, nrow = 2,
                          labels = c(paste("PDE solution: t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1),
                                     paste("FFBS emulation: t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1)
                          ),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          hjust = -0.5,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 35,
         units = "cm",
         dpi = 100
  )
  
}

{
  ### quantile heat map for sigma----
  mc_sigma <- res_ffbs$Sigma
  sigma_low <- apply(mc_sigma, MARGIN = c(1, 2), FUN = quantile, probs = 0.025, na.rm = T)
  sigma_med <- apply(mc_sigma, MARGIN = c(1, 2), FUN = quantile, probs = 0.5, na.rm = T)
  sigma_high <- apply(mc_sigma, MARGIN = c(1, 2), FUN = quantile, probs = 0.975, na.rm = T)
  
  # min_sigma <- min(mc_sigma)
  min_sigma <- 0
  max_sigma <- max(mc_sigma)

  dat_med <- data.frame(row = rep(1:125, each = 125), col = rep(1:125, times = 125), sol = as.vector(sigma_med))%>% 
    as.data.frame()
  
  plot_heat_sigma_med <- dat_med %>%  ggplot(aes(x = col, y = row, fill = sol)) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(min_sigma, max_sigma),
                         oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "Value")+
    # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
    theme(text = element_text(size=28),
          legend.text = element_text(size = 28),
          legend.key.size = unit(1.5, "cm"))
  # plot_heat_sigma_med
  
  ggsave(filename = "plot_heat_sigma_med.png",
         path = path_fig,
         plot = plot_heat_sigma_med,
         device = "png",
         width = 50,
         height = 45,
         units = "cm",
         dpi = 300
  )
  
  dat_low <- data.frame(row = rep(1:125, each = 125), col = rep(1:125, times = 125), sol = as.vector(sigma_low))%>% 
    as.data.frame()
  
  plot_heat_sigma_low <- dat_low %>%  ggplot(aes(x = col, y = row, fill = sol)) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(min_sigma, max_sigma),
                         oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "Value")+
    # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
    theme(text = element_text(size=28),
          legend.text = element_text(size = 28),
          legend.key.size = unit(1.5, "cm"))
  # plot_heat_sigma_low
  
  ggsave(filename = "plot_heat_sigma_low.png",
         path = path_fig,
         plot = plot_heat_sigma_low,
         device = "png",
         width = 50,
         height = 45,
         units = "cm",
         dpi = 300
  )
  
  dat_high <- data.frame(row = rep(1:125, each = 125), col = rep(1:125, times = 125), sol = as.vector(sigma_high))%>% 
    as.data.frame()
  
  plot_heat_sigma_high <- dat_high %>%  ggplot(aes(x = col, y = row, fill = sol)) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(min_sigma, max_sigma),
                         oob = scales::squish) +
    labs(x = NULL, y = NULL, fill = "Value")+
    # scale_x_continuous(limits = c(-123.8, -114.2), expand = c(0, 0)) +
    # scale_y_continuous(limits = c(32.15, 42.04), expand = c(0, 0)) +
    theme(text = element_text(size=28),
          legend.text = element_text(size = 28),
          legend.key.size = unit(1.5, "cm"))
  # plot_heat_sigma_high
  
  ggsave(filename = "plot_heat_sigma_high.png",
         path = path_fig,
         plot = plot_heat_sigma_high,
         device = "png",
         width = 50,
         height = 45,
         units = "cm",
         dpi = 300
  )
  
  plot_heat_sigma_panel <-ggarrange(plot_heat_sigma_low, plot_heat_sigma_med, plot_heat_sigma_high,
            ncol = 3, nrow = 1,
            labels = c(paste("2.5% Quantile"),
                       paste("50% Quantile"),
                       paste("97.5% Quantile")
            ),
            font.label = list(size = 28),
            vjust = 1.2,
            hjust = -0.5,
            align = "hv",
            common.legend = T,
            legend = "right"
  )
  
  ggsave(filename = paste("plot_heat_sigma_panel.png", sep = ""),
         path = path_fig,
         plot = plot_heat_sigma_panel,
         device = "png",
         width = 60,
         height = 20,
         units = "cm",
         dpi = 300
  )
  
}


{
  ### Error plot----
  #### specific time, location, all inputs (point)----
  {
    alpha_error <- 0.1
    res_pre <- res_pre_MC
    
    time_num <- 15
    sp_num <- 2
    y_true <- dt_pde_test[[time_num]][,sp_num]
    y_pre <- res_pre[[time_num]][,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    error_width <- (max(y_true) - min(y_true)) / 30 # denominator is for aesthetic 
    y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
      geom_abline(col = "red")
  }
  
  #### specific time, input, all locations (field)----
  {
    time_num <- 20
    sp_num <- c(1:N_sp)
    # input_num <- 1
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    error_width <- (max(y_pre_stat["y_true"]) - min(y_pre_stat["y_true"])) / 30
    plot_error_1 <- y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS prediction")
    
    ggsave(filename = paste("plot_error_1_", as.numeric(Sys.time()), ".png", sep = ""),
           path = path_fig,
           plot = plot_error_1,
           device = "png",
           width = 42,
           height = 35,
           units = "cm",
           dpi = 100
    ) 
  }
  
  # panel
  time_p <- tstamp[c(4, 6, 8)]
  sp_num <- c(1:N_sp)
  plot_error_comp_ls <- list()
  plot_band_ls <- list()
  y_pre_stat_error_comp_ls <- list()
  for (t in 1:length(time_p)) {
    time_num <- time_p[t]
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    y_pre_stat_error_comp_ls[[t]] <- y_pre_stat
  }
  
  for (t in 1:length(time_p)) {
    # scatter plot
    error_width <- (max(y_pre_stat_error_comp_ls[[t]]["y_true"]) - min(y_pre_stat_error_comp_ls[[t]]["y_true"])) / 30
    plot_error_1 <- y_pre_stat_error_comp_ls[[t]] %>% ggplot(aes(x = y_true, y = med)) + 
      # geom_linerange(aes(ymin = lower, ymax = upper)) 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error / 3)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error / 3) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS emulation")
    plot_error_comp_ls[[t]] <- plot_error_1
    
    # error band plot
    # fsize <- 34
    plot_band_predict <- y_pre_stat_error_comp_ls[[t]] %>%
      ggplot(aes(x = y_true, y = med)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, x = y_true), fill = "#A6CEE3", alpha = 1) +
      geom_point(aes(y = med), color = "#1F78B4", size = 0.7, alpha = 1) +
      geom_abline(color = "#E31A1C", linewidth = 1) + 
      labs(x = "PDE solution", y = "FFBS emulation") +
      # theme(
      #   axis.title.x = element_text(size = fsize),
      #   axis.title.y = element_text(size = fsize),
      #   axis.text = element_text(size = fsize)  # Adjust axis tick labels as needed
      # ) + 
      theme_minimal()
    plot_band_ls[[t]] <- plot_band_predict
  }
  
  ggsave(filename = paste("plot_band_compare_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = ggarrange(plot_band_ls[[1]], 
                          plot_band_ls[[2]],
                          plot_band_ls[[3]],
                          ncol = 3, nrow = 1,
                          labels = c(paste("t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1)),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          hjust = -2,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 20,
         units = "cm",
         dpi = 300
  )
  
  ggsave(filename = paste("plot_error_compare_", as.numeric(Sys.time()), ".png", sep = ""),
         path = path_fig,
         plot = ggarrange(plot_error_comp_ls[[1]], 
                          plot_error_comp_ls[[2]],
                          plot_error_comp_ls[[3]],
                          ncol = 3, nrow = 1,
                          labels = c(paste("t =", tstamp[4]-1),
                                     paste("t =", tstamp[6]-1),
                                     paste("t =", tstamp[8]-1)),
                          font.label = list(size = 28),
                          vjust = 1.2,
                          hjust = -2,
                          align = "hv",
                          common.legend = T,
                          legend = "right"
                          
         ),
         device = "png",
         width = 60,
         height = 20,
         units = "cm",
         dpi = 300
  )
  
  
}

if(F){
  # #### specific time, all inputs, locations (field panel)----
  # {
  # time_num <- c(10, 20, 30, 40)
  # for (t in 1:length(time_num)) {
  #   if(t == 1){
  #     plot_ls <- list()
  #     y_pre_stat_full_ls <- list()
  #   }
  #   y_pre_stat_full <- data.frame()
  #   for (i in 1:N_sp) {
  #     y_true <- dt_pde_test[[time_num[t]]][,i]
  #     y_pre <- res_pre[[time_num[t]]][,i,]
  #     y_pre_stat <- data.frame(y_true = y_true,
  #                              med = apply(X = y_pre, MARGIN = 1, FUN = median),
  #                              lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
  #                              upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
  #     y_pre_stat_full <- rbind(y_pre_stat_full, y_pre_stat)
  #   }
  #   
  #   y_pre_stat_full <- y_pre_stat_full / N_people
  #   y_pre_stat_full_ls[[t]] <- y_pre_stat_full
  #   error_width <- (max(y_pre_stat_full["y_true"]) - min(y_pre_stat_full["y_true"])) / 30
  #   y_true_max <- max(y_pre_stat_full["y_true"])
  #   
  #   p <- y_pre_stat_full %>% 
  #     dplyr::filter(lower >= -1, upper <= 1) %>% 
  #     ggplot(aes(x = y_true, y = med)) + 
  #     # geom_smooth() + 
  #     geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = 0.05)+
  #     geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = 0.05) + 
  #     geom_abline(col = "red") + 
  #     # lims(y = c(-y_true_max + 0.1, y_true_max + 0.1)) +
  #     labs(x = "PDE solution", y = "FFBS prediction")
  #   plot_ls[[t]] <- p
  # }
  # 
  # ggsave(filename = paste("plot_panel_errorbar4_", as.numeric(Sys.time()), ".png", sep = ""),
  #        path = path_fig,
  #        plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
  #                         plot_ls[[4]], 
  #                         ncol = 2, nrow = 2,
  #                         labels = c(paste("t =", time_num[1]), 
  #                                    paste("t =", time_num[2]),
  #                                    paste("t =", time_num[3]),
  #                                    paste("t =", time_num[4])),
  #                         font.label = list(size = 28),
  #                         vjust = 1.2,
  #                         # hjust = -1,
  #                         align = "hv",
  #                         common.legend = T,
  #                         legend = "right"
  #                         
  #        ),
  #        device = "png",
  #        width = 42,
  #        height = 35,
  #        units = "cm",
  #        dpi = 100
  # ) 
  # }
  # 
  # 
  # #### specific input, all times, locations (all field)----
  # # Idea: if the new input is closer to the training input, the prediction may look better
  # # list: input; col: (space)_time_1-(space)_time_nT; row: replication
  # res_reformback_all <- list()
  # for (i in 1:n_input_new) {
  #   res_reformback <- c()
  #   for (j in 1:nT_ori) {
  #     res_reformback <- rbind(res_reformback, res_pre[[j]][i,,]) # each row includes replications for a location
  #   }
  #   # write.csv(res_reformback, file = paste(path_data, "/res_reformback_", i, ".csv", sep = ""))
  #   res_reformback_all[[i]] <- res_reformback
  # }
  # 
  # # choose input not far away
  # input_full <- rbind(input, input_new)
  # dist_input_full <- as.matrix(stats::dist(input_full, method = "euclidean", diag = T, upper = T))
  # dist_to_train <- t(dist_input_full[1:N,(N+1):dim(dist_input_full)[2]])
  # dist_to_train_mean <- rowMeans(dist_to_train)
  # order_input <- order(dist_to_train_mean)
  # sort(dist_to_train_mean)
  # order_select <- order_input[1:9]
  # 
  # for (i in 1:length(order_select)) {
  #   print(i)
  #   if(i == 1){
  #     plot_ls <- list()
  #     plot_ls_pt <- list()
  #   }
  #   y_true_mat <- dat_pde[[(n_train + order_select[i])]]
  #   y_true <- matrix(t(y_true_mat), ncol = 1)
  #   dt_pre <- cal_errorbar(X = res_reformback_all[[order_select[i]]])
  #   dt_pre <- cbind(y_true, dt_pre) / N_people
  #   
  #   error_width <- (max(dt_pre["y_true"]) - min(dt_pre["y_true"])) / 50
  #   p <- dt_pre %>% 
  #     dplyr::filter(lower >= -1, upper <= 1) %>% 
  #     ggplot(aes(x = y_true, y = med)) + 
  #     geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = 0.01)+
  #     geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = 0.01) + 
  #     geom_abline(color = "red") +
  #     labs(x = "PDE solution", y = "FFBS prediction")
  #     # lims(x = c(0, 1), y = c(-1, 1))
  #   ggsave(filename = paste("plot_error_all", as.numeric(Sys.time()), ".png", sep = ""),
  #          path = path_fig,
  #          plot = p,
  #          device = "png",
  #          width = 42,
  #          height = 35,
  #          units = "cm",
  #          dpi = 100
  #   ) 
  #   
  #   
  #   plot_ls[[i]] <- p
  #   
  #   p2 <- dt_pre %>% 
  #     dplyr::filter(lower >= -1, upper <= 1) %>% 
  #     ggplot(aes(x = y_true, y = med)) + 
  #     geom_point(aes(alpha = 0.01)) + 
  #     geom_smooth(se = F) +
  #     # geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
  #     # geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
  #     geom_abline(color = "red") +
  #     labs(x = "PDE solution", y = "FFBS prediction")
  #   # lims(x = c(0, 1), y = c(-1, 1))
  #   plot_ls_pt[[i]] <- p2
  # }
  # 
  # ggsave(filename = paste("plot_panel_errorbar9_", as.numeric(Sys.time()), ".png", sep = ""),
  #        path = path_fig,
  #        plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
  #                         plot_ls[[4]], plot_ls[[5]], plot_ls[[6]], 
  #                         plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
  #                         ncol = 3, nrow = 3,
  #                         # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
  #                         #            "(4) MBI", "(5) MBA", "(6) Graphical",
  #                         #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
  #                         font.label = list(size = 28),
  #                         vjust = 1.2,
  #                         # hjust = -1,
  #                         align = "hv",
  #                         common.legend = T,
  #                         legend = "right"
  #                         
  #        ),
  #        device = "png",
  #        width = 60,
  #        height = 50,
  #        units = "cm",
  #        dpi = 100
  # ) 
  # 
  # ggsave(filename = paste("plot_panel_pt9_", as.numeric(Sys.time()), ".png", sep = ""),
  #        path = path_fig,
  #        plot = ggarrange(plot_ls_pt[[1]], plot_ls_pt[[2]], plot_ls_pt[[3]],
  #                         plot_ls_pt[[4]], plot_ls_pt[[5]], plot_ls_pt[[6]], 
  #                         plot_ls_pt[[7]], plot_ls_pt[[8]], plot_ls_pt[[9]],
  #                         ncol = 3, nrow = 3,
  #                         # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
  #                         #            "(4) MBI", "(5) MBA", "(6) Graphical",
  #                         #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
  #                         font.label = list(size = 28),
  #                         vjust = 1.2,
  #                         # hjust = -1,
  #                         align = "hv",
  #                         common.legend = T,
  #                         legend = "right"
  #                         
  #        ),
  #        device = "png",
  #        width = 60,
  #        height = 50,
  #        units = "cm",
  #        dpi = 100
  # ) 
}


# # quick check for the plot
# quick_save("plot_eb1", path_fig = path_fig, plot = plot_ls[[1]])
# quick_save("plot_eb2", path_fig = path_fig, plot = plot_ls[[2]])
# quick_save("plot_pt1", path_fig = path_fig, plot = plot_ls_pt[[1]])
# quick_save("plot_pt2", path_fig = path_fig, plot = plot_ls_pt[[2]])
