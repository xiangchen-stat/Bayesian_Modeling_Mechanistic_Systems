## Load packages and data ----
# Check if pacman package is installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  # If not installed, install pacman
  install.packages("pacman")
}

# using pacman to load other packages
library(pacman)
options(tidyverse.quiet = TRUE)
p_load(ggpubr) # for table and plots
# p_load(MBA, reshape2, ggmap, sf) # for spatial analysis

# Set ggplot theme
theme_set(theme_minimal(base_size = 22))
col_epa <- c("#00e400", "#ffff00", "#ff7e00", "#ff0000", "#99004c", "#7e0023")
col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
             "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429",
             "#a62021", "#871b1c")

## Plot PDE results ----
Nx <- 50
Ny <- 50
n_input <- 1
nT <- 51
ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
sir <- list()
path <- "Dan_pde_nopermute_train_50x50_51_1_1"
path_data <- file.path("..", "data", path)
path_fig <- file.path("..", "figures", path)
for (i in 1:n_input) {
  sir[[i]] <- read.csv(file = paste0(path_data, "/pde_solution_", i, ".csv", sep = ""))
}

plot_ls <- list() 
tstamp <- as.integer(seq(1, nT, length.out = 9))
file_num <- 1
ind_plot <- 1

# set max limit for all plots
max_y <- max(sir[[file_num]])
max_lny <- log(1 + max_y)
# Plot log(1 + y)
for (j in tstamp) {
  if(j == 1){
    ind_plot <- 1 # start counting
  }
  temp <- t(sir[[file_num]][j,])
  rownames(temp) <- NULL
  colnames(temp) <- NULL
  dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol =  temp)%>% 
    as.data.frame()
  
  p <- ggplot(dt, aes(x = col, y = row, fill = log(1+sol))) +
    geom_raster() +
    scale_fill_gradientn(colours = col_bgr,
                         limits = c(0, max_lny),
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


# ggsave(filename = paste0("plot_sir", file_num, "_", j, ".png"),
#        path = path_fig,
#        plot = plot_ls,
#        device = "png",
#        width = 60,
#        height = 50,
#        units = "cm",
#        dpi = 50
# )

ggsave(filename = "plot_panel_sir.png",
       path = path_fig,
       plot = ggarrange(plot_ls[[1]], plot_ls[[2]], plot_ls[[3]],
                        plot_ls[[4]], plot_ls[[5]], plot_ls[[6]],
                        plot_ls[[7]], plot_ls[[8]], plot_ls[[9]],
                        ncol = 3, nrow = 3,
                        # labels = c("(1) Data", "(2) Gaussian Process", "(3) Bilinear", 
                        #            "(4) MBI", "(5) MBA", "(6) Graphical",
                        #            "(7) Inverse Graphical", "(8) Nearest Neighbor"),
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


