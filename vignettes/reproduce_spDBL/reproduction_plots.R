# Figure code is separate from inference cache fingerprints.
plot_emulation_reproduction <- function(r, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  m <- r$models$MNIW
  # Exact palette and global colour limits in R/func_plot.R and the source
  # experiment. The maximum is over ALL ten saved held-out trajectories.
  col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
    "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429", "#a62021", "#871b1c")
  root <- reproduction_root()
  paths <- file.path(root, "data/emulation", paste0("pde_solution_", 41:50, ".csv"))
  max_y <- max(vapply(paths, function(path) max(as.matrix(utils::read.csv(path))), numeric(1)))
  theme <- ggplot2::theme_minimal(base_size = 22)
  time_labels <- paste("t =", m$plot_times - 1L)
  heatmap <- function(values) {
    df <- data.frame(m$grid, sol = values)
    ggplot2::ggplot(df, ggplot2::aes(col, row, fill = sol)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colours = col_bgr, limits = c(0, max_y), oob = scales::squish) +
      ggplot2::labs(x = "x", y = "y", fill = "Value") + theme +
      ggplot2::theme(text = ggplot2::element_text(size = 28),
        legend.text = ggplot2::element_text(size = 28), legend.key.size = grid::unit(1.5, "cm"),
        # Leave room for the source panel labels above the top y-axis tick.
        plot.margin = ggplot2::margin(t = 25, r = 5.5, b = 5.5, l = 5.5))
  }
  panels <- c(lapply(m$truth, heatmap),
    lapply(m$plot_times, function(t) heatmap(m$prediction[[t]][2, ])))
  heat <- ggpubr::ggarrange(plotlist = panels, ncol = 3, nrow = 2,
    labels = c(paste("PDE solution:", time_labels[1]), time_labels[2:3],
      paste("FFBS emulation:", time_labels[1]), time_labels[2:3]),
    font.label = list(size = 28), vjust = 1.2, hjust = -0.5, align = "hv",
    common.legend = TRUE, legend = "right")
  intervals <- lapply(1:3, function(i) {
    a <- m$mean_draws[[i]][2, , ]
    q <- t(apply(a, 1, stats::quantile, c(0.025, 0.5, 0.975)))
    data.frame(actual = m$truth[[i]] / 10000, lower = q[, 1] / 10000,
      median = q[, 2] / 10000, upper = q[, 3] / 10000,
      time = time_labels[i])
  })
  error_panels <- lapply(intervals, function(df) {
    error_width <- diff(range(df$actual)) / 30
    ggplot2::ggplot(df, ggplot2::aes(actual, median)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), size = 0.2, alpha = 0.1 / 3) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = error_width, alpha = 0.1 / 3) +
      ggplot2::geom_abline(colour = "red") +
      ggplot2::labs(x = "PDE solution", y = "FFBS emulation") + theme
  })
  error <- ggpubr::ggarrange(plotlist = error_panels, ncol = 3, nrow = 1,
    labels = time_labels, font.label = list(size = 28), vjust = 1.2, hjust = -2,
    align = "hv", common.legend = TRUE, legend = "right")
  comparison <- do.call(rbind, lapply(r$models, function(m)
    data.frame(m$pointwise, model = m$model)))
  compare <- ggpubr::ggarrange(
    ggplot2::ggplot(comparison, ggplot2::aes(episode, WAIC, colour = model)) +
      ggplot2::geom_line() + theme + ggplot2::labs(x = "episode", y = "waic"),
    ggplot2::ggplot(comparison, ggplot2::aes(episode, lppd, colour = model)) +
      ggplot2::geom_line() + theme + ggplot2::labs(x = "episode"),
    ncol = 2, nrow = 1, labels = c("WAIC", "lppd"), font.label = list(size = 28),
    vjust = 1.2, hjust = -4, align = "hv", common.legend = TRUE, legend = "right")
  plots <- list(heat_compare = heat, error_compare = error, model_comparison = compare)
  exports <- data.frame(name = names(plots), width_cm = c(60, 60, 40),
    height_cm = c(35, 20, 20), dpi = c(100, 300, 100))
  for (i in seq_len(nrow(exports))) ggplot2::ggsave(
    file.path(directory, paste0(exports$name[i], ".png")), plots[[exports$name[i]]],
    width = exports$width_cm[i], height = exports$height_cm[i], units = "cm",
    dpi = exports$dpi[i], bg = "white")
  saveRDS(list(palette = col_bgr, fill_limits = c(0, max_y), input_id = 42L,
    times = m$plot_times - 1L, legend = "right", base_size = 22,
    error_alpha = 0.1 / 3, errorbar_width_divisor = 30, exports = exports,
    source_md5 = tools::md5sum(reproduction_reference_sources(root, "emulation"))),
    file.path(directory, "plot-settings.rds"))
  invisible(plots)
}

plot_calibration_reproduction <- function(r, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  root <- reproduction_root()
  theme <- ggplot2::theme_minimal(base_size = 22)
  active <- which(r$manifest$settings$upper > r$manifest$settings$lower)
  colours <- c("orange", "lightgreen", "purple", "lightyellow", "orange")
  eta_panels <- lapply(active, function(j) {
    df <- data.frame(value = r$eta_draws[, j])
    limits <- c(r$manifest$settings$lower[j], r$manifest$settings$upper[j]) +
      c(-1, 1) * if (j == 1L) 0.1 else 0.01
    ggplot2::ggplot(df, ggplot2::aes(value)) +
      ggplot2::geom_histogram(bins = 30, fill = colours[j], colour = "black") +
      ggplot2::geom_vline(xintercept = r$eta$true[j], linewidth = 1, colour = "red") +
      ggplot2::geom_vline(xintercept = r$eta$lower[j], linetype = "dashed", colour = "black") +
      ggplot2::geom_vline(xintercept = r$eta$median[j], linetype = "solid", colour = "black") +
      ggplot2::geom_vline(xintercept = r$eta$upper[j], linetype = "dashed", colour = "black") +
      ggplot2::lims(x = limits) + ggplot2::labs(x = NULL, y = NULL, title = paste("eta", j)) + theme
  })
  eta_plot <- ggpubr::ggarrange(plotlist = eta_panels, ncol = length(active), nrow = 1,
    font.label = list(size = 28), vjust = 1, hjust = -0.3, align = "hv",
    common.legend = TRUE, legend = "right")
  tau <- data.frame(time = 0:25, true = r$field$tau2,
    lower = r$tau2[1, ], median = r$tau2[2, ], upper = r$tau2[3, ])
  tau_plot <- ggplot2::ggplot(tau, ggplot2::aes(time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#A6CEE3", alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = median), colour = "#1F78B4", linewidth = 1, linetype = "dotdash") +
    ggplot2::geom_line(ggplot2::aes(y = true), colour = "#E31A1C", linewidth = 1) +
    ggplot2::labs(x = "Time", y = "Value", title = expression("Variance (" * tau^2 * ")")) +
    ggplot2::ylim(0, max(tau$upper)) + theme
  recovery <- function(name) {
    actual <- if (name == "u") r$field$u else r$field$z
    a <- r$quantiles[[name]]
    df <- data.frame(actual = as.vector(actual[, r$field$locations]),
      lower = as.vector(a[1, , ]), median = as.vector(a[2, , ]), upper = as.vector(a[3, , ]))
    df <- df[order(df$actual), ]
    ggplot2::ggplot(df, ggplot2::aes(actual, median)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#A6CEE3", alpha = 1) +
      ggplot2::geom_point(colour = "#1F78B4", size = 0.7, alpha = 1) +
      ggplot2::geom_abline(colour = "#E31A1C", linewidth = 1) +
      ggplot2::labs(x = if (name == "u") "True spatial bias" else "True field value",
        y = "Posterior prediction", title = if (name == "u")
          expression("Spatial Bias (" * u * ")") else expression("Field Value (" * z * ")")) + theme
  }
  panel <- ggpubr::ggarrange(tau_plot, recovery("u"), recovery("z_rep"), ncol = 3, nrow = 1,
    font.label = list(size = 24), align = "hv", common.legend = TRUE, legend = "right")

  # Match the original parameter-quantile panel: its rows evaluate parameter
  # vectors, rather than pointwise quantiles of recovered y. Use GP means in
  # place of its PDE calls, retaining the saved simulator truth in row one.
  fit_file <- file.path(dirname(directory), "calibration.rds")
  if (!file.exists(fit_file)) stop("The GP heat panel requires the saved fit: ", fit_file)
  fit_md5 <- unname(tools::md5sum(fit_file))
  eta_vectors <- rbind(lower = r$eta$lower, upper = r$eta$upper)
  colnames(eta_vectors) <- r$eta$parameter
  heat_cache_file <- file.path(directory, "plot_heat_y_panel-data.rds")
  cached <- if (file.exists(heat_cache_file)) tryCatch(readRDS(heat_cache_file),
    error = function(e) NULL) else NULL
  if (!is.null(cached) && identical(cached$fit_md5, fit_md5) &&
      identical(cached$eta, eta_vectors) && identical(cached$fields$truth, r$field$y)) {
    heat_fields <- cached$fields
  } else {
    fit <- readRDS(fit_file)
    emulated <- spDBL::emulator_predict(fit$emulator, eta_vectors)
    heat_fields <- list(truth = r$field$y,
      lower = do.call(rbind, lapply(emulated, function(y) y[1, ])),
      upper = do.call(rbind, lapply(emulated, function(y) y[2, ])))
  }
  stopifnot(all(vapply(heat_fields, function(y)
    identical(dim(y), dim(r$field$y)) && all(is.finite(y)), logical(1))))
  tstamp <- as.integer(seq(1, nrow(r$field$y), length.out = 9))
  times <- tstamp[c(2, 5, 8)]
  col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
    "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429", "#a62021", "#871b1c")
  paths <- file.path(root, "data/calibration", paste0("pde_solution_", 41:50, ".csv"))
  max_y <- max(vapply(paths, function(path) max(as.matrix(utils::read.csv(path))), numeric(1)))
  heatmap <- function(y, t) {
    df <- data.frame(r$field$grid, sol = y[t, ])
    ggplot2::ggplot(df, ggplot2::aes(col, row, fill = sol)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_gradientn(colours = col_bgr, limits = c(0, max_y), oob = scales::squish) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Value") + theme +
      ggplot2::theme(text = ggplot2::element_text(size = 28),
        legend.text = ggplot2::element_text(size = 28), legend.key.size = grid::unit(1.5, "cm"),
        plot.margin = ggplot2::margin(t = 25, r = 5.5, b = 5.5, l = 5.5))
  }
  heat_panels <- unlist(lapply(heat_fields, function(y) lapply(times, function(t)
    heatmap(y, t))), recursive = FALSE)
  labels <- unlist(lapply(c("True eta:", "2.5% eta:", "97.5% eta:"), function(label)
    c(paste(label, "t =", times[1] - 1L), paste("t =", times[2:3] - 1L))))
  heat_panel <- ggpubr::ggarrange(plotlist = heat_panels, ncol = 3, nrow = 3,
    labels = labels, font.label = list(size = 28), vjust = 1.2, hjust = -0.3,
    align = "hv", common.legend = TRUE, legend = "right")
  heat_panel <- ggpubr::annotate_figure(heat_panel,
    left = ggpubr::text_grob("Y", rot = 90, vjust = 1, size = 28),
    bottom = ggpubr::text_grob("X", vjust = 0, size = 28))
  plots <- list(eta_histograms = eta_plot, tau2_u_z_panel = panel, plot_heat_y_panel = heat_panel)
  exports <- data.frame(name = names(plots), width_cm = c(45, 45, 40),
    height_cm = c(15, 15, 35), dpi = 300)
  for (i in seq_len(nrow(exports))) ggplot2::ggsave(
    file.path(directory, paste0(exports$name[i], ".png")), plots[[exports$name[i]]],
    width = exports$width_cm[i], height = exports$height_cm[i], units = "cm",
    dpi = exports$dpi[i], bg = "white")
  invisible(file.copy(file.path(directory, "eta_histograms.png"),
    file.path(directory, "plot_eta_panel.png"), overwrite = TRUE))
  invisible(file.copy(file.path(directory, "tau2_u_z_panel.png"),
    file.path(directory, "plot_panel_tau_u_z.png"), overwrite = TRUE))
  saveRDS(list(fit_md5 = fit_md5, eta = eta_vectors, fields = heat_fields, time_indices = times,
    times = times - 1L, fill_limits = c(0, max_y),
    interpretation = "Saved simulator truth; GP posterior means at componentwise parameter quantile vectors, not pointwise field credible bounds."),
    file.path(directory, "plot_heat_y_panel-data.rds"))
  saveRDS(list(histogram_fills = colours[active], histogram_bins = 30,
    histogram_limits = lapply(active, function(j)
      c(r$manifest$settings$lower[j], r$manifest$settings$upper[j]) + c(-1, 1) * if (j == 1L) 0.1 else 0.01),
    ribbon = "#A6CEE3", estimate = "#1F78B4", truth = "#E31A1C",
    palette = col_bgr, base_size = 22, exports = exports,
    source_md5 = tools::md5sum(reproduction_reference_sources(root, "calibration"))),
    file.path(directory, "plot-settings.rds"))
  invisible(plots)
}
