# Shared by the two research vignettes. All model fitting uses installed spDBL.
.reproduction_directory <- local({
  files <- vapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L &&
        basename(path) == "reproduction_helpers.R") path else NA_character_
  }, character(1))
  files <- files[!is.na(files)]
  if (!length(files)) stop("Source reproduction_helpers.R from its file")
  dirname(normalizePath(utils::tail(files, 1L), mustWork = TRUE))
})

reproduction_directory <- function() .reproduction_directory

reproduction_results_directory <- function(experiment, profile = "quick") {
  name <- if (experiment == "calibration" && profile == "full") "calibration-full" else experiment
  file.path(reproduction_directory(), "_results", name)
}

reproduction_root <- function(start = reproduction_directory()) {
  start <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(start, "data/emulation/pde_para.csv")) ||
        file.exists(file.path(start, "data/calibration/pde_para.csv"))) return(start)
    parent <- dirname(start)
    if (identical(parent, start)) stop("Cannot locate data/emulation or data/calibration above the reproduction folder")
    start <- parent
  }
}

reproduction_setup <- function(root = reproduction_root()) {
  local_library <- file.path(root, "spDBL/.r-library")
  if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
  if (!requireNamespace("spDBL", quietly = TRUE) ||
      utils::packageVersion("spDBL") < "1.1.1")
    stop("Install current spDBL into spDBL/.r-library before running")
  library(spDBL)
  root
}

reproduction_manifest <- function(root, experiment, settings) {
  paths <- c(file.path(root, "data", experiment, "pde_para.csv"),
    file.path(root, "data", experiment, paste0("pde_solution_", 1:50, ".csv")))
  if (!all(file.exists(paths)))
    stop("Missing simulator input files: ", paste(paths[!file.exists(paths)], collapse = ", "))
  list(schema_version = 2L, settings = settings,
    package_version = as.character(utils::packageVersion("spDBL")),
    helper_md5 = unname(tools::md5sum(file.path(reproduction_directory(), "reproduction_helpers.R"))),
    input_md5 = setNames(unname(tools::md5sum(paths)), substring(paths, nchar(root) + 2)))
}

reproduction_reference_sources <- function(root, experiment) {
  paths <- c(file.path(reproduction_directory(), "reproduction_plots.R"),
    file.path(root, "R/func_plot.R"), file.path(root, "vignettes/reproduce",
      if (experiment == "emulation") "emulate_EP_MNIW.R" else "calibrate_EP_sigma2R.R"))
  paths[file.exists(paths)]
}

reproduction_cached_results <- function(root, experiment, profile = "quick") {
  out <- reproduction_results_directory(experiment, profile)
  files <- file.path(out, c("manifest.rds", "results.rds"))
  if (!all(file.exists(files))) stop("No completed ", experiment, " results in ", out)
  saved <- readRDS(files[1])
  expected <- reproduction_manifest(root, experiment, saved$settings)
  if (experiment == "calibration" && !identical(saved$settings, calibration_reproduction_settings(profile)))
    stop("Saved calibration settings do not match profile ", profile)
  if (!identical(saved, expected)) stop("Saved results do not match the current data, package, or helper code: ", out)
  readRDS(files[2])
}

read_reproduction_data <- function(root, experiment) {
  Nx <- Ny <- if (experiment == "emulation") 100L else 6L
  T <- if (experiment == "emulation") 51L else 26L
  directory <- file.path(root, "data", experiment)
  input <- as.matrix(utils::read.csv(file.path(directory, "pde_para.csv")))
  stopifnot(identical(dim(input), c(50L, 5L)))
  train <- lapply(seq_len(T), function(t) matrix(NA_real_, 40L, Nx * Ny))
  test <- lapply(seq_len(T), function(t) matrix(NA_real_, 10L, Nx * Ny))
  for (i in 1:50) {
    y <- as.matrix(utils::read.csv(file.path(directory, paste0("pde_solution_", i, ".csv"))))
    stopifnot(identical(dim(y), c(T, Nx * Ny)), all(is.finite(y)))
    for (t in seq_len(T)) {
      if (i <= 40L) train[[t]][i, ] <- y[t, ] else test[[t]][i - 40L, ] <- y[t, ]
    }
    if (i %% 10L == 0L) message("Read ", experiment, " simulator run ", i, "/50")
  }
  list(input = input, train = train, test = test, Nx = Nx, Ny = Ny, T = T,
    grid = data.frame(row = rep(seq_len(Ny), Nx), col = rep(seq_len(Nx), each = Ny)))
}

episode_setup <- function(d, window) {
  C <- d$Ny * window
  B <- d$Nx / window
  index <- spDBL::generate.grid.rowsnake(40L, d$Nx * d$Ny, 40L, C)
  index_test <- spDBL::generate.grid.rowsnake(10L, d$Nx * d$Ny, 10L, C)
  # The original scripts build both AR designs before drawing coefficients.
  set.seed(1234L)
  F <- spDBL::gen_F_ls_AR2_EP(d$train, d$T, B, index)
  F_test <- spDBL::gen_F_ls_AR2_EP(d$test, d$T, B, index_test)
  block_cols <- t(vapply(seq_len(B), function(b) seq.int(index[b, 4], index[b, 5]), integer(C)))
  Y <- unlist(lapply(d$train, function(y) lapply(seq_len(B), function(b)
    y[, block_cols[b, ], drop = FALSE])), recursive = FALSE)
  phi <- 3 / (0.5 * max(stats::dist(d$input[1:40, ])))
  phi_sp <- 3 / (0.5 * max(stats::dist(d$grid[seq_len(C), ])))
  list(F = F, F_test = F_test, Y = Y, C = C, B = B, E = d$T * B,
    P = ncol(F[[1]]), block_cols = block_cols, phi = phi, phi_sp = phi_sp,
    V = spDBL::gen_exp_kernel(d$input[1:40, ], phi, 1.1, 1e-4),
    R = spDBL::gen_gp_kernel(d$grid[seq_len(C), ], phi_sp, 1.1, 1e-4),
    rng_before_sampling = .Random.seed)
}

prediction_scores <- function(prediction, truth, method) {
  do.call(rbind, lapply(seq_len(nrow(truth[[1]])), function(i) {
    actual <- unlist(lapply(truth, function(y) y[i, ]), use.names = FALSE)
    predicted <- unlist(lapply(prediction, function(y) y[i, ]), use.names = FALSE)
    data.frame(method = method, input_id = 40L + i,
      RMSE = sqrt(mean((predicted - actual)^2)), MAE = mean(abs(predicted - actual)),
      max_abs_error = max(abs(predicted - actual)))
  }))
}

log_mean_exp <- function(x) { a <- max(x); a + log(mean(exp(x - a))) }

emulation_model <- function(model, d, s, out) {
  file <- file.path(out, paste0(model, ".rds"))
  if (file.exists(file)) { message("Reuse completed ", model); return(readRDS(file)) }
  args <- list(Y = s$Y, F_ls = s$F, G_ls = diag(s$P), W_ls = diag(s$P),
    V_ls = s$V, m0 = matrix(1, s$P, s$C), M0 = diag(s$P), nT = s$E,
    delta = 1, verbose = TRUE)
  message("Fit ", model, ": ", s$E, " episodes; coefficient dimension ", s$P, " x ", s$C)
  fit_time <- system.time({
    fit <- switch(model,
      MNIW = do.call(spDBL::FFBS, c(args, list(n0 = s$P + 3, D0 = diag(s$C)))),
      MNIG = do.call(spDBL::FFBS_sigma2R, c(args, list(n0 = s$P + 3, D0 = s$P + s$C, R = s$R))),
      I = do.call(spDBL::FFBS_I, args))
  })[["elapsed"]]
  message(model, " FFBS finished in ", round(fit_time, 1), " seconds")
  sampling_time <- system.time({
    draws <- switch(model,
      MNIW = spDBL::FFBS_sampling(10L, fit, s$F, diag(s$P), s$E, sampling = "marginal"),
      MNIG = spDBL::FFBS_sampling_sigma2R(10L, fit, s$F, diag(s$P), s$E, s$R, sampling = "marginal"),
      I = spDBL::FFBS_sampling_I(10L, fit, s$F, diag(s$P), s$E))
  })[["elapsed"]]
  # Check the optimized log-density against the package study's density.
  Sigma <- function(m) switch(model, MNIW = draws$Sigma[, , m],
    MNIG = draws$Sigma[m] * s$R, I = diag(s$C))
  L_V <- t(chol(s$V))
  L_C <- lapply(1:10, function(m) t(chol(Sigma(m))))
  constant <- vapply(L_C, function(L) -40 * s$C / 2 * log(2 * pi) -
    s$C * sum(log(diag(L_V))) - 40 * sum(log(diag(L))), numeric(1))
  density <- function(e, m) {
    q <- forwardsolve(L_V, s$Y[[e]] - s$F[[e]] %*% draws[[e]][, , m])
    q <- forwardsolve(L_C[[m]], t(q))
    constant[m] - sum(q^2) / 2
  }
  density_check <- mniw::dMNorm(s$Y[[1]], s$F[[1]] %*% draws[[1]][, , 1],
    SigmaR = s$V, SigmaC = Sigma(1), log = TRUE)
  stopifnot(isTRUE(all.equal(density(1, 1), as.numeric(density_check), tolerance = 1e-8)))
  log_lik <- matrix(NA_real_, 10L, s$E)
  analytic <- G <- trace_S <- numeric(s$E)
  terminal <- fit$ff[[s$E]]
  score_time <- system.time(for (e in seq_len(s$E)) {
    for (m in 1:10) log_lik[m, e] <- density(e, m)
    st <- fit$bs$st[, , e]
    St <- fit$bs$St[, , e]
    analytic[e] <- switch(model,
      MNIW = spDBL::lppd_IW_1t(s$Y[[e]], s$F[[e]], s$V, st, St, terminal$nt, terminal$Dt),
      MNIG = spDBL::lppd_IG_1t(s$Y[[e]], s$F[[e]], s$V, st, St, terminal$nt, terminal$Dt, s$R),
      I = spDBL::lppd_id_1t(s$Y[[e]], s$F[[e]], s$V, st, St))
    G[e] <- sum((s$Y[[e]] - s$F[[e]] %*% st)^2)
    trace_S[e] <- sum(diag(St))
    if (e %% 510L == 0L) message(model, " scoring ", e, "/", s$E)
  })[["elapsed"]]
  lppd <- apply(log_lik, 2, log_mean_exp)
  pwaic <- apply(log_lik, 2, stats::var)
  P <- switch(model,
    MNIW = sum(diag(terminal$Dt)) / (terminal$nt - 2),
    MNIG = sum(diag(s$R)) * terminal$Dt / (terminal$nt - 1), I = s$C) * sum(trace_S)
  result <- list(model = model,
    summary = data.frame(model = model, lppd = sum(lppd), lppd_analytic = sum(analytic),
      p_waic = sum(pwaic), WAIC = -2 * sum(lppd - pwaic),
      G = sum(G), P = P, D = sum(G) + P,
      FFBS_seconds = fit_time, sampling_seconds = sampling_time, scoring_seconds = score_time),
    pointwise = data.frame(episode = seq_len(s$E), lppd = lppd,
      lppd_analytic = analytic, p_waic = pwaic, WAIC = -2 * (lppd - pwaic)),
    log_lik = log_lik, density_check_absolute_error = abs(density(1, 1) - density_check))
  if (model == "MNIW") {
    message("Predict all ten held-out 100x100 fields")
    episodic <- spDBL::FFBS_predict_exact(s$Y, fit, s$F, s$F_test,
      d$input[1:40, ], d$input[41:50, ], s$E, s$phi)
    result$prediction <- spDBL::recover_from_EP_exact(episodic, d$T, s$E, s$block_cols)
    result$validation <- prediction_scores(result$prediction, d$test, "AR2 with saved held-out lags")
    # Only the three plotted times need the MC arrays; retain every spatial site.
    result$plot_times <- as.integer(seq(1, d$T, length.out = 9))[c(4, 6, 8)]
    selected <- unlist(lapply(result$plot_times, function(t) (t - 1L) * s$B + seq_len(s$B)))
    mc <- spDBL::FFBS_predict_MC(10L, s$Y[selected], draws[selected], s$F[selected],
      s$F_test[selected], d$input[1:40, ], d$input[41:50, ], length(selected), s$phi)
    result$mean_draws <- spDBL::recover_from_EP_MC(mc, 3L, length(selected), 10L)
    result$truth <- lapply(result$plot_times, function(t) d$test[[t]][2, ])
    result$grid <- d$grid
  }
  stopifnot(all(is.finite(as.matrix(result$summary[, -1]))), all(is.finite(log_lik)))
  saveRDS(result, file, compress = "gzip")
  message("Saved ", normalizePath(file))
  result
}

run_emulation_reproduction <- function(root = reproduction_root()) {
  root <- reproduction_setup(root)
  out <- reproduction_results_directory("emulation")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  settings <- list(Nx = 100L, Ny = 100L, times = 0:50, train = 1:40,
    validation = 41:50, window = 1L, nsam = 10L, AR = 2L, seed = 1234L,
    gp_tune = 0.5, gp_sigma2 = 1.1, gp_tau2 = 1e-4, sampling = "marginal")
  manifest <- reproduction_manifest(root, "emulation", settings)
  manifest_file <- file.path(out, "manifest.rds")
  if (file.exists(manifest_file) && !identical(readRDS(manifest_file), manifest))
    stop("Emulation cache is stale; move ", out, " aside before rerunning")
  saveRDS(manifest, manifest_file)
  final_file <- file.path(out, "results.rds")
  if (file.exists(final_file)) return(readRDS(final_file))
  d <- read_reproduction_data(root, "emulation")
  s <- episode_setup(d, 1L)
  models <- list()
  for (model in c("MNIW", "MNIG", "I")) {
    models[[model]] <- emulation_model(model, d, s, out)
    gc()
  }
  results <- list(manifest = manifest, models = models,
    comparison = do.call(rbind, lapply(models, `[[`, "summary")))
  utils::write.csv(results$comparison, file.path(out, "model-comparison.csv"), row.names = FALSE)
  utils::write.csv(models$MNIW$validation, file.path(out, "heldout-scores.csv"), row.names = FALSE)
  saveRDS(results, final_file, compress = "gzip")
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  results
}

calibration_reproduction_settings <- function(profile = c("full", "quick")) {
  profile <- match.arg(profile)
  settings <- list(Nx = 6L, Ny = 6L, times = 0:25, train = 1:40,
    validation = 41:50, window = 2L,
    nsam = if (profile == "full") 30000L else 300L, AR = 2L, seed = 1234L,
    gp_tune = 0.5, gp_sigma2 = 1.1, gp_tau2 = 1e-4, sampling = "marginal",
    lower = c(2, 0.2, 0, 0, 0), upper = c(4, 0.4, 0, 0.2, 0),
    iter = 1L, warmup = 0L, chains = 1L, keep = 1L, calibration_seed = 12345678L,
    n0 = 30, d0 = 2, discount = 1, discrepancy_covariance = "original unnormalised kernel")
  settings
}

run_calibration_reproduction <- function(root = reproduction_root(), profile = c("full", "quick")) {
  root <- reproduction_setup(root)
  profile <- match.arg(profile)
  settings <- calibration_reproduction_settings(profile)
  out <- reproduction_results_directory("calibration", profile)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  manifest <- reproduction_manifest(root, "calibration", settings)
  manifest_file <- file.path(out, "manifest.rds")
  if (file.exists(manifest_file) && !identical(readRDS(manifest_file), manifest))
    stop("Calibration cache is stale; move ", out, " aside before rerunning")
  saveRDS(manifest, manifest_file)
  final_file <- file.path(out, "results.rds")
  if (file.exists(final_file)) return(readRDS(final_file))
  d <- read_reproduction_data(root, "calibration")
  s <- episode_setup(d, 2L)
  # Leave the seed stream after BOTH original AR designs, as in the source.
  em_time <- system.time({
    em <- spDBL::emulator_learn(d$input[1:40, ], d$train, d$Nx, d$Ny,
      F_ls_train = s$F, covariance = "sigma2R", episode_window = 2L,
      nsam = settings$nsam, sampling = "marginal", seed = NULL)
  })[["elapsed"]]
  stopifnot(identical(em$setup$block_cols, s$block_cols))
  # Reconstruct the original stream explicitly: MNIG_sampler draws gamma
  # variates, then exactly E * P * C * M normal variates; GP means consume none.
  # This also checks that a future sampler change cannot silently change z.
  rng_after_sampling <- .Random.seed
  assign(".Random.seed", s$rng_before_sampling, envir = .GlobalEnv)
  terminal <- em$para_ffbs$ff[[s$E]]
  invisible(invgamma::rinvgamma(settings$nsam, shape = terminal$nt, rate = terminal$Dt))
  for (e in seq_len(s$E)) invisible(stats::rnorm(s$P * s$C * settings$nsam))
  stopifnot(identical(rng_after_sampling, .Random.seed))
  phi_field <- 3 / (0.5 * max(stats::dist(d$grid)))
  U <- spDBL::gen_gp_kernel(d$grid, phi_field, 1.1, 1e-4)
  y_true <- do.call(rbind, lapply(d$test, function(y) y[10, ]))
  generated <- spDBL::gen_calibrate_data(y_true,
    list(n0 = 30, d0 = 2, b = 1, u0 = rep(0, 36L)), U)
  sites <- sample(seq_len(36L))
  field <- c(generated, list(y = y_true, eta = d$input[50, ], U = U,
    locations = sites, grid = d$grid, u0 = rep(0, 36L),
    field_rng_before_generation = rng_after_sampling))
  saveRDS(field, file.path(out, "field.rds"))
  cal <- spDBL::calibration_prepare(em, settings$lower, settings$upper)
  direct <- spDBL::emulator_predict(em, d$input[41:50, ], F_new_ls = s$F_test)
  offline <- spDBL::emulator_predict(cal, d$input[41:50, ])
  validation <- rbind(prediction_scores(direct, d$test, "Original saved-lag AR2 means"),
    prediction_scores(offline, d$test, "PDE-free offline AR2 GP means"))
  utils::write.csv(validation, file.path(out, "heldout-scores.csv"), row.names = FALSE)
  # Loading a known compatible cache avoids recompilation, when available.
  shared_cache <- file.path(root, "output/calibration-vignette-marginal-oneiter500/stan-cache")
  cache <- if (dir.exists(shared_cache)) shared_cache else file.path(out, "stan-cache")
  fit_time <- system.time({
    fit <- spDBL::calibration_fit(cal, field$z[, sites], field$U,
      locations = sites, seed = 12345678L, iter = 1L, warmup = 0L,
      chains = 1L, draws_per_emulator = 1L, n0 = 30, d0 = 2,
      discount = 1, refresh = 0L, checkpoint_dir = file.path(out, "conditionals"),
      cache_dir = cache)
  })[["elapsed"]]
  predictions <- spDBL::calibration_predict(fit, seed = 12345679L)
  saveRDS(fit, file.path(out, "calibration.rds"), compress = "gzip")
  saveRDS(predictions, file.path(out, "predictions.rds"), compress = "gzip")
  q <- function(a) apply(a, c(2, 3), stats::quantile, c(0.025, 0.5, 0.975))
  tau_q <- apply(fit$mcmc_tau2, 2, stats::quantile, c(0.025, 0.5, 0.975))
  eta_q <- t(apply(fit$mcmc_eta, 2, stats::quantile, c(0.025, 0.5, 0.975)))
  eta <- data.frame(parameter = colnames(fit$mcmc_eta), true = field$eta,
    lower = eta_q[, 1], median = eta_q[, 2], upper = eta_q[, 3], row.names = NULL)
  quantiles <- lapply(predictions, q)
  truth <- list(y = field$y[, sites], u = field$u[, sites], z_rep = field$z[, sites])
  scores <- do.call(rbind, lapply(names(quantiles), function(name) {
    a <- quantiles[[name]]; actual <- truth[[name]]
    data.frame(quantity = name, RMSE = sqrt(mean((a[2, , ] - actual)^2)),
      interval_fraction = mean(a[1, , ] <= actual & actual <= a[3, , ]),
      mean_interval_width = mean(a[3, , ] - a[1, , ]))
  }))
  results <- list(manifest = manifest, field = field, validation = validation,
    eta = eta, eta_draws = fit$mcmc_eta, tau2 = tau_q, quantiles = quantiles,
    scores = scores, diagnostics = fit$diagnostics,
    settings = fit$settings, execution_only = fit$execution_only,
    active_parameters = cal$parameter_names[cal$predictor$upper > cal$predictor$lower],
    fixed_parameters = stats::setNames(cal$predictor$lower[cal$predictor$upper == cal$predictor$lower],
      cal$parameter_names[cal$predictor$upper == cal$predictor$lower]),
    timing = data.frame(stage = c("emulation", "calibration"), seconds = c(em_time, fit_time)))
  utils::write.csv(eta, file.path(out, "eta-summary.csv"), row.names = FALSE)
  utils::write.csv(scores, file.path(out, "recovery-scores.csv"), row.names = FALSE)
  saveRDS(results, final_file, compress = "gzip")
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  results
}

