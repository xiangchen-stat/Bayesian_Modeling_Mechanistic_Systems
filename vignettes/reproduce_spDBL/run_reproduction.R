#!/usr/bin/env Rscript
# Rscript vignettes/reproduce_spDBL/run_reproduction.R [emulation|calibration|all] [full|quick] [--dry-run|--render-only]
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[startsWith(args, "--file=")][1])
directory <- dirname(normalizePath(this_file, mustWork = TRUE))
source(file.path(directory, "reproduction_helpers.R"))
root <- reproduction_setup(reproduction_root(directory))
args <- commandArgs(trailingOnly = TRUE)
flags <- args[startsWith(args, "--")]
positional <- args[!startsWith(args, "--")]
if (length(positional) > 2L) stop("Too many positional arguments")
mode <- if (length(positional)) match.arg(positional[1], c("emulation", "calibration", "all")) else "all"
profile <- if (length(positional) >= 2L) match.arg(positional[2], c("full", "quick")) else "full"
if (any(!flags %in% c("--dry-run", "--render-only")) || anyDuplicated(flags)) stop("Unknown arguments")
dry_run <- "--dry-run" %in% flags
render_only <- "--render-only" %in% flags
if (dry_run) {
  if (mode != "emulation") {
    print(calibration_reproduction_settings(profile))
    cat("Calibration output:", reproduction_results_directory("calibration", profile), "\n")
    cat("Report:", if (profile == "full") "calibration-full.html" else "calibration.html", "\n")
  }
  if (mode != "calibration") cat("Emulation: 100 x 100 grid, 51 times, 10 draws; output:",
    reproduction_results_directory("emulation"), "\n")
  quit(save = "no", status = 0L)
}
for (experiment in if (mode == "all") c("emulation", "calibration") else mode) {
  if (render_only) invisible(reproduction_cached_results(root, experiment, profile)) else
    if (experiment == "emulation") run_emulation_reproduction(root) else run_calibration_reproduction(root, profile)
  input <- file.path(directory, paste0(experiment, ".Rmd"))
  if (experiment == "calibration") {
    rmarkdown::render(input, quiet = TRUE, params = list(profile = profile, render_only = render_only),
      output_file = if (profile == "full") "calibration-full.html" else "calibration.html",
      envir = new.env(parent = globalenv()))
  } else rmarkdown::render(input, quiet = TRUE, params = list(render_only = render_only),
    envir = new.env(parent = globalenv()))
}
