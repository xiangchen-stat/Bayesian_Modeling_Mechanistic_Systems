This folder reproduces the emulation and calibration experiments using the
installed spDBL package and the repository's saved simulator CSVs. HTML reports,
figures, and results are written here, under vignettes/reproduce_spDBL.

From the repository root:

```sh
# Install the current source into the optional project library:
mkdir -p spDBL/.r-library
R CMD INSTALL --library=spDBL/.r-library spDBL

# Inspect paths and settings without fitting:
Rscript vignettes/reproduce_spDBL/run_reproduction.R all full --dry-run

# Use completed results and regenerate reports; missing/stale results cause an error:
Rscript vignettes/reproduce_spDBL/run_reproduction.R emulation --render-only
Rscript vignettes/reproduce_spDBL/run_reproduction.R all quick --render-only
Rscript vignettes/reproduce_spDBL/run_reproduction.R calibration full --render-only

# Fit when needed and then render:
Rscript vignettes/reproduce_spDBL/run_reproduction.R all quick
Rscript vignettes/reproduce_spDBL/run_reproduction.R calibration full
```

The runner works from any working directory when invoked using its path.
It defaults to all full. Full calibration uses 30,000 emulator draws; quick
uses 300. Both use one chain, one transition per draw, zero warmup, and one
retained sample per draw, with alpha1 and alpha3 fixed at zero. One-transition
quantiles do not establish posterior convergence.

Files needed
- run_reproduction.R: command-line runner.
- reproduction_helpers.R: root discovery, data import, model fitting and caches.
- reproduction_plots.R: figures, including plot_heat_y_panel.
- emulation.Rmd and calibration.Rmd: report sources.
- Installed current spDBL and its dependencies; knitr, rmarkdown, and Pandoc.
  Fresh calibration also needs rstan and a C++ toolchain on first compilation.
- For emulation: data/emulation/pde_para.csv and pde_solution_1.csv through
  pde_solution_50.csv.
- For calibration: the equivalent 51 files in data/calibration/.

Root discovery looks for these data folders above this directory. It does
not require the old package vignette folder, the development reference
scripts, or a source-package DESCRIPTION file. Reference scripts and
R/func_plot.R are included in figure provenance when present. An existing
spDBL/.r-library is preferred; otherwise the ordinary installed package is used.

Reports and completed results
- emulation.html and _results/emulation/
- calibration.html and _results/calibration/ for quick
- calibration-full.html and _results/calibration-full/ for full

Completed runs have been brought over from the previous location, with their
data checksums, settings, package version, and numerical calculations verified.
Their import-provenance.rds files record their origins. The original results
remain unchanged. Local result files are independent copies (APFS clones when
available), so future writes here do not update the original files.

To rerender a completed experiment, keep manifest.rds and results.rds in its
result directory. Calibration also needs calibration.rds for the GP heat panel.
The small plot_heat_y_panel-data.rds cache avoids recalculating its GP means
when the fit checksum and parameter quantile vectors match. Raw CSV files are
still verified against the manifest. Cached results are invalidated when the
input data, numerical helper code, settings, or package version changes.

Direct R Markdown rendering is also supported. calibration.Rmd defaults to
quick; use the runner above for the full report, or explicitly supply
params = list(profile = "full", render_only = TRUE) and
output_file = "calibration-full.html". emulation.Rmd accepts render_only too.

Experiment setup and figures
Emulation retains all 50 saved simulator runs, the 40/10 split, the 100x100
grid, 51 times, ten coefficient draws, and all three covariance models.
Calibration retains the 6x6 grid, 26 times, and original covariance scales.
Synthetic-field generation follows coefficient sampling in the original RNG
stream, so full and quick use different realizations from the same generating
model and saved simulator truth.

Figures retain the original colours, layouts, labels, fonts, and export
settings. plot_heat_y_panel uses GP means at the componentwise parameter
quantile vectors in place of PDE evaluations. Its lower/upper rows are not
pointwise field credible bounds.

This folder lives outside the spDBL package source and is not included by
R CMD build spDBL. Generated results, archives, logs, and HTML are ignored by git.
