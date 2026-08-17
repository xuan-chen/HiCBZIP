# Toy HiCBZIP tutorial.
# Run from the code_clean folder:
# source("examples/run_toy_all_HiCBZIP_variants.R")

source("_common/project_paths.R")
require_packages(c("dplyr", "magrittr", "ggplot2", "gridExtra"))
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(gridExtra)
})

# Load the same core function files used by the analysis scripts.
source_hicbzip_core()

# 1. Load toy diagonal data: 20 bins x 8 cells, 210 x 8 long matrix.
toy_file <- path_here("data", "toy_chr5_diag20_cells8_lambda0.2.rds")
if (!file.exists(toy_file)) {
  stop("Toy data not found: ", toy_file, call. = FALSE)
}
toy <- readRDS(toy_file)
sim_y <- toy$sim_y
coverage <- toy$metadata$lambda

# Display the first cell's observed 20 x 20 contact matrix.
print(round(matrix_long_to_matrix2D(sim_y[, 1]), 1))

# 2. Run all three HiCBZIP variants.
muS_gbnb <- run_BZIP_GB_NB(
  sim_y = sim_y,
  coverage = coverage,
  include_diag = TRUE,
  return_hyper = TRUE
)$muS

muS_ngs <- run_HiCBZIP_NGS(
  sim_y = sim_y,
  coverage = coverage,
  stan_file = path_here("HiCBZIP", "BHZIP_GM_nocov.stan"),
  include_diag = TRUE,
  silent_stan = TRUE
)

muS_nm <- run_HiCBZIP_NM(
  sim_y = sim_y,
  coverage = coverage,
  stan_file = path_here("HiCBZIP", "BHZIP_match_normal.stan"),
  include_diag = TRUE,
  silent_stan = TRUE
)

# 3. Save one output file and one 8-cell heatmap per method.
out_dir <- path_here("examples")
ensure_dir(out_dir)

settings <- list(
  coverage = coverage,
  stan_mcmc = "toy defaults: chains=2, warmup=100, sampling=100",
  silent_stan = TRUE,
  include_diag = TRUE
)

save_cells8_heatmap <- function(mat_long, file_name, title) {
  png(file.path(out_dir, file_name), width = 1800, height = 900, res = 160)
  mat_long_to_heatmaps(
    mat_long,
    n_plot = ncol(mat_long),
    ncol = 4,
    off_diag = FALSE,
    title = title,
    normalize = function(x) log1p(x),
    tick_by = 5
  )
  dev.off()
}

save_cells8_heatmap(
  sim_y,
  "toy_raw_observed_cells8_heatmaps.png",
  "Toy chr5 raw observed counts: 8 cells"
)

saveRDS(
  list(method = "HiCBZIP-GB(NB)", muS = muS_gbnb, metadata = toy$metadata, settings = settings),
  file.path(out_dir, "toy_HiCBZIP_GB_NB_output.rds")
)
save_cells8_heatmap(
  muS_gbnb,
  "toy_HiCBZIP_GB_NB_cells8_heatmaps.png",
  "Toy chr5 HiCBZIP-GB(NB) imputed maps: 8 cells"
)

saveRDS(
  list(method = "HiCBZIP-N(GS)", muS = muS_ngs, metadata = toy$metadata, settings = settings),
  file.path(out_dir, "toy_HiCBZIP_NGS_output.rds")
)
save_cells8_heatmap(
  muS_ngs,
  "toy_HiCBZIP_NGS_cells8_heatmaps.png",
  "Toy chr5 HiCBZIP-N(GS) imputed maps: 8 cells"
)

saveRDS(
  list(method = "HiCBZIP-N(M)", muS = muS_nm, metadata = toy$metadata, settings = settings),
  file.path(out_dir, "toy_HiCBZIP_NM_output.rds")
)
save_cells8_heatmap(
  muS_nm,
  "toy_HiCBZIP_NM_cells8_heatmaps.png",
  "Toy chr5 HiCBZIP-N(M) imputed maps: 8 cells"
)
