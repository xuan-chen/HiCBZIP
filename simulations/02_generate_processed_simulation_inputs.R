# Generate processed simulation inputs from SCL 3D coordinates.
#
# This script implements the final simulation-generation step described in the
# manuscript: for each chromosome window, use SCL-derived 3D coordinates from
# three source cell lines to simulate 30 underlying cells, with 10 replicate
# cells from each source. The resulting `true_muS` matrices are the processed
# inputs consumed by the manuscript HiCBZIP simulation workflows.
#
# Expected input:
#   data/processed/simulation/scl_outputs/<region>/GSM..._<region>.txt
#
# Output:
#   data/processed/simulation/input/Simulation_snm3Cseq_human_brain_astrocytes_50k_chr*_..._K3X10.RData

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("data.table", "dplyr", "tibble"))
source_hicbzip_core()

scl_root <- Sys.getenv(
  "HICBZIP_SIM_SCL_OUTPUT_DIR",
  unset = path_here("data", "processed", "simulation", "scl_outputs")
)
out_dir <- Sys.getenv(
  "HICBZIP_SIM_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "input")
)
ensure_dir(out_dir)

resolution <- 50 * 1e3
max_dist <- 5e6
n_rep <- 10L

parse_M_to_bp <- function(x) {
  x <- as.character(x)
  out <- numeric(length(x))
  is_M <- grepl("M$", x, ignore.case = TRUE)
  is_K <- grepl("K$", x, ignore.case = TRUE)
  is_num <- !(is_M | is_K)
  out[is_M] <- as.numeric(sub("M$", "", x[is_M], ignore.case = TRUE)) * 1e6
  out[is_K] <- as.numeric(sub("K$", "", x[is_K], ignore.case = TRUE)) * 1e3
  out[is_num] <- as.numeric(x[is_num])
  out
}

simulation_regions <- tibble::tribble(
  ~chr, ~start, ~end,
   1L,  "50M",  "55M",
   4L, "140M", "145M",
   5L, "100M", "105M",
   8L, "135M", "140M",
  10L,  "25M",  "30M",
  11L,  "85M",  "90M",
  15L,  "80M",  "85M",
  16L,  "45M",  "50M",
  17L,  "55M",  "60M",
  22L,  "10M",  "15M"
) |>
  dplyr::mutate(
    start_bp = parse_M_to_bp(start),
    end_bp = parse_M_to_bp(end),
    chr_name = paste0("chr", chr),
    region_id = paste0(chr_name, "_", start, "_", end)
  )

read_scl_coords <- function(path) {
  coords <- data.table::fread(path, header = FALSE)
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  if (ncol(coords) < 3) {
    stop("SCL coordinate file must contain at least 3 columns: ", path, call. = FALSE)
  }
  if (ncol(coords) >= 4) {
    coords[, 2:4, drop = FALSE]
  } else {
    coords[, seq_len(3), drop = FALSE]
  }
}

for (r in seq_len(nrow(simulation_regions))) {
  reg <- simulation_regions[r, ]
  region_dir <- file.path(scl_root, reg$region_id)
  scl_files <- list.files(region_dir, pattern = "^GSM.*\\.txt$", full.names = TRUE)
  if (length(scl_files) != 3L) {
    stop("Expected 3 SCL output files in ", region_dir, "; found ", length(scl_files), call. = FALSE)
  }

  message("Generating processed simulation input for ", reg$region_id)
  coords_list <- lapply(sort(scl_files), read_scl_coords)
  max_row <- min(vapply(coords_list, nrow, integer(1)))
  coords_list <- lapply(coords_list, function(x) x[seq_len(max_row), , drop = FALSE])

  repeated_coords <- rep(coords_list, each = n_rep)
  cell_source <- rep(seq_along(coords_list), each = n_rep)
  cell_rep <- rep(seq_len(n_rep), times = length(coords_list))
  cell_names <- paste0("source", cell_source, "_rep", cell_rep)

  set.seed(123456)
  true_mu <- vapply(repeated_coords, simulate_mu_from_3Dcoords, numeric(max_row * (max_row - 1) / 2))
  colnames(true_mu) <- cell_names
  true_mu[!is.finite(true_mu)] <- 0

  N <- nrow(true_mu)
  K <- ncol(true_mu)
  true_pi_mean <- 1 - exp(true_mu) / (1 + exp(true_mu))
  eps <- 0.3
  plus <- 0.2
  true_pi <- matrix(
    runif(N * K, min = pmax(0, true_pi_mean - eps + plus), max = pmin(1, true_pi_mean + eps + plus)),
    N,
    K
  )
  colnames(true_pi) <- cell_names
  true_S <- matrix(rbinom(N * K, 1, true_pi), N, K)
  colnames(true_S) <- cell_names
  true_muS <- true_mu * (1 - true_S)
  colnames(true_muS) <- cell_names

  coverage <- 0.3
  sim_y <- matrix(rpois(N * K, lambda = coverage * as.numeric(true_muS)), N, K)
  colnames(sim_y) <- cell_names

  out_file <- file.path(
    out_dir,
    paste0(
      "Simulation_snm3Cseq_human_brain_astrocytes_",
      resolution / 1000, "k_chr", reg$chr, "_",
      reg$start_bp / 1e6, "M_", reg$end_bp / 1e6, "M_md",
      max_dist / 1e6, "M_K3X10.RData"
    )
  )

  save(
    resolution, max_dist, n_rep, N, K, coverage,
    true_mu, true_pi, true_S, true_muS, sim_y,
    cell_names, cell_source, cell_rep,
    file = out_file
  )
  message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
}

message("Processed simulation input generation complete.")
