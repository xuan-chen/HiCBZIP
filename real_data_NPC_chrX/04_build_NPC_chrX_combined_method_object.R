# Build the unified NPC chrX object used by the manuscript summary scripts.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "tibble"))
source_hicbzip_core()
source(file.path(script_dir, "NPC_chrX_workflow_helpers.R"))

input_file <- Sys.getenv(
  "HICBZIP_REAL1_INPUT_RDATA",
  unset = path_here("data", "processed", "NPC_chrX", "data_NPC250k_0h_X_full.RData")
)
sim_y_file <- Sys.getenv(
  "HICBZIP_REAL1_SIM_Y_LIST",
  unset = path_here("results", "NPC_chrX", "benchmark_inputs", "sim_y_list_NPC_chrX.RData")
)
gb_file <- path_here("results", "NPC_chrX", "HiCBZIP_GB_NB", "o.bziphic_gbnb.muS.full.AllCoverage.RData")
hicimpute_file <- path_here("results", "NPC_chrX", "HiCImpute", "o.hicimpute.muS.full.AllCoverage.RData")
schicluster_file <- path_here("results", "NPC_chrX", "scHiCluster", "o.hicluster.muS.full.AllCoverage.RData")
higashi_file <- path_here("results", "NPC_chrX", "Higashi", "higashi_unified_chrX_NPC250k_diag.RData")
fasthigashi_file <- path_here("results", "NPC_chrX", "FastHigashi", "fasthigashi_unified_chrX_NPC250k_diag.RData")
out_file <- Sys.getenv(
  "HICBZIP_REAL1_UNIFIED_OUTPUT",
  unset = path_here("data", "processed", "NPC_chrX", "list_muS_unified_260614.RData")
)

require_files(c(input_file, sim_y_file, gb_file, hicimpute_file, schicluster_file, higashi_file, fasthigashi_file),
              label = "NPC chrX combined workflow input")

load(input_file)
load(sim_y_file)
load(gb_file)
load(hicimpute_file)
load(schicluster_file)
load(higashi_file)
load(fasthigashi_file)

if (!exists("true_muS")) true_muS <- replicate(K, bulk)
coverage_labels <- npc_chrX_coverage_labels()
lambda_vals <- npc_chrX_coverage_values()
lambda_list <- paste0("ld=", coverage_labels)
resolution <- if (exists("resolution")) resolution else 250000L
N <- nrow(true_muS)
K <- ncol(true_muS)
n <- if (exists("n")) n else npc_chrX_infer_n_bins_from_diag_long(N)

raw_list <- npc_chrX_standardize_method_list(sim_y_list, coverage_labels, "Raw")
gb_list <- npc_chrX_standardize_method_list(muS_bzip_gb_nb, coverage_labels, "HiCBZIP-GB(NB)")
hicimpute_list <- npc_chrX_standardize_method_list(o.hicimpute.muS.full, coverage_labels, "HiCImpute")
schicluster_list <- npc_chrX_standardize_method_list(o.hicluster.muS.full, coverage_labels, "scHiCluster")

nm_dir <- path_here("results", "NPC_chrX", "HiCBZIP_NM")
nm_files <- file.path(nm_dir, paste0("CMDSTAN_mclapply_BHZIP_matchN_NPC250k_0h_X_full_", coverage_labels, ".RData"))
require_files(nm_files, label = "HiCBZIP-N(M) per-coverage results")
nm_list <- lapply(nm_files, function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("res_list", envir = env)) stop("Missing res_list in: ", path, call. = FALSE)
  do.call(rbind, lapply(env$res_list, function(x) x$mu_tilde))
})
nm_list <- npc_chrX_standardize_method_list(nm_list, coverage_labels, "HiCBZIP-N(M)")

extract_unified <- function(tbl, model_name) {
  rows <- tbl[tbl$model == model_name, , drop = FALSE]
  rows$coverage <- as.character(rows$coverage)
  lapply(coverage_labels, function(cov_label) {
    hit <- rows[rows$coverage == cov_label, , drop = FALSE]
    if (nrow(hit) != 1) stop("Expected one ", model_name, " row for coverage ", cov_label, call. = FALSE)
    hit$muS[[1]]
  }) |>
    npc_chrX_standardize_method_list(coverage_labels = coverage_labels, method_name = model_name)
}

higashi_list <- extract_unified(higashi_unified_chrX, "Higashi(nbr5)")
fasthigashi_list <- extract_unified(fasthigashi_unified_chrX, "Fast-Higashi")

list_muS_combine <- c(raw_list, gb_list, nm_list, schicluster_list, hicimpute_list, higashi_list, fasthigashi_list)

expected_dim <- c(N, K)
bad <- names(list_muS_combine)[!vapply(list_muS_combine, function(x) identical(dim(x), expected_dim), logical(1))]
if (length(bad) > 0) {
  stop("These matrices do not match expected dimension ", paste(expected_dim, collapse = " x "), ":\n- ",
       paste(bad, collapse = "\n- "), call. = FALSE)
}

scc_muS <- npc_chrX_build_scc_table(list_muS_combine, true_muS, K, lambda_list, resolution)
IS_vector_plot <- npc_chrX_build_is_table(list_muS_combine, true_muS, K, lambda_list, resolution)
IS_vector_plot2 <- IS_vector_plot
sparsity_list <- coverage_labels
sparsity_labels <- lambda_list
bulk <- if (exists("bulk")) bulk else rowMeans(true_muS)

ensure_dir(dirname(out_file))
save(list_muS_combine, true_muS, bulk, lambda_list, lambda_vals, n, N, resolution,
     sparsity_labels, sparsity_list, scc_muS, IS_vector_plot, IS_vector_plot2,
     file = out_file)
message("Saved unified NPC chrX object: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
