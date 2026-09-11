# Run HiCImpute for NPC chrX manuscript downsampling levels.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("HiCImpute"))
source_hicbzip_core()
source(file.path(script_dir, "..", "NPC_chrX_workflow_helpers.R"))

input_file <- Sys.getenv(
  "HICBZIP_REAL1_INPUT_RDATA",
  unset = path_here("data", "processed", "NPC_chrX", "data_NPC250k_0h_X_full.RData")
)
sim_y_file <- Sys.getenv(
  "HICBZIP_REAL1_SIM_Y_LIST",
  unset = path_here("results", "NPC_chrX", "benchmark_inputs", "sim_y_list_NPC_chrX.RData")
)
out_dir <- path_here("results", "NPC_chrX", "HiCImpute")
out_file <- file.path(out_dir, "o.hicimpute.muS.full.AllCoverage.RData")

require_files(input_file, label = "NPC chrX processed input")
load(input_file)
if (file.exists(sim_y_file)) {
  load(sim_y_file)
} else {
  if (!exists("true_muS")) true_muS <- replicate(K, bulk)
  sim_y_list <- npc_chrX_make_sim_y_list(true_muS)
  lambda_vals <- npc_chrX_coverage_values()
}

n_bins <- if (exists("n")) n else npc_chrX_infer_n_bins_from_diag_long(nrow(sim_y_list[[1]]))
ensure_dir(out_dir)

o.hicimpute.muS.full <- lapply(seq_along(sim_y_list), function(ii) {
  coverage <- as.numeric(lambda_vals[[ii]])
  scHiC <- sim_y_list[[ii]]
  HiCImpute::MCMCImpute(
    scHiC = scHiC,
    startval = c(10, 1, 10, 1, 10, 0.5, 10, 0.5, 0, replicate(ncol(scHiC), 1)),
    n = n_bins,
    mc.cores = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_WORKERS", unset = "1")),
    cutoff = 0.5,
    niter = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_NITER", unset = "1000")),
    burnin = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_BURNIN", unset = "500"))
  )
})
names(o.hicimpute.muS.full) <- paste0("HiCImpute,ld=", npc_chrX_coverage_labels())

save(o.hicimpute.muS.full, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
