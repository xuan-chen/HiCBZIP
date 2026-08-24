# Create the processed NPC chrX input object used by the NPC chrX recovery study.
#
# Expected input:
#   data/processed/NPC_chrX/NPC250k_0h_X.mat
#
# Output:
#   data/processed/NPC_chrX/data_NPC250k_0h_X_full.RData

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr"))
source_hicbzip_core()

input_mat <- Sys.getenv(
  "HICBZIP_REAL1_NPC_MAT",
  unset = path_here("data", "processed", "NPC_chrX", "NPC250k_0h_X.mat")
)
out_dir <- path_here("data", "processed", "NPC_chrX")
out_file <- file.path(out_dir, "data_NPC250k_0h_X_full.RData")

require_files(input_mat, label = "NPC chrX contact matrix")
ensure_dir(out_dir)

pop0h <- read.table(input_mat)
mat <- as.matrix(pop0h[, -1])
mat[is.na(mat)] <- 0

bulk <- matrix2D_to_matrix_long(mat, include.diag = TRUE)
N <- length(bulk)
K <- as.integer(Sys.getenv("HICBZIP_REAL1_N_CELLS", unset = "30"))
n <- nrow(mat)

set.seed(123456)
coverage <- as.numeric(Sys.getenv("HICBZIP_REAL1_EXAMPLE_COVERAGE", unset = "0.01"))
true_muS <- replicate(K, bulk)
sim_y <- matrix(rpois(N * K, lambda = coverage * true_muS), nrow = N, ncol = K)

save(pop0h, mat, bulk, N, K, n, coverage, true_muS, sim_y, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
