# Minimal HiCBZIP-GB(NB) example using a tiny synthetic count matrix.
#
# Run from the repository root after sourcing the core method code:
# source("code_clean/HiCBZIP/BHZIP_EB.R")
# source("code_clean/examples/run_minimal_HiCBZIP_example.R")

set.seed(1)

library(dplyr)
library(magrittr)

# Six locus pairs from a 3-bin contact map, across four cells.
example_counts <- matrix(
  c(
    0, 1, 0, 2,
    3, 2, 4, 1,
    0, 0, 1, 0,
    5, 4, 3, 6,
    0, 0, 0, 1,
    2, 1, 0, 2
  ),
  nrow = 6,
  byrow = TRUE
)

colnames(example_counts) <- paste0("cell", seq_len(ncol(example_counts)))

imputed <- run_BZIP_GB_NB(
  sim_y = example_counts,
  coverage = 0.2,
  r = 1,
  threshold = 3,
  B = 10,
  include_diag = TRUE,
  return_hyper = TRUE
)

print(imputed$muS)
