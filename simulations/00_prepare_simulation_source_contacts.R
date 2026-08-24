# Extract chromosome-window contact pairs for SCL.
#
# This script prepares the SCL input files used to build the biologically
# informed simulation. It expects three indexed-contact files, one per source
# cell line, under `data/raw/simulation/Raw_Data/`.
#
# Expected raw files:
#   GSM3749700_190305_snm3Cseq_hs_29yr_BA10_UMB5580_1_UMB5580_2_A11_AD004_indexed_contacts.txt
#   GSM3750251_190315_snm3Cseq_hs_21yr_BA10_UMB5577_3_UMB5577_4_A9_AD010_indexed_contacts.txt
#   GSM3751478_190315_snm3Cseq_hs_29yr_BA10_UMB5580_3_UMB5580_4_C10_AD006_indexed_contacts.txt
#
# Output:
#   data/processed/simulation/scl_inputs/<region>/<GSM...>_<region>.txt

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("data.table", "dplyr", "tibble"))

raw_dir <- Sys.getenv(
  "HICBZIP_SIM_RAW_DIR",
  unset = path_here("data", "raw", "simulation", "Raw_Data")
)
out_root <- Sys.getenv(
  "HICBZIP_SIM_SCL_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "scl_inputs")
)

source_files <- c(
  "GSM3749700_190305_snm3Cseq_hs_29yr_BA10_UMB5580_1_UMB5580_2_A11_AD004_indexed_contacts.txt",
  "GSM3750251_190315_snm3Cseq_hs_21yr_BA10_UMB5577_3_UMB5577_4_A9_AD010_indexed_contacts.txt",
  "GSM3751478_190315_snm3Cseq_hs_29yr_BA10_UMB5580_3_UMB5580_4_C10_AD006_indexed_contacts.txt"
)

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

require_files(file.path(raw_dir, source_files), label = "simulation raw indexed-contact file")
ensure_dir(out_root)

for (src in source_files) {
  raw_path <- file.path(raw_dir, src)
  gsm <- sub("^(GSM[0-9]+).*", "\\1", src)
  message("Reading raw contact file: ", basename(raw_path))

  contacts <- data.table::fread(
    raw_path,
    header = FALSE,
    select = c(2, 3, 4, 5),
    col.names = c("chrom1", "pos1", "chrom2", "pos2")
  )

  for (r in seq_len(nrow(simulation_regions))) {
    reg <- simulation_regions[r, ]
    out_dir <- file.path(out_root, reg$region_id)
    ensure_dir(out_dir)

    region_pairs <- contacts |>
      dplyr::filter(
        chrom1 == reg$chr_name,
        chrom2 == reg$chr_name,
        pos1 >= reg$start_bp,
        pos1 <= reg$end_bp,
        pos2 >= reg$start_bp,
        pos2 <= reg$end_bp
      ) |>
      dplyr::transmute(pos1 = as.integer(pos1), pos2 = as.integer(pos2))

    out_file <- file.path(out_dir, paste0(gsm, "_", reg$region_id, ".txt"))
    data.table::fwrite(region_pairs, out_file, sep = "\t", col.names = FALSE)
    message("Wrote ", nrow(region_pairs), " pairs: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  }
}

message("SCL input extraction complete. Next run simulation/01_run_scl_for_simulation.sh.")
