#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT="${1:-results/simulation/benchmark_methods/inputs/scHiCluster}"
RES="${2:-50000}"
CHROM_SIZES="${3:-}"

if [[ -z "${CHROM_SIZES}" ]]; then
  echo "Usage: bash simulation/benchmark_methods/03_run_scHiCluster_simulation.sh <scHiCluster_input_root> <resolution_bp> <chrom_sizes_file>"
  echo "Example: bash simulation/benchmark_methods/03_run_scHiCluster_simulation.sh results/simulation/benchmark_methods/inputs/scHiCluster 50000 /path/to/hg19.chrom.sizes.txt"
  exit 1
fi

if ! command -v hicluster >/dev/null 2>&1; then
  echo "Missing executable: hicluster"
  echo "Activate the conda environment that contains the schicluster package before running this script."
  exit 1
fi

for cov_dir in "${ROOT}"/sim_HBA3_random_chr_*; do
  [[ -d "${cov_dir}/contact_pair" ]] || continue
  contact_dir="${cov_dir}/contact_pair"
  out_dir="${cov_dir}/imputed_matrix"
  mkdir -p "${out_dir}"
  echo "Coverage directory: $(basename "${cov_dir}")"

  for f in "${contact_dir}"/*.txt; do
    base="$(basename "${f}")"
    cell="${base%_chr*.txt}"
    chrom="$(echo "${base}" | sed -E 's/.*_chr([0-9XY]+)\.txt/\1/')"
    out_h5="${out_dir}/${cell}_chr${chrom}.hdf5"
    if [[ -s "${out_h5}" ]]; then
      echo "  skip existing: ${cell} chr${chrom}"
      continue
    fi
    echo "  run: ${cell} chr${chrom}"
    hicluster impute-cell \
      --indir "${contact_dir}" \
      --outdir "${out_dir}" \
      --cell "${cell}" \
      --chrom "${chrom}" \
      --res "${RES}" \
      --chrom_file "${CHROM_SIZES}"
  done
done
