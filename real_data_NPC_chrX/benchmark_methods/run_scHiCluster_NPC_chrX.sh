#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT="${1:-results/NPC_chrX/benchmark_inputs/scHiCluster}"
RES="${2:-250000}"
CHROM_SIZES="${3:-}"
CHROM="X"

if [[ -z "${CHROM_SIZES}" ]]; then
  echo "Usage: bash real_data_NPC_chrX/benchmark_methods/run_scHiCluster_NPC_chrX.sh <scHiCluster_input_root> <resolution_bp> <chrom_sizes_file>"
  echo "Example: bash real_data_NPC_chrX/benchmark_methods/run_scHiCluster_NPC_chrX.sh results/NPC_chrX/benchmark_inputs/scHiCluster 250000 /path/to/hg19.chrom.sizes.txt"
  exit 1
fi

if ! command -v hicluster >/dev/null 2>&1; then
  echo "Missing executable: hicluster"
  echo "Activate the conda environment that contains the schicluster package before running this script."
  exit 1
fi

echo "Input root: ${ROOT}"
echo "Resolution: ${RES}"
echo "Chrom sizes: ${CHROM_SIZES}"
echo "hicluster: $(command -v hicluster)"

cov_dirs=( "${ROOT}"/NPC250k_0h_X_full_* )
if (( ${#cov_dirs[@]} == 0 )); then
  echo "No coverage directories matched: ${ROOT}/NPC250k_0h_X_full_*"
  exit 1
fi

for cov_dir in "${cov_dirs[@]}"; do
  contact_dir="${cov_dir}/contact_pair"
  out_dir="${cov_dir}/imputed_matrix"
  mkdir -p "${out_dir}"
  files=( "${contact_dir}"/*.txt )
  echo "Coverage directory: $(basename "${cov_dir}")"
  echo "Contact files: ${#files[@]}"
  if (( ${#files[@]} == 0 )); then
    echo "Skipping empty directory: ${contact_dir}"
    continue
  fi

  for f in "${files[@]}"; do
    base="$(basename "${f}")"
    cell="${base%_chrX.txt}"
    out_h5="${out_dir}/${cell}_chrX.hdf5"
    if [[ -s "${out_h5}" ]]; then
      echo "  skip existing: ${cell}"
      continue
    fi
    echo "  run: ${cell}"
    hicluster impute-cell \
      --indir "${contact_dir}" \
      --outdir "${out_dir}" \
      --cell "${cell}" \
      --chrom "${CHROM}" \
      --res "${RES}" \
      --chrom_file "${CHROM_SIZES}"
  done
done
