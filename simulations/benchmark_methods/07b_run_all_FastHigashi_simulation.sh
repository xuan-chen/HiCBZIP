#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

RUNS_ROOT="${1:-results/simulation/benchmark_methods/inputs/FastHigashi/runs}"
LOG_DIR="${2:-results/simulation/benchmark_methods/logs/FastHigashi}"
RANK="${FAST_HIGASHI_RANK:-3}"
mkdir -p "${LOG_DIR}"

for run_dir in "${RUNS_ROOT}"/chr*_cov*; do
  [[ -d "${run_dir}" ]] || continue
  run_name="$(basename "${run_dir}")"
  log_file="${LOG_DIR}/${run_name}.log"
  if [[ -s "${run_dir}/impute_prwr.hdf5" ]]; then
    echo "skip existing Fast-Higashi output: ${run_name}"
    continue
  fi
  echo "run Fast-Higashi: ${run_name}"
  python simulation/benchmark_methods/07_run_FastHigashi_simulation.py --run_dir "${run_dir}" --rank "${RANK}" > "${log_file}" 2>&1
done
