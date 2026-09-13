#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

RUNS_ROOT="${1:-results/simulation/benchmark_methods/inputs/FastHigashi/runs}"
LOG_DIR="${2:-results/simulation/benchmark_methods/logs/FastHigashi_export}"
mkdir -p "${LOG_DIR}"

for run_dir in "${RUNS_ROOT}"/chr*_cov*; do
  [[ -d "${run_dir}" ]] || continue
  run_name="$(basename "${run_dir}")"
  log_file="${LOG_DIR}/${run_name}.log"
  if [[ ! -s "${run_dir}/impute_prwr.hdf5" ]]; then
    echo "skip missing Fast-Higashi output: ${run_name}"
    continue
  fi
  if compgen -G "${run_dir}/export_long/fasthigashi_prwr_*.csv" >/dev/null; then
    echo "skip existing export: ${run_name}"
    continue
  fi
  echo "export Fast-Higashi: ${run_name}"
  python simulation/benchmark_methods/08_export_FastHigashi_simulation.py --run_dir "${run_dir}" > "${log_file}" 2>&1
done
