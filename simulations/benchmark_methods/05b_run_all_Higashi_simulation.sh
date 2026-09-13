#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

RUNS_ROOT="${1:-results/simulation/benchmark_methods/inputs/Higashi/runs}"
LOG_DIR="${2:-results/simulation/benchmark_methods/logs/Higashi}"
mkdir -p "${LOG_DIR}"

for run_dir in "${RUNS_ROOT}"/chr*_cov*; do
  [[ -d "${run_dir}" ]] || continue
  run_name="$(basename "${run_dir}")"
  log_file="${LOG_DIR}/${run_name}.log"
  if compgen -G "${run_dir}/export_long/ori_*.csv" >/dev/null; then
    echo "skip existing export: ${run_name}"
    continue
  fi
  echo "run Higashi: ${run_name}"
  python simulation/benchmark_methods/05_run_Higashi_simulation.py --run_dir "${run_dir}" > "${log_file}" 2>&1
done
