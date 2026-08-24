#!/usr/bin/env bash
set -euo pipefail

# Template for running SCL on the simulation chromosome-window pair files.
#
# Run from the repository root after `00_prepare_simulation_source_contacts.R`.
# Set SCL_BIN to the SCL executable. The SCL resolution argument is 0.05,
# matching 50 kb windows in the original manuscript workflow.
#
# Example:
#   SCL_BIN=/path/to/scl bash simulation/01_run_scl_for_simulation.sh

SCL_BIN="${SCL_BIN:-scl}"
INPUT_ROOT="${HICBZIP_SIM_SCL_INPUT_DIR:-data/processed/simulation/scl_inputs}"
OUTPUT_ROOT="${HICBZIP_SIM_SCL_OUTPUT_DIR:-data/processed/simulation/scl_outputs}"
RESOLUTION="${HICBZIP_SCL_RESOLUTION:-0.05}"

mkdir -p "${OUTPUT_ROOT}"

find "${INPUT_ROOT}" -type f -name "GSM*.txt" | sort | while read -r input_file; do
  rel="${input_file#${INPUT_ROOT}/}"
  out_file="${OUTPUT_ROOT}/${rel}"
  mkdir -p "$(dirname "${out_file}")"
  echo "Running SCL: ${input_file} -> ${out_file}"
  "${SCL_BIN}" -i "${input_file}" -o "${out_file}" -res "${RESOLUTION}"
done

echo "SCL runs complete. Next run simulation/02_generate_processed_simulation_inputs.R."
