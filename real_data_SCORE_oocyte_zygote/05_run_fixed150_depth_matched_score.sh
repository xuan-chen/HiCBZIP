#!/usr/bin/env bash
# Run the finalized imputation pipelines with the original SCORE flags and seeds.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${script_dir}/.." && pwd)"
out_dir="${HICBZIP_SCORE_FIXED150_OUTPUT:-${root}/results/SCORE_oocyte_zygote/fixed150}"
out_dir="$(cd "${out_dir}" && pwd)"
source_dir="${HICBZIP_SCORE_SOURCE:-${root}/data/processed/SCORE_oocyte_zygote/source_inputs}"
source_dir="$(cd "${source_dir}" && pwd)"
score_bin="$(command -v "${SCORE_BIN:-score}")"
score_bin="$(cd "$(dirname "${score_bin}")" && pwd)/$(basename "${score_bin}")"
anchors="${source_dir}/mm10.genome_split_1M"
n_runs="${HICBZIP_FIXED150_RUNS:-10}"
[[ -s "${anchors}" && "${n_runs}" =~ ^([1-9]|10)$ ]] || {
  echo "Missing anchor file or invalid run count (expected 1 through 10)." >&2
  exit 1
}
if (( $# == 0 )); then set -- gbnb nm schicimpute; fi

for model in "$@"; do
  case "${model}" in
    gbnb|nm)
      cohort=fixed150
      reference="${out_dir}/metadata/oocyte_zygote_ref_fixed150.tsv"
      ;;
    schicimpute)
      cohort=complete146
      reference="${out_dir}/metadata/oocyte_zygote_ref_schicimpute_complete146.tsv"
      ;;
    *) echo "Unknown model: ${model}" >&2; exit 1 ;;
  esac
  [[ -s "${reference}" ]]
  if [[ "${model}" == schicimpute ]]; then
    dataset=oocyte_zygote_schicimpute_complete_case_146_depth_matched
  else
    dataset="oocyte_zygote_${model}_depth_matched_fixed150"
  fi
  pairs="${out_dir}/depth_matched_inputs/pair_files"
  scool="${out_dir}/depth_matched_inputs/scool/${dataset}_1M.scool"
  mkdir -p "$(dirname "${scool}")"
  # Use the same score cooler input as the final SnapATAC run.
  if [[ ! -s "${scool}" ]]; then
    "${score_bin}" cooler --dset "${dataset}" --resolution 1M \
      --data_dir "${pairs}/depth_matched_integer/${model}" \
      --anchor_file "${anchors}" --reference "${reference}" --out "${scool}"
  fi
  [[ -s "${scool}" ]]

  for workflow in innerproduct snapatac; do
    work="${out_dir}/outputs/${cohort}_score/external_imputation/${workflow}/depth_matched/${model}"
    mkdir -p "${work}/results" "${work}/logs"
    if [[ "${workflow}" == innerproduct ]]; then
      args=(--data_dir "${pairs}/depth_matched_continuous/${model}"
            --anchor_file "${anchors}" --embedding_algs InnerProduct --no_cache)
    else
      args=(--scool "${scool}" --embedding_algs SnapATAC --no_viz --snapatac_no_idf)
    fi

    for ((run=1; run<=n_runs; run++)); do
      printf -v run_id '%02d' "${run}"
      dset="oocyte_zygote_${model}_depth_matched_${cohort}_${workflow}_r${run_id}"
      metrics="${work}/results/${dset}/1M/${workflow}/${workflow}.json"
      if [[ -s "${metrics}" ]]; then continue; fi
      (cd "${work}" && "${score_bin}" embed \
        --dset "${dset}" --out "${work}/results" --resolution 1M \
        "${args[@]}" --reference "${reference}" --min_depth 5000 \
        --ignore_chr_filter --use_xy --seed "$((2025 + run))") \
        > "${work}/logs/run_${run_id}.stdout.log" \
        2> "${work}/logs/run_${run_id}.stderr.log"
      [[ -s "${metrics}" ]] || { echo "Missing SCORE metrics: ${metrics}" >&2; exit 1; }
    done
  done
done
