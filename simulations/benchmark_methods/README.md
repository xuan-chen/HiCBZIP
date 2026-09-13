# Simulation Benchmark Methods

This folder contains study-specific workflows for benchmark methods used in the HiCBZIP simulation study.

The main simulation workflow in `simulation/` focuses on HiCBZIP method code and manuscript-level summaries. These scripts document how the benchmark outputs can be regenerated from the processed simulation inputs when needed.

## Inputs

Place processed simulation input files under:

```text
data/processed/simulation/input/
```

Expected file pattern:

```text
Simulation_snm3Cseq_human_brain_astrocytes_50k_chr*_K3X10.RData
```

Each file must contain `true_muS`, an off-diagonal long-vector matrix with cells in columns.

## Scripts

| Script | Purpose |
| --- | --- |
| `01_prepare_simulation_benchmark_inputs.R` | Creates per-coverage inputs for scHiCluster, Higashi, and Fast-Higashi. |
| `02_run_HiCImpute_simulation.R` | Runs HiCImpute across the manuscript chromosomes and coverage levels. |
| `03_run_scHiCluster_simulation.sh` | Runs scHiCluster from prepared contact-pair files. |
| `04_collect_scHiCluster_simulation.R` | Collects scHiCluster HDF5 outputs into one RDS object. |
| `05_run_Higashi_simulation.py` | Runs Higashi for one prepared chromosome/coverage directory. |
| `05b_run_all_Higashi_simulation.sh` | Runs Higashi across all prepared simulation directories. |
| `06_collect_Higashi_simulation.R` | Collects Higashi CSV outputs into one RData object. |
| `07_run_FastHigashi_simulation.py` | Runs Fast-Higashi for one prepared chromosome/coverage directory. |
| `07b_run_all_FastHigashi_simulation.sh` | Runs Fast-Higashi across all prepared simulation directories. |
| `08_export_FastHigashi_simulation.py` | Exports Fast-Higashi HDF5 output to a long-vector CSV file. |
| `08b_export_all_FastHigashi_simulation.sh` | Exports all completed Fast-Higashi simulation runs. |
| `09_collect_FastHigashi_simulation.R` | Collects Fast-Higashi CSV outputs into one RData object. |

External method versions and environments are documented in `environment/` and in the top-level `benchmarks/README.md`.

## Notes

The archived manuscript data object can be used for final figure and metric reproduction without rerunning these benchmark methods. Rerunning all benchmark methods is computationally heavier and may require method-specific conda environments.
