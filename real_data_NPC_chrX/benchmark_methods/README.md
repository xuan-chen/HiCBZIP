# NPC chrX Benchmark Methods

These scripts document how benchmark-method outputs for the NPC chrX study can be regenerated from prepared per-coverage inputs.

They are not required for the main public workflow when starting from the archived unified object:

`data/processed/NPC_chrX/list_muS_unified_260614.RData`

The main repository workflow keeps the focus on HiCBZIP method code, the archived manuscript-scale data object, and the final summary scripts. These study-specific benchmark-method scripts are included as provenance for reviewers who want to inspect the benchmark execution pattern.

## Contents

| Script | Purpose |
| --- | --- |
| `prepare_NPC_chrX_benchmark_inputs.R` | Creates per-coverage inputs for HiCImpute, scHiCluster, Higashi, and Fast-Higashi from `data_NPC250k_0h_X_full.RData`. |
| `run_HiCImpute_NPC_chrX.R` | Runs HiCImpute across the reported downsampling levels. |
| `run_scHiCluster_NPC_chrX.sh` | Command template for running scHiCluster through the `hicluster impute-cell` interface. |
| `collect_scHiCluster_NPC_chrX.R` | Converts scHiCluster HDF5 outputs into the long-vector format used by the manuscript summaries. |
| `run_Higashi_NPC_chrX.py` | Runs Higashi for one prepared coverage directory. |
| `collect_Higashi_NPC_chrX.R` | Collects exported Higashi CSV matrices into one RData object. |
| `run_FastHigashi_NPC_chrX.py` | Runs Fast-Higashi for one prepared coverage directory. |
| `export_FastHigashi_NPC_chrX.py` | Exports Fast-Higashi HDF5 outputs to diagonal-included long-vector CSV files. |
| `collect_FastHigashi_NPC_chrX.R` | Collects exported Fast-Higashi CSV matrices into one RData object. |

Software versions and environment records are documented under `environment/` and `benchmarks/README.md`.
