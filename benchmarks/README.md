# External Benchmark Methods

This directory documents the external methods used for comparison in the HiCBZIP manuscript. The repository focuses on the HiCBZIP implementation and final reproducibility workflows; external packages are run through their own software distributions and consumed as processed outputs by the summary scripts.

## Provenance Summary

| Method | Confirmed software record | Consumed by workflow |
| --- | --- | --- |
| HiCImpute | R package version 1.0. | Simulation summaries and NPC chrX summaries. |
| Higashi | Higashi 0.1.0a0 in the recorded `server_higashi_env_simulation` environment. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |
| Fast-Higashi | Fast-Higashi 0.1.1a0 in the recorded `server_fast_higashi_env_simulation` environment. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |
| SCORE | Python environment records are provided for SCORE and SCORE/Higashi-compatible runs. | Real-data study 2. |
| scHiCluster | Python package `schicluster` 1.3.5.dev22+gd566046 in the recorded `schicluster` conda environment. No `scHiCluster` or `schicluster` command-line executable was present on PATH in that environment. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |

## Simulation Benchmark Inputs and Outputs

The simulation summary workflows consume processed matrix objects rather than full external-method working directories:

| Method | Input to public summary workflow | Notes |
| --- | --- | --- |
| HiCImpute | `data/processed/simulation/muS_combined_plot_260124.RData` | Generated from the processed simulation inputs using `HiCImpute::MCMCImpute()`. Manuscript settings used `n = 100`, `cutoff = 0.5`, `niter = 1000`, `burnin = 500`, and `mc.cores = 1`. |
| scHiCluster | `data/processed/simulation/muS_combined_plot_260124.RData` | Contact-pair files were generated for each chromosome window, coverage level, and cell; scHiCluster imputed HDF5 outputs were converted back to long-form matrices for summary. |
| Higashi | `data/processed/simulation/higashi_unified_allchr_260208.RData` | Higashi outputs were exported as long-form matrices and merged with the simulation summary workflow. The manuscript comparison used the `Higashi(nbr5)` output as `Higashi`. |
| Fast-Higashi | `data/processed/simulation/muS_combined_plot_both_fasthigashi_260613.RData` | Fast-Higashi was run with partial RWR enabled, rank 3, and convolution enabled only for coverage values below 0.1 in the validated simulation workflow. |

## Simulation Command Templates

HiCImpute simulation runs used the following R-level call pattern:

```r
HiCImpute::MCMCImpute(
  scHiC = scHiC,
  bulk = bulk_vec,
  expected = NULL,
  n = 100,
  mc.cores = 1,
  cutoff = 0.5,
  niter = 1000,
  burnin = 500
)
```

The scHiCluster simulation workflow used per-cell contact-pair text files as input and collected imputed `.hdf5` matrices from each coverage-specific output directory. The recorded server environment check was:

```bash
conda activate schicluster
which scHiCluster
which schicluster
conda list | grep -i "schicluster\|hicluster"
pip list | grep -i "schicluster\|hicluster"
```

This confirmed package `schicluster` version `1.3.5.dev22+gd566046`; the executable command was not available on PATH.

Higashi simulation inputs were prepared as `data.txt`, `label.txt`, and `config.json` run folders for each chromosome and coverage. The export step collected `ori`, `nbr0`, and `nbr5` matrices from each run; the manuscript summaries use `nbr5`.

Fast-Higashi simulation runs used the following command pattern after preparing per-run folders:

```bash
bash run_all_fasthigashi.sh <runs_root> <log_dir>
bash export_all_fasthigashi.sh <runs_root> <log_dir>
```

The run wrapper called `FastHigashi(..., do_rwr = TRUE, do_conv = coverage < 0.1, filter = FALSE, do_col = FALSE)` followed by `fast_process_data()`, `prep_dataset(batch_norm = TRUE)`, `run_model(rank = 3)`, and `only_partial_rwr()`.

## SCORE command template

```bash
score cooler \
  --dset <dataset_name> \
  --data_dir <pair_file_directory> \
  --anchor_file <genome_split_anchor_file> \
  --reference <cell_reference_tsv> \
  --resolution 1M \
  --out <output.scool>

score embed \
  --dset <dataset_name> \
  --resolution 1M \
  --scool <input.scool> \
  --reference <cell_reference_tsv> \
  --embedding_algs InnerProduct \
  --min_depth 5000 \
  --use_xy \
  --seed <seed>
```

The scripted SCORE calls used by the reproducibility workflow are provided in:

- `real_data_SCORE_oocyte_zygote/run_SCORE_innerproduct_four_inputs_10runs.R`
- `real_data_SCORE_oocyte_zygote/run_SCORE_snapatac_four_inputs_10runs.R`
