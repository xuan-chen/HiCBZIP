# Simulation study

Ordered public workflow:

0. To regenerate processed simulation inputs from raw indexed-contact files:
   - Place the three source indexed-contact files under `data/raw/simulation/Raw_Data/`.
   - Run `00_prepare_simulation_source_contacts.R` to extract the ten chromosome windows.
   - Run SCL with `01_run_scl_for_simulation.sh`.
   - Run `02_generate_processed_simulation_inputs.R` to create the processed `true_muS` inputs.
1. Alternatively, place archived processed simulation RData inputs under `data/processed/simulation/input/`.
2. Run `03_run_HiCBZIP_GB_NB_simulation.R`.
3. Run `04_run_HiCBZIP_NGS_simulation.R`.
4. Run `05_run_HiCBZIP_NM_simulation.R`.
5. Add processed non-regenerated external benchmark outputs under `data/processed/simulation/`.
6. Run `06_prepare_simulation_manuscript_summaries.R`.

Included manuscript summary notebooks:

- `summarize_simulation_metrics.Rmd`
- `summarize_simulation_clustering.Rmd`
- `make_simulation_heatmap_panels.Rmd`

## SCL preprocessing

SCL was used only to infer 3D chromosome coordinates for simulation input generation. The source folder used for this project was `SCL1.0_source_code`, from the public SCL paper site. Its README gives the build command:

```bash
cd SCL1.0_source_code/scripts/
g++ -o scl scl.cpp
```

The manuscript workflow uses:

```bash
SCL_BIN=/path/to/SCL1.0_source_code/scripts/scl bash simulation/01_run_scl_for_simulation.sh
```

The SCL source/executable record and hashes are stored in `environment/scl_source_record.txt`.
