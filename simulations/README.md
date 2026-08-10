# Simulation study

Ordered public workflow:

1. Place processed simulation RData inputs under `data/processed/simulation/input/`.
2. Run `01_run_HiCBZIP_GB_NB_simulation.R`.
3. Add processed non-regenerated method outputs under `data/processed/simulation/`.
4. Run `02_prepare_simulation_manuscript_summaries.R`.

Included manuscript summary notebooks:

- `summarize_simulation_metrics.Rmd`
- `summarize_simulation_clustering.Rmd`
- `make_simulation_heatmap_panels.Rmd`

Important unresolved item: the old simulation N(M)/GM runner references `BHZIP_GM_nocov_250704.stan`, which is not present in the candidate code folder. Add that Stan model if it produced a reported manuscript result.
