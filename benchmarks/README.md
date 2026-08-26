# External Benchmark Methods

This directory documents the external methods used for comparison in the HiCBZIP manuscript. The repository focuses on the HiCBZIP implementation and final reproducibility workflows; external packages are run through their own software distributions and consumed as processed outputs by the summary scripts.

## Provenance Summary

| Method | Documentation included here | Consumed by workflow |
| --- | --- | --- |
| scHiCluster | Input/output format and command template for generated matrix or JSON summaries. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |
| HiCImpute | Input/output format and command template for imputed contact matrices. | Simulation summaries and NPC chrX summaries. |
| Higashi | Environment record, configuration fields, and export/conversion template. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |
| Fast-Higashi | Environment record, configuration fields, and export/conversion template. | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder. |
| SCORE | Environment record and command template for `score cooler`, `score embed InnerProduct`, and `score embed SnapATAC`. | Real-data study 2. |

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
