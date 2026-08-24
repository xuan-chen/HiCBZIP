# External benchmark methods

The public repository documents external benchmark methods by version and command notes. Full local runner folders are not included unless a reviewer specifically asks for them.

## Required notes to fill before release

| Method | What to document | Consumed by clean workflow |
| --- | --- | --- |
| scHiCluster | Version/commit, contact-pair input format, imputed matrix output format, command template | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder |
| HiCImpute | Version/commit, input matrix format, imputed output format, command template | Simulation summaries, NPC chrX summaries |
| Higashi | Version/commit, config fields, command template, export/conversion notes | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder |
| Fast-Higashi | Version/commit, config fields, command template, export/conversion notes | Simulation summaries, NPC chrX summaries, SCORE integrated-result JSON folder |
| SCORE | Version/commit, Python environment, `score cooler`, `score embed InnerProduct`, and `score embed SnapATAC` commands | Real-data study 2 |

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

Use `real_data_SCORE_oocyte_zygote/run_SCORE_innerproduct_four_inputs_10runs.R` and `real_data_SCORE_oocyte_zygote/run_SCORE_snapatac_four_inputs_10runs.R` as the exact scripted SCORE calls for the public workflow.
