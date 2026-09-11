#!/usr/bin/env python
"""Run Higashi for one prepared NPC chrX coverage directory."""

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from higashi.Higashi_wrapper import Higashi


def matrix2d_to_long_diag(mtx):
    return np.nan_to_num(mtx[np.tril_indices_from(mtx, k=0)])


def ensure_label_info_pickle(run_dir):
    data_txt = run_dir / "data.txt"
    out_pkl = run_dir / "label_info.pickle"
    if out_pkl.exists():
        return out_pkl

    df = pd.read_csv(data_txt, sep="\t")
    cells = df[["cell_name", "cell_id"]].drop_duplicates().sort_values("cell_id").reset_index(drop=True)
    cells["label"] = 0
    with out_pkl.open("wb") as handle:
        pickle.dump(cells, handle, protocol=4)
    return out_pkl


def patch_runtime_config(run_dir):
    config_path = run_dir / "config.json"
    runtime_config = run_dir / "config.runtime.json"
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)

    data_txt = run_dir / "data.txt"
    df = pd.read_csv(data_txt, sep="\t")
    chrom = cfg["chrom_list"][0]
    resolution = int(cfg["resolution"])
    if df.shape[0] == 0:
        max_pos = 0
        needed_size = 2 * resolution
    else:
        max_pos = int(np.max([df["pos1"].max(), df["pos2"].max()]))
        needed_size = ((max_pos // resolution) + 2) * resolution

    local_sizes = run_dir / "chrom.sizes.local.txt"
    with local_sizes.open("w", encoding="utf-8") as handle:
        handle.write(f"{chrom}\t{needed_size}\n")

    cfg["data_dir"] = str(run_dir.as_posix())
    cfg["temp_dir"] = str((run_dir / "temp").as_posix())
    cfg["genome_reference_path"] = str(local_sizes.as_posix())
    with runtime_config.open("w", encoding="utf-8") as handle:
        json.dump(cfg, handle, indent=2)

    print(f"Runtime config: {runtime_config}")
    print(f"Local chrom size for {chrom}: {needed_size} bp (max observed position {max_pos})")
    return runtime_config, cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_dir", required=True, help="Prepared directory with config.json and data.txt")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    export_dir = run_dir / "export_long_diag"
    export_dir.mkdir(parents=True, exist_ok=True)
    ensure_label_info_pickle(run_dir)
    runtime_config, cfg = patch_runtime_config(run_dir)
    chrom = cfg["chrom_list"][0]
    n_cells = int(json.load((run_dir / "run_metadata.json").open())["n_cells"])

    model = Higashi(str(runtime_config))
    model.process_data()
    model.prep_model()
    model.train_for_embeddings()
    model.train_for_imputation_nbr_0()
    model.impute_no_nbr()
    model.train_for_imputation_with_nbr()
    model.impute_with_nbr()

    ori0, _, _ = model.fetch_map(chrom, 0)
    vec_len = len(matrix2d_to_long_diag(ori0.toarray()))
    ori_mat = np.zeros((vec_len, n_cells), dtype=float)
    nbr0_mat = np.zeros((vec_len, n_cells), dtype=float)
    nbr5_mat = np.zeros((vec_len, n_cells), dtype=float)

    for cell_id in range(n_cells):
        ori, nbr0, nbr5 = model.fetch_map(chrom, cell_id)
        ori_mat[:, cell_id] = matrix2d_to_long_diag(ori.toarray())
        nbr0_mat[:, cell_id] = matrix2d_to_long_diag(nbr0.toarray())
        nbr5_mat[:, cell_id] = matrix2d_to_long_diag(nbr5.toarray())

    np.savetxt(export_dir / f"ori_diag_{vec_len}x{n_cells}.csv", ori_mat, delimiter=",")
    np.savetxt(export_dir / f"nbr0_diag_{vec_len}x{n_cells}.csv", nbr0_mat, delimiter=",")
    np.savetxt(export_dir / f"nbr5_diag_{vec_len}x{n_cells}.csv", nbr5_mat, delimiter=",")
    print(f"Saved Higashi exports under: {export_dir}")


if __name__ == "__main__":
    main()
