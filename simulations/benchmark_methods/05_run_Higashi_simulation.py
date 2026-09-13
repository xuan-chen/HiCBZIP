#!/usr/bin/env python
"""Run Higashi for one prepared simulation chromosome/coverage directory."""

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from higashi.Higashi_wrapper import Higashi


def matrix2d_to_long_offdiag(mtx):
    return np.nan_to_num(mtx[np.tril_indices_from(mtx, k=-1)])


def ensure_label_info_pickle(run_dir):
    out_pkl = run_dir / "label_info.pickle"
    if out_pkl.exists():
        return out_pkl
    df = pd.read_csv(run_dir / "data.txt", sep="\t")
    cells = df[["cell_name", "cell_id"]].drop_duplicates().sort_values("cell_id").reset_index(drop=True)
    cells["label"] = 0
    with out_pkl.open("wb") as handle:
        pickle.dump(cells, handle, protocol=4)
    return out_pkl


def build_runtime_config(run_dir):
    with (run_dir / "config.json").open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    cfg["data_dir"] = str(run_dir.as_posix())
    cfg["temp_dir"] = str((run_dir / "temp").as_posix())
    runtime_config = run_dir / "config.runtime.json"
    with runtime_config.open("w", encoding="utf-8") as handle:
        json.dump(cfg, handle, indent=2)
    return runtime_config, cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_dir", required=True)
    args = parser.parse_args()
    run_dir = Path(args.run_dir).resolve()
    export_dir = run_dir / "export_long"
    export_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "run_metadata.json").open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    ensure_label_info_pickle(run_dir)
    runtime_config, cfg = build_runtime_config(run_dir)
    chrom = metadata["chr_name"]
    n_cells = int(metadata["n_cells"])
    n_bins = int(metadata["n_bins"])
    start_bin = int(metadata["start_bp"]) // int(metadata["resolution"])
    end_bin = start_bin + n_bins

    model = Higashi(str(runtime_config))
    model.process_data()
    model.prep_model()
    model.train_for_embeddings()
    model.train_for_imputation_nbr_0()
    model.impute_no_nbr()
    model.train_for_imputation_with_nbr()
    model.impute_with_nbr()

    vec_len = int(metadata["vec_len"])
    mats = {name: np.zeros((vec_len, n_cells), dtype=float) for name in ("ori", "nbr0", "nbr5")}
    for cell_id in range(n_cells):
        ori, nbr0, nbr5 = model.fetch_map(chrom, cell_id)
        for name, sparse_mtx in zip(("ori", "nbr0", "nbr5"), (ori, nbr0, nbr5)):
            dense = sparse_mtx.toarray()
            region = dense[start_bin:end_bin, start_bin:end_bin]
            mats[name][:, cell_id] = matrix2d_to_long_offdiag(region)

    for name, mat in mats.items():
        np.savetxt(export_dir / f"{name}_{vec_len}x{n_cells}.csv", mat, delimiter=",")
    print(f"Saved Higashi simulation exports under: {export_dir}")


if __name__ == "__main__":
    main()
