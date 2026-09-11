#!/usr/bin/env python
"""Export Fast-Higashi NPC chrX imputed maps as diagonal-included long vectors."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def matrix2d_to_long_diag(mtx):
    return np.nan_to_num(mtx[np.tril_indices_from(mtx, k=0)])


def valid_bins_from_data(run_dir, metadata):
    df = pd.read_csv(run_dir / "data.txt", sep="\t")
    n_bins = int(metadata["n_bins"])
    resolution = int(metadata["resolution"])
    valid = np.zeros(n_bins, dtype=bool)
    if df.shape[0] == 0:
        return valid
    bin1 = df["pos1"].to_numpy(dtype=int) // resolution
    bin2 = df["pos2"].to_numpy(dtype=int) // resolution
    valid[bin1[(bin1 >= 0) & (bin1 < n_bins)]] = True
    valid[bin2[(bin2 >= 0) & (bin2 < n_bins)]] = True
    return valid


def expand_if_needed(mtx, metadata, valid_bins):
    n_bins = int(metadata["n_bins"])
    if mtx.shape == (n_bins, n_bins):
        return mtx
    kept = np.flatnonzero(valid_bins)
    if mtx.shape != (len(kept), len(kept)):
        raise ValueError(f"Matrix shape {mtx.shape} does not match full bins or kept-bin subset.")
    full = np.zeros((n_bins, n_bins), dtype=float)
    full[np.ix_(kept, kept)] = mtx
    return full


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_dir", required=True)
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    with (run_dir / "run_metadata.json").open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    h5_path = run_dir / "impute_prwr.hdf5"
    if not h5_path.exists():
        raise FileNotFoundError(f"Missing Fast-Higashi output: {h5_path}")

    export_dir = run_dir / "export_long_diag"
    export_dir.mkdir(parents=True, exist_ok=True)
    valid_bins = valid_bins_from_data(run_dir, metadata)
    chrom = metadata["chr_name"]

    vectors = []
    with h5py.File(h5_path, "r") as handle:
        if chrom not in handle:
            raise KeyError(f"Chromosome group {chrom} not found in {h5_path}")
        group = handle[chrom]
        cell_keys = sorted((x for x in group.keys() if x != "shape"), key=lambda x: int(x))
        for key in cell_keys:
            mtx = np.asarray(group[key], dtype=float)
            full_mtx = expand_if_needed(mtx, metadata, valid_bins)
            vectors.append(matrix2d_to_long_diag(full_mtx))

    mat = np.column_stack(vectors)
    expected = (int(metadata["vec_len"]), int(metadata["n_cells"]))
    if mat.shape != expected:
        raise ValueError(f"Exported matrix shape {mat.shape} does not match expected {expected}")

    out_csv = export_dir / f"fasthigashi_prwr_diag_{expected[0]}x{expected[1]}.csv"
    np.savetxt(out_csv, mat, delimiter=",")
    qc = {
        "run_dir": str(run_dir.as_posix()),
        "h5_path": str(h5_path.as_posix()),
        "coverage": metadata["coverage"],
        "final_shape": [int(mat.shape[0]), int(mat.shape[1])],
        "diagonal_included": True,
    }
    with (export_dir / f"fasthigashi_prwr_diag_{expected[0]}x{expected[1]}_qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2)
    print(f"Saved Fast-Higashi export: {out_csv}")


if __name__ == "__main__":
    main()
