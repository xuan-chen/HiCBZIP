#!/usr/bin/env python
"""Export Fast-Higashi simulation HDF5 output as an off-diagonal long-vector CSV."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def matrix2d_to_long_offdiag(mtx):
    return np.nan_to_num(mtx[np.tril_indices_from(mtx, k=-1)])


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

    export_dir = run_dir / "export_long"
    export_dir.mkdir(parents=True, exist_ok=True)
    chrom = metadata["chr_name"]
    n_cells = int(metadata["n_cells"])
    vec_len = int(metadata["vec_len"])
    n_bins = int(metadata["n_bins"])
    start_bin = int(metadata["start_bp"]) // int(metadata["resolution"])
    end_bin = start_bin + n_bins

    vectors = []
    with h5py.File(h5_path, "r") as handle:
        if chrom not in handle:
            raise KeyError(f"Chromosome group {chrom} not found in {h5_path}")
        group = handle[chrom]
        cell_keys = sorted((x for x in group.keys() if x != "shape"), key=lambda x: int(x))
        if len(cell_keys) != n_cells:
            raise ValueError(f"Expected {n_cells} cells, found {len(cell_keys)}")
        for key in cell_keys:
            dense = np.asarray(group[key], dtype=float)
            if dense.shape == (n_bins, n_bins):
                region = dense
            else:
                region = dense[start_bin:end_bin, start_bin:end_bin]
            vectors.append(matrix2d_to_long_offdiag(region))

    mat = np.column_stack(vectors)
    if mat.shape != (vec_len, n_cells):
        raise ValueError(f"Exported shape {mat.shape} does not match expected {(vec_len, n_cells)}")

    out_csv = export_dir / f"fasthigashi_prwr_{vec_len}x{n_cells}.csv"
    np.savetxt(out_csv, mat, delimiter=",")
    with (export_dir / f"fasthigashi_prwr_{vec_len}x{n_cells}_qc.json").open("w", encoding="utf-8") as handle:
        json.dump({"run_dir": str(run_dir.as_posix()), "coverage": metadata["coverage"], "final_shape": [vec_len, n_cells]}, handle, indent=2)
    print(f"Saved Fast-Higashi export: {out_csv}")


if __name__ == "__main__":
    main()
