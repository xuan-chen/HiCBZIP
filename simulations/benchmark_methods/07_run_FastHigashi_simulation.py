#!/usr/bin/env python
"""Run Fast-Higashi for one prepared simulation chromosome/coverage directory."""

import argparse
import json
import pickle
from pathlib import Path

import pandas as pd
from fasthigashi.FastHigashi_Wrapper import FastHigashi


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
    cfg.setdefault("input_format", "higashi_v1")
    cfg.setdefault("structured", True)
    cfg.setdefault("resolution_fh", [cfg["resolution"]])
    runtime_config = run_dir / "config.runtime.json"
    with runtime_config.open("w", encoding="utf-8") as handle:
        json.dump(cfg, handle, indent=2)
    return runtime_config, cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--rank", type=int, default=3)
    parser.add_argument("--dim1", type=float, default=0.6)
    parser.add_argument("--tol", type=float, default=2e-5)
    parser.add_argument("--off_diag", type=int, default=100)
    parser.add_argument("--conv_threshold", type=float, default=0.1)
    parser.add_argument("--skip_partial_rwr", action="store_true")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    with (run_dir / "run_metadata.json").open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    coverage = float(metadata["coverage"])
    do_conv = coverage < args.conv_threshold
    ensure_label_info_pickle(run_dir)
    runtime_config, cfg = build_runtime_config(run_dir)

    settings = {
        "coverage": metadata["coverage"],
        "rank": args.rank,
        "off_diag": args.off_diag,
        "conv_threshold": args.conv_threshold,
        "do_rwr": True,
        "do_conv": do_conv,
        "diagonal_included_export": False,
    }
    with (run_dir / "run_runtime_settings.json").open("w", encoding="utf-8") as handle:
        json.dump(settings, handle, indent=2)

    wrapper = FastHigashi(
        config_path=str(runtime_config),
        path2input_cache=str(run_dir),
        path2result_dir=str(run_dir),
        off_diag=args.off_diag,
        filter=False,
        do_conv=do_conv,
        do_rwr=True,
        do_col=False,
        no_col=False,
    )

    raw_cache = run_dir / "temp" / "raw" / f"{cfg['chrom_list'][0]}_sparse_adj.npy"
    if not raw_cache.exists():
        wrapper.fast_process_data()
    wrapper.prep_dataset(batch_norm=True)
    wrapper.run_model(dim1=args.dim1, rank=args.rank, n_iter_parafac=1, tol=args.tol, extra=f"{run_dir.name}_rank{args.rank}")
    if not args.skip_partial_rwr:
        wrapper.only_partial_rwr()
    print("Fast-Higashi simulation run completed.")


if __name__ == "__main__":
    main()
