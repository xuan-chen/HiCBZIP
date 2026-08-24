# Environment records

This directory contains the minimum environment records included with the public HiCBZIP repository. These files document the R, CmdStan, SCORE, Higashi, and Fast-Higashi environments needed to understand and reproduce the manuscript analyses.

Extra diagnostic files are kept locally in `optional/`, but they are not part of the required public release.

## Included records

| Record group | Purpose | Key confirmed versions |
| --- | --- | --- |
| `local_plotting` | Local R environment used for manuscript summaries and plotting. | See `sessionInfo_local_plotting.txt`. |
| `server_r_env_stan_workflow` | R environment used to run the HiCBZIP Stan workflow. | R 4.5.3; `cmdstanr` 0.8.1. |
| `server_stan_cmdstan_toolchain` | CmdStan compiler/toolchain used by `cmdstanr`. | CmdStan 2.36.0 at `/home/xchen/miniconda3/envs/stan/bin/cmdstan`. |
| `server_score` | SCORE environment used for SCORE metric summaries. | Python 3.10.19; numpy 1.26.4; scipy 1.15.3; scikit-learn 1.7.2; scanpy 1.10.4; cooler 0.9.3. |
| `server_score_higashi_py39_clean` | SCORE/Higashi-compatible Python environment used for integrated benchmark evaluation. | Python 3.9.25; torch 2.8.0+cu128; numpy 1.26.4; scipy 1.13.1; scanpy 1.10.3; cooler 0.9.3. |
| `server_higashi_env_simulation` | Higashi environment used for simulation benchmark runs. | Python 3.9.23; Higashi 0.1.0a0; torch 1.12.1; numpy 1.21.5; scipy 1.7.3. |
| `server_fast_higashi_env_simulation` | Fast-Higashi environment used for simulation benchmark runs. | Python 3.9.25; Fast-Higashi 0.1.1a0; torch 2.8.0+cu128; numpy 1.26.4; scipy 1.13.1. |
| `scl_source_record` | SCL source and executable record for simulation input generation. | SCL1.0 source folder; build command `g++ -o scl scl.cpp`; SHA256 hashes recorded. |

## Files

- `sessionInfo_local_plotting.txt`
- `sessionInfo_server_r_env_stan_workflow.txt`
- `R_version_server_r_env.txt`
- `conda_list_server_r_env.txt`
- `R_version_server_stan_cmdstan_toolchain.txt`
- `conda_list_server_stan_cmdstan_toolchain.txt`
- `conda_list_server_score.txt`
- `conda_list_server_score_higashi_py39_clean.txt`
- `conda_list_server_higashi_env_simulation.txt`
- `conda_list_server_fast_higashi_env_simulation.txt`
- `pip_show_higashi_env_simulation.txt`
- `fasthigashi_metadata_search.txt`
- `scl_source_record.txt`

Some Python packages do not expose a reliable `__version__` attribute. For those cases, package metadata was captured with `pip show` or `importlib.metadata`.
