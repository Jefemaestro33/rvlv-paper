# Environment

Pure-Python, CPU-only (no GPU needed). Developed on Ubuntu 22.04, Python 3.10.12.

```bash
python3.10 -m venv .venv
source .venv/bin/activate
pip install -r env/requirements.txt
```

## Pinned versions (as used)
- python 3.10.12
- numpy 2.2.6
- scipy 1.15.3
- nibabel 5.4.2
- matplotlib 3.10.9

## Paths
By default the scripts look for repo-local `raw/` and `masks/` directories and write to
repo-local `results/`. You can override all three locations:
```bash
export PLI_RAW=/path/to/raw          # contains zenodo_11243763/ and zenodo_12636019/
export PLI_MASKS=/path/to/masks      # contains totalseg_11243763/ and totalseg_chambers_11243763/
export PLI_OUT=/path/to/results      # optional; defaults to <repo>/results
```
Outputs are written under `results/` and `paper/figures/`.

On the project VM, the validated paths are:
```bash
export PLI_RAW=/mnt/pli_scratch/active_raw/bern
export PLI_MASKS=/mnt/pli_scratch/work
```

## Compute notes
- The two list-mode files are ~17 GB and ~12.5 GB; scripts stream them via `np.memmap` and
  sample/chunk — peak RAM stays well under 32 GB.
- A full human pipeline run is a few minutes (one 17 GB streaming pass is cached to
  `results/heart_region_events.npy`; subsequent steps reuse it).
