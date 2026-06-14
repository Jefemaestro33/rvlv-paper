# Current Status

This repository is a clean, paper-facing package for the single-subject RV/LV
positronium-lifetime re-analysis.

## Ready

- Validated code is grouped under `code/`.
- Public data pointers and integrity checks are in `data/DATA.md`.
- Environment and path configuration are in `env/` and `code/config.py`.
- Paper-ready figures and tables are in `paper/figures/` and `results/`.
- `code/run_all.sh` resolves paths from its own location and can be launched from
  the repository root as `bash code/run_all.sh`.

## Local data layout

By default, scripts look for:

- `raw/zenodo_11243763/`
- `raw/zenodo_12636019/`
- `masks/totalseg_11243763/`
- `masks/totalseg_chambers_11243763/`

Set `PLI_RAW`, `PLI_MASKS`, and optionally `PLI_OUT` to use another location.

## Validation performed

- All Python scripts compile.
- `run_all.sh` syntax and path dry-run pass.
- Smoke tests with explicit `PLI_RAW`, `PLI_MASKS`, and temporary `PLI_OUT` pass:
  `01_verify_inputs.py` reports 529,820,881 human events and
  `mat01_inspect_82rb.py` reports 390,780,205 82Rb material events.
- The orientation result (`histoimage_corr_best = 0.998`) was previously
  reproduced before path parametrization; the path changes do not alter its logic.

## Not included

The raw list-mode files and TotalSegmentator masks are intentionally not tracked.
The full development audit tree is retained outside this repository as provenance.
