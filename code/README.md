# Code

Validated analysis scripts grouped by phase. Each script is **self-contained** (no cross-imports)
so it can be run and audited in isolation. Superseded/exploratory scripts (`04`,`05`,`06`,
`mat02`, `mapcheck`, `mapfix`) are intentionally **excluded** — they live in the working archive,
not the paper repo.

| script | produces | paper element |
|---|---|---|
| `1_human_rvlv/01_verify_inputs.py` | input manifest + t-spectrum | Methods §data |
| `1_human_rvlv/02_derive_mapping.py` | listmode→histo mapping (confirmed by `orient.py`) | Methods §mapping |
| `1_human_rvlv/03_extract_heart_region.py` | heart-bbox event cache (17 GB pass) | — |
| `1_human_rvlv/06b_sweep_freeEMG.py` | baseline τ + registration sweep | Results 3.1 |
| `1_human_rvlv/07_plots.py` | fig1–fig4 | Figures |
| `1_human_rvlv/08_gradient_control.py` | intra-chamber gradient | Results 3.x |
| `1_human_rvlv/09_core_vs_rim.py` | blood-core vs wall | Results 3.2 |
| `1_human_rvlv/10_oxygenation_ladder.py` | compartment ladder + matched control | Results 3.3, fig5 |
| `1_human_rvlv/11_ladder_regression.py` | matched-control stats | Results 3.3 |
| `1_human_rvlv/12_lung_air_confound.py` | lung-air check | Results 3.4, fig6 |
| `1_human_rvlv/13_verify_codex_rebuttal.py` | wall-myocardium PSF + HU | Results 3.4 |
| `2_orientation/orient.py` | **chamber identity** (histoimage 0.998 + organ panel) | Methods §identity, Addendum C |
| `3_instrument_82rb/mat01_inspect_82rb.py` | ⁸²Rb format/geometry | Methods §instrument |
| `3_instrument_82rb/mat03_clean_core.py` | IRF + clean-core τ | Results §instrument |
| `3_instrument_82rb/mat04_gradient.py` | positional envelope | Results §instrument |

## Paths
All scripts resolve paths through `code/config.py`.

- `PLI_RAW` points to the downloaded Zenodo records. Default: repo-local `raw/`.
- `PLI_MASKS` points to TotalSegmentator outputs. Default: repo-local `masks/`.
- `PLI_OUT` points to generated tables/caches/figures. Default: repo-local `results/`.

A fresh clone can either export those variables or create the repo-local `raw/` and `masks/`
directories described in `data/DATA.md` and `env/ENVIRONMENT.md`.
