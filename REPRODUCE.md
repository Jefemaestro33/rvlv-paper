# Reproduce

Set up `env/`, download the two Zenodo records (`data/DATA.md`), then run the validated pipeline.
Each step lists the **expected headline numbers** so you can diff an independent run.

```bash
source .venv/bin/activate
export PLI_RAW=/path/to/raw          # contains zenodo_11243763/ and zenodo_12636019/
export PLI_MASKS=/path/to/masks      # contains totalseg_11243763/ and totalseg_chambers_11243763/
bash code/run_all.sh          # or run the steps below individually
```

On the project VM, use:
```bash
export PLI_RAW=/mnt/pli_scratch/active_raw/bern
export PLI_MASKS=/mnt/pli_scratch/work
```

## 1 — Human RV/LV (`code/1_human_rvlv/`)
| step | script | expected key output |
|---|---|---|
| inputs | `01_verify_inputs.py` | 529,820,881 events; prompt peak ~t0; mask voxel counts |
| mapping | `02_derive_mapping.py` | mapping `(−x,+z,−y)`; 16× heart enrichment |
| extract (17 GB pass) | `03_extract_heart_region.py` | ~41.8M heart-bbox events cached |
| baseline + registration sweep | `06b_sweep_freeEMG.py` | RV 1.700±0.172, LV 1.396±0.130, **RV−LV +0.304 (1.4σ)**; 100% sign-stable over ±16 mm, σ 0.09 ns |
| figures | `07_plots.py` | fig1–fig4 |
| intra-chamber gradient | `08_gradient_control.py` | within-LV swing +0.43 ns (cautionary) |
| blood vs wall | `09_core_vs_rim.py` | core Δ +0.43, rim Δ +0.03 (in blood, not wall) |
| oxygenation ladder | `10_oxygenation_ladder.py` | **pulm. artery − aorta = −0.01 ± 0.34 ns (null)**; fig5 |
| ladder regression | `11_ladder_regression.py` | matched-control summary |
| lung-air confound | `12_lung_air_confound.py` | corr(τ,air) −0.04; fig6 |
| wall-myocardium check | `13_verify_codex_rebuttal.py` | LV more myo-PSF yet lower τ → wall leakage refuted |

## 2 — Orientation / chamber identity (`code/2_orientation/`, τ-free)
| script | expected key output |
|---|---|
| `orient.py` | **histoimage reproduction: mapping `(−x,+z,−y)` → corr 0.998; all 47 alternatives ≤ 0.60.** Joint organ panel (heart+liver+spleen+kidneys) won by the same mapping → **RV is the right ventricle, no flip** |

## 3 — ⁸²Rb instrument characterization (`code/3_instrument_82rb/`)
| script | expected key output |
|---|---|
| `mat01_inspect_82rb.py` | 390,780,205 events; compact uniform-quartz source |
| `mat03_clean_core.py` | **IRF: t0 = 0.0672, σ = 0.0879** (≈ human 0.066/0.108); τ_oPs(quartz) 1.713 ± 0.022 ns (paper Bayesian 1.589 — model/VOI); within-core spread 0.118 ns |
| `mat04_gradient.py` | **positional τ envelope ~0.19 (±25 mm core) to ~0.35 ns (±15 mm) over 43 mm** — order-of-magnitude *systematic bound*, not a calibrated gradient |

## Headline numbers (also in `results/key_numbers.json`)
- RV−LV (human, baseline) = **+0.304 ns (~1.4–1.5σ)**, blood-borne, registration-stable, chamber-identity-confirmed.
- Matched-vessel control (pulm. artery vs aorta) = **−0.01 ± 0.34 ns** → rules out O₂ of the RV/LV *magnitude*; small physiological O₂ (≪0.3 ns) not excluded.
- Instrument positional envelope (uniform ⁸²Rb quartz) = **~0.1–0.35 ns over 40 mm** → the effect falls within it.
- → **not a clean oxygenation biomarker** on this subject.

> Superseded/exploratory scripts (`04`,`05`,`06`, `mat02`, `mapcheck`, `mapfix`) are **not** in this
> repo; they live in the working archive. `02_derive_mapping.py`'s mapping is *confirmed* by
> `orient.py` (the histoimage test), which is the citable basis for chamber identity.
