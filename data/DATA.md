# Data availability

All inputs are public, CC-BY-4.0, deposited by Mercolli/Steinberger (Dept. of Nuclear Medicine,
Inselspital, Bern University Hospital) for positronium-lifetime studies on a Biograph Vision
Quadra. **The raw list-mode is not stored in this repository** (≈30 GB total); download it from
Zenodo into a local `raw/` directory and point the scripts at it (see `env/ENVIRONMENT.md`).

## Records used
| Zenodo | role | size | key file |
|---|---|---|---|
| [11243763](https://doi.org/10.5281/zenodo.11243763) | **in-vivo human [⁸²Rb]Cl** (cardiac) | 17.2 GB | `20230606_positronium_patient_evaluated_data_Histo_Out.l` + AC-CT + histoimage |
| [12636019](https://doi.org/10.5281/zenodo.12636019) | **⁸²Rb reference material** (uniform quartz disks) | 23.4 GB | `Rb82_coins_Histo_Out.l` |
| [13443797](https://doi.org/10.5281/zenodo.13443797) | ¹²⁴I phantom (context only) | 3.3 GB | `Histo_Out_separated.l` |

## Download
```bash
# human subject (only the listmode + CT + histoimage are needed)
wget -c "https://zenodo.org/records/11243763/files/20230606_positronium_patient_evaluated_data_Histo_Out.l?download=1" -O raw/human/Histo_Out.l
wget -c "https://zenodo.org/records/11243763/files/20230606_lm_2_AC_CT_WB_1.65mm.nii.gz?download=1" -O raw/human/CT.nii.gz
wget -c "https://zenodo.org/records/11243763/files/20230606_histoimage3d_lm.nii.gz?download=1" -O raw/human/histoimage.nii.gz
# 82Rb reference material (instrument characterization)
wget -c "https://zenodo.org/records/12636019/files/Rb82_coins_Histo_Out.l?download=1" -O raw/material/Rb82_coins_Histo_Out.l
```

## Integrity (verify after download)
| file | bytes | events (bytes/32) |
|---|---|---|
| human `Histo_Out.l` | 16,954,268,192 | 529,820,881 |
| material `Rb82_coins_Histo_Out.l` | 12,504,966,560 | 390,780,205 |

`size % 32 == 0` (32 bytes/event = 4× float64). Human listmode sha256(head‖tail) =
`bacde203a82e224fde372e82ae587436a6eca3b779c9883c7c44d22071163e7b`.

## List-mode format & coordinate convention (important)
Each event is `float64 (x, y, z, t)`, 32 bytes; positions in mm (scanner frame), `t` = the
prompt→annihilation time difference in ns. Sentinel/invalid `t` values appear at ±655360 ns
(drop with `|t| < 50`).

**Coordinate mapping** (resolved empirically; see `code/2_orientation/`): the mapping that
reproduces the provider's `histoimage` (the authors' own reconstruction) at coarse-block
correlation **0.998** — and that wins a joint multi-organ alignment panel — is
`world = (−x, +z, −y)` of the stored columns. This **differs** from the documented Zenodo
`read_singles_binary` helper (`(+x, −z, −y)`, which appears to target a CT-aligned frame); the
0.998 histoimage reproduction is the decisive ground truth and fixes RV/LV identity with no flip.

## Anatomical masks
Organ/chamber masks are generated from the human CT with **TotalSegmentator** (`total` task →
organs incl. liver/spleen/kidneys; `heartchambers_highres` task → RV, LV, RA, LA, aorta,
pulmonary artery, myocardium). The `heartchambers` task needs a free academic license. See
`code/1_human_rvlv/` for the exact mask paths expected.
