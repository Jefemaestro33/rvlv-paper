# Cardiac positronium-lifetime imaging is not a clean in-vivo oxygenation biomarker
### A confound-controlled re-analysis of public long-axial-FOV PET/CT data

> **Status: medRxiv preprint — single subject.** This repository reproduces a cautionary
> methodological result on one public human dataset. DOI:
> https://doi.org/10.64898/2026.06.14.26355630. The conclusion is a *null with controls*,
> not a positive claim; multi-subject replication is the natural next step (see `AUDIT.md` §next).

## TL;DR
Ortho-positronium (o-Ps) lifetime is increasingly proposed as an in-vivo biomarker of tissue
oxygenation. Using the natural arterial/venous oxygenation contrast of the human heart in a
public long-axial-FOV PET/CT [⁸²Rb] dataset, we find a reproducible right-ventricle >
left-ventricle o-Ps lifetime contrast (Δτ ≈ +0.30 ns) **in the oxygen-expected direction** —
but it is **not attributable to oxygenation**:

- a **structure- and position-matched control** (pulmonary artery, deoxygenated, vs aorta,
  oxygenated) shows **no detectable difference** (−0.01 ± 0.34 ns);
- the contrast is **registration-stable** and lives in the **blood core** (not the wall), and the
  **chamber identity is rigorously confirmed** (the mapping reproduces the provider's histoimage
  at corr 0.998; no flip);
- it falls **within the instrument's own positional envelope** measured in an **isotope-matched
  ⁸²Rb uniform-quartz reference** (τ varies ~0.1–0.35 ns over ~40 mm with no tissue difference);
- lung-air and wall-myocardium partial-volume are individually disfavoured.

At the achievable precision (~1.5σ on one subject) the cardiac o-Ps contrast does not provide a
clean oxygenation readout. We contribute a **reusable confound-control battery** (matched-vessel
test + exclusion checks + isotope-matched instrument reference) for evaluating any future o-Ps
oxygenation claim.

## Repository map
| path | what |
|---|---|
| `paper/manuscript.md` | the manuscript draft |
| `paper/figures/` | publication figures (fig1–fig6) |
| `code/1_human_rvlv/` | human RV/LV pipeline (mapping → fits → registration sweep → oxygenation ladder → confounds) |
| `code/2_orientation/` | coordinate-convention & RV/LV chamber-identity determination (τ-free) |
| `code/3_instrument_82rb/` | ⁸²Rb reference-material instrument characterization (IRF + positional envelope) |
| `code/run_all.sh` | driver: runs the validated pipeline end-to-end |
| `results/` | tables the paper cites (`*.csv`) + `key_numbers.json` |
| `data/DATA.md` | Zenodo pointers, download commands, checksums, coordinate-convention note |
| `env/` | `requirements.txt` + setup |
| `REPRODUCE.md` | step-by-step reproduction with expected numbers |
| `AUDIT.md` | the independent dual-audit trail (Claude ↔ Codex) and the corrections it forced |

## Data
All inputs are public (CC-BY-4.0), Bern/Inselspital on Zenodo. The raw list-mode (≈17 GB human,
≈12.5 GB ⁸²Rb material) is **not** included here — see `data/DATA.md` to download. No new data
were acquired; no IRB.

## Reproduce
See `REPRODUCE.md`. In short: set up `env/`, download the two Zenodo records, run `code/run_all.sh`.
Expected headline numbers are listed per step for diffing.

## Reproducibility / integrity
Every result was produced by **two independent analysis pipelines** that converged; the second
pipeline (Codex) caught and corrected two over-claims during development (see `AUDIT.md`). This
adversarial cross-check is part of the method, not an afterthought.

## License & citation
Code: MIT. Text & figures: CC-BY-4.0. Citation: medRxiv preprint,
doi: https://doi.org/10.64898/2026.06.14.26355630.
Built on the public data of Mercolli, Steinberger et al. (Inselspital Bern) — please cite their
original deposits (Zenodo 11243763, 12636019; EJNMMI Phys 2024/2025) when using this work.
