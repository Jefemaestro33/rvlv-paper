#!/usr/bin/env bash
# Driver for the validated pipeline. Activate the venv and set PLI_RAW / PLI_MASKS
# if the data and masks are not under repo-local raw/ and masks/ (see env/ENVIRONMENT.md).
set -euo pipefail
cd "$(dirname "$0")"          # run from code/ regardless of invocation dir
PY="${PY:-python}"

echo "== 1. Human RV/LV =="
$PY 1_human_rvlv/01_verify_inputs.py
$PY 1_human_rvlv/02_derive_mapping.py
$PY 1_human_rvlv/03_extract_heart_region.py     # one 17 GB streaming pass -> cache
$PY 1_human_rvlv/06b_sweep_freeEMG.py           # baseline + registration sweep
$PY 1_human_rvlv/07_plots.py
$PY 1_human_rvlv/08_gradient_control.py
$PY 1_human_rvlv/09_core_vs_rim.py
$PY 1_human_rvlv/10_oxygenation_ladder.py        # full-listmode pass -> oxygenation ladder
$PY 1_human_rvlv/11_ladder_regression.py
$PY 1_human_rvlv/12_lung_air_confound.py
$PY 1_human_rvlv/13_verify_codex_rebuttal.py

echo "== 2. Orientation / chamber identity (tau-free) =="
$PY 2_orientation/orient.py                      # histoimage 0.998 + organ panel

echo "== 3. 82Rb instrument characterization =="
$PY 3_instrument_82rb/mat01_inspect_82rb.py
$PY 3_instrument_82rb/mat03_clean_core.py
$PY 3_instrument_82rb/mat04_gradient.py

echo "== copying regenerated figures to paper/figures/ =="
cp ../results/plots/*.png ../paper/figures/ 2>/dev/null || true

echo "== done. See results/ and paper/figures/. Expected numbers in REPRODUCE.md =="
