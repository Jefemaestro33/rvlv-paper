# Paper skeleton — confound-controlled test of positronium lifetime as an in-vivo oxygenation biomarker

> Status: scaffold for a solo, no-IRB preprint built on public data. Honest framing:
> the claim is the **falsification + the confound-control method**, NOT a positive mechanism.

## Title (options, honest — emphasize the null + the method)
1. *An intracardiac positronium-lifetime contrast in human PET is not explained by blood oxygenation: a confound-controlled re-analysis*
2. *A matched-vessel control argues against positronium lifetime as an in-vivo oxygenation biomarker*

## Positioning / authorship
- Solo-able (Darell, lead). Public data → **no IRB**. Aligns with the program strategy:
  first win = a non-IRB validation/skeptic artifact (preprint), biology↔physics bridge.
- Respectful re-analysis of Mercolli et al. (Bern/Inselspital) public 82Rb data. Credit
  them for demonstrating in-vivo o-Ps *measurability*; address the **field's** biomarker
  aspiration (J-PET/Moskal hypoxia/tissue-marker program), not their claims.
- Methodological strength to state plainly: results independently re-derived by **two
  separate analysis pipelines** that converged (and one caught an over-claim in the other —
  worth a sentence as a reproducibility safeguard).

## Target venues
- Preprint: **arXiv physics.med-ph** (+ medRxiv mirror).
- Journal: **EJNMMI Physics** (same community/data), **Phys. Med. Biol.**, or **IEEE TRPMS**.

## Abstract (draft, ~200 words)
> Ortho-positronium (o-Ps) lifetime is increasingly proposed as an in-vivo biomarker of
> tissue oxygenation/hypoxia. We test this directly using the natural arterial/venous
> oxygenation contrast of the human heart, re-analysing a public long-axial-FOV PET/CT
> [⁸²Rb] dataset (5.3×10⁸ list-mode events). We independently derive the list-mode↔CT
> coordinate mapping, extract per-chamber o-Ps lifetimes (two-component EMG model), and
> find a right-ventricle > left-ventricle contrast (Δτ ≈ +0.30 ns) in the oxygen-expected
> direction. The contrast is stable to ±16 mm registration perturbation and resides in the
> blood core rather than the chamber wall. However, a structure- and position-matched
> control — the pulmonary artery (deoxygenated) versus the aorta (oxygenated) — shows no
> detectable difference (Δτ = −0.01 ± 0.34 ns), ruling out an oxygenation effect of the
> observed (~0.3 ns) magnitude though not a small physiological one; lung-air partial volume,
> wall-myocardium leakage and registration are disfavoured as drivers, and the effect falls
> within an isotope-matched ⁸²Rb uniform-quartz instrumental positional envelope (~0.1–0.35 ns
> over 40 mm). The residual right-heart elevation (~1.5σ) is
> not explained by oxygenation and its cause is unresolved. We conclude that, on this dataset,
> in-vivo o-Ps lifetime does not yield a clean oxygenation readout, and we provide a reusable
> confound-control battery (matched-vessel test + exclusion checks) for evaluating future o-Ps
> oxygenation claims. Multi-subject replication with anatomy decoupled from oxygenation is required.

## 1. Introduction
- o-Ps lifetime physics in matter (pick-off + paramagnetic O₂ ortho→para conversion).
- The in-vivo biomarker aspiration (hypoxia/tissue characterisation); why it matters clinically.
- The identifiability problem: τ is dominated by tissue structure/free volume; the O₂ term is
  small (tens of ps for physiological dissolved-O₂ differences) and confounded.
- The heart as a natural in-vivo oxygenation test (RV deoxy vs LV oxy), and its built-in
  control: pulmonary artery (deoxy) vs pulmonary vein/aorta (oxy) decouple O₂ from vessel type.
- Contribution: a confound-controlled test on public data; a cautionary null + a reusable method.

## 2. Methods
- **Data:** Zenodo 11243763, [⁸²Rb]Cl human, Biograph Vision Quadra; 529,820,881 events
  (x,y,z,t float64); AC-CT + reconstructed histo-image.
- **Coordinate mapping & chamber identity (derived + rigorously confirmed, τ-free):** 48-way
  axis-permutation/sign search; the chosen mapping (−LMx,+LMz,−LMy) **reproduces the provider's
  histoimage at corr 0.998** (all 47 alternatives ≤0.60; coarse-block, flip-sensitive) and wins a
  **joint multi-organ panel** (heart+liver+spleen+both kidneys; a single organ is not decisive) —
  confirming RV/LV identity with **no flip, zero lifetime information**. CT/histo/masks LPS-consistent.
  (Footnote: differs from the documented Zenodo `read_singles_binary`, which likely targets a
  CT-aligned frame; histoimage 0.998 is decisive.)
- **Compartment masks:** TotalSegmentator (heart-chambers + total tasks): RV, LV, RA, LA,
  aorta, pulmonary artery, pulmonary vein, SVC, IVC, myocardium, lung lobes.
- **Lifetime model:** background + 2-component EMG (Gaussian IRF ⊗ short + o-Ps exponentials),
  fit over [−3,12] ns, Poisson-weighted; tail single-exponential as model-independent check.
- **Registration stress test:** global translational sweep (±3 vox dense grid, ±10 vox per-axis).
- **Confound battery:** core-vs-rim erosion; lung-air PSF field; myocardium PSF field; CT HU;
  **isotope-matched ⁸²Rb reference** (uniform quartz, Zenodo 12636019) for the instrument
  positional envelope + IRF validation.
- **Reproducibility:** two independent pipelines; scripts + intermediate products archived.

## 3. Results
- **3.0 Chamber identity (τ-free).** The mapping reproduces the provider's histoimage at
  **corr 0.998** (vs ≤0.60 for all 47 alternative orientations) and wins a joint
  heart+liver+spleen+kidney organ panel → **RV/LV identity confirmed, no flip**, using no
  lifetime information (independently reproduced by the second pipeline).
- **3.1 Baseline + registration stability.** RV 1.700±0.172, LV 1.396±0.130 ns; Δ +0.304
  (1.4σ); 100% sign-stable across 342 shifts (σ 0.09 ns); per-axis positive to ±16 mm.
  → *Fig 1 (spectra), Fig 2 (grid), Fig 3 (per-axis).*
- **3.2 Blood not wall.** Core (erode-2) Δ +0.43; rim Δ +0.03 (RV 1.54 ≈ LV 1.51).
- **3.3 Oxygenation ladder + matched control.** Table of compartments; **pulm. artery − aorta
  = −0.01 ± 0.34 ns (no detectable difference; ±0.34 ns ≫ expected physiological O₂ effect)**;
  venous>arterial trend carried by right-heart chambers. → *Fig 5 (ladder + τ-vs-SO₂).*
- **3.4 Confound checks.** Lung-air (corr −0.04, air≈0 in cores → disfavoured); wall-myocardium
  external leakage (LV more myo-PSF, lower τ; corr +0.06 → refuted); HU non-separability (86%
  overlap → uninformative); **isotope-matched ⁸²Rb uniform-quartz envelope** — positional τ
  varies ~0.1–0.35 ns over 40 mm (systematic bound, *not* a calibrated gradient; IRF validated
  t0 0.067≈human), and the RV-LV effect falls within it. → *Fig 6 (lung-air), confound-battery table.*
- **3.5 Open mechanism.** Right-heart elevation real (~1.5σ) but cause unresolved; candidate
  hypotheses (internal trabeculae + ⁸²Rb avidity; regional scatter; bolus kinetics) — none proven.

## 4. Discussion
- Why the magnitude (0.3–0.4 ns) is ~10× the dissolved-O₂ prediction → it was never O₂.
- Oxygenation/structure degeneracy in cardiac anatomy; the matched control as the way to break it.
- Deoxyhemoglobin (paramagnetic) would push the *opposite* way → not dominant either.
- Implications for the field's o-Ps hypoxia-imaging program: natural in-vivo contrasts are
  confounded; quantitative O₂ claims need matched controls + much higher effective statistics.
- The reusable confound battery as the portable contribution.

## 5. Limitations (state plainly)
- **Single subject; contrasts ~1.5σ** — not a detection of RV>LV. The robust results are the
  *no-detectable-difference* of the matched control and the disfavouring of named confounds —
  not a positive claim of zero oxygenation effect.
- **Precision floor:** the matched-control uncertainty (±0.34 ns) exceeds the expected
  physiological dissolved-O₂ effect (tens of ps), so a *small* O₂ contribution cannot be excluded.
- **Instrumental:** an isotope-matched ⁸²Rb uniform-quartz reference shows a positional τ
  envelope ~0.1–0.35 ns over 40 mm (systematic bound, **not** a calibrated gradient — direction
  unstable across cores); the RV-LV effect falls within it. IRF validated (t0 0.067≈human 0.066).
  Patient-specific scatter/attenuation/randoms/motion still not individually modelled.
- o-Ps effective statistics modest (small fraction on a large flat background); 2-EMG model
  systematic ~0.1 ns; χ²/ndf 2–3.
- TotalSegmentator masks on CT; thin vessels (carotid/subclavian/pulm. vein) too noisy to use.
- Mechanism of right-heart elevation unresolved (internal trabeculae / scatter / bolus untested).

## 6. Conclusion
- In-vivo o-Ps lifetime is not a clean oxygenation biomarker on this dataset; apparent cardiac
  contrasts are anatomy-confounded. The matched-vessel control + exclusion battery should be
  standard before any o-Ps oxygenation claim.

## Figures (in `paper/figures/`)
- Fig 1 `fig1_spectra.png` — RV/LV Δt spectra + EMG fits.
- Fig 2 `fig2_grid_distribution.png` — RV−LV over ±3-vox registration grid.
- Fig 3 `fig3_per_axis.png` — RV−LV vs single-axis shift.
- Fig 4 `fig4_vs_radius.png` — RV−LV vs |shift|.
- Fig 5 `fig5_oxygenation_ladder.png` — ladder + τ-vs-SO₂ (the key null).
- Fig 6 `fig6_lung_air_confound.png` — τ vs lung-air fraction.
- (Optional) Fig 7 — confound-battery summary schematic.

## Data & code availability
- Public data: Zenodo 11243763 (and 13443797 phantom, 12636019 reference materials).
- Analysis scripts and paper-ready tables/figures: this repository. The full development audit tree
  is retained separately as provenance and can be archived as supplementary material if needed.

## What's needed before submission (gaps / decisions)
1. **Framing decision:** ship as single-subject *caution* now, OR upgrade to multi-subject.
2. **Multi-subject (the strengthening lever):** request the Ga68-DOTATOC / Ga68-PSMA listmode
   from Bern (Mercolli/Steinberger). Whole-body → aorta + vena cava present → replicate the
   arterial/venous + matched-vessel test in independent anatomy/positioning. Decisive upgrade.
3. **Pulmonary-vein tie-breaker:** still noisy (±0.78); a better PV measurement would add a
   second crossover control.
4. Prose, references (Mercolli 2024/2025; Moskal positronium imaging; Shibuya O₂ quenching;
   SIMPLE-Moment recon), author/affiliation, figure polish.
