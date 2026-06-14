# Audit trail (independent dual-pipeline cross-check)

Every result here was produced by **two independent analysis pipelines** (developed separately
and run against the same public data) that converged. The second pipeline acted adversarially —
it twice caught the first over-claiming, forcing the rigorous tests that the paper now rests on.
We record this because the corrections *are* the methodology.

## Corrections forced during development
1. **"Mechanism = internal trabecular myocardium" → retracted.** An interim reading argued the
   right-heart elevation was muscle, because RV core τ ≈ myocardium τ. A follow-up refuted the
   *wall*-myocardium version (a myocardium-PSF field shows LV has 3× more nearby myocardium yet
   *lower* τ) and showed CT Hounsfield units cannot separate blood from muscle (86% overlap). The
   positive mechanism is left **open**; only what it is *not* is claimed.
2. **"Chamber identity proved by the liver anchor + canonical τ" → replaced.** Liver enrichment
   alone is not decisive (free translation lets several orientations hit the liver), and the
   canonical-registration τ had LV pinned at a fit bound. Identity was re-established the right
   way — **τ-free**: the mapping reproduces the provider's histoimage at corr 0.998 (vs ≤0.60 for
   all 47 alternatives) and wins a *joint* multi-organ panel (heart+liver+spleen+kidneys). No flip.

## What is robust vs weak (stated plainly)
- **Robust:** the matched-vessel *null* (−0.01 ± 0.34 ns); registration stability; blood-not-wall;
  chamber identity (histoimage 0.998); the isotope-matched instrument positional envelope; IRF
  validation (t0 0.067 ≈ human 0.066).
- **Weak / bounded:** the RV−LV effect itself is ~1.5σ on one subject; absolute τ is
  model-dependent (~0.1 ns; our quartz 1.71 vs the reference 1.589); the instrument envelope is a
  *systematic bound* (~0.1–0.35 ns), not a calibrated gradient (direction unstable across cores).

## Next (to strengthen beyond a single-subject caution)
- **Multi-subject replication** — the only published in-vivo o-Ps list-mode is this one ⁸²Rb
  subject; two further human subjects (Ga68) and a higher-quality ⁴⁴Sc acquisition exist but are
  not yet public (requestable from the Bern authors).
- **Ground-truth oxygenation** — no public dataset pairs in-vivo o-Ps with independent pO₂; the
  cheapest unlock is bench-top PALS (²²Na + scintillators) on O₂-controlled samples, not more PET.

*(The full working archive, including the exploratory and superseded scripts and the per-step
re-run logs from the second pipeline, is kept separately from this clean paper repository.)*
