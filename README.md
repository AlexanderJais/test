# Underwater microphone data vs. electron plasma oscillations

An evidence-based investigation of whether electron plasma (Langmuir)
oscillations can be identified and located in publicly available underwater
microphone (hydrophone) data.

**Read the findings: [REPORT.md](REPORT.md)** (part 1 — can plasma oscillations appear in hydrophone data?)
and **[PLASMA_CRAFT_REPORT.md](PLASMA_CRAFT_REPORT.md)** (part 2 — can public hydrophone *arrays* find and
locate a plasma-propelled underwater craft? Includes a TDOA localization pipeline validated to 25 km
against a USGS ground-truth event).

Short version: public hydrophone data exists and is excellent (we analyze
MBARI's open 256 kHz MARS archive directly from AWS); electron plasma
oscillations cannot appear in it for three independent physical reasons
(wrong medium, wrong field, wrong frequency band); an empirical ~1 Hz
resolution survey of six calibrated slices spanning 11 months confirms
every narrowband oscillation-like signal is engineering, shipping, or
biology. The report also covers the nearest *real* plasma↔ocean-acoustics
phenomena (sea-surface lightning, sonoluminescent cavitation) and how those
are actually localized.

## Layout

- `src/fetch_data.py` — byte-range fetch of 60 s slices from the public S3 archive
- `src/spectral_analysis.py` — calibration, Welch PSDs, narrowband line detection, spectrograms
- `src/line_catalog.py` — cross-slice line persistence, classification, physics band-mismatch chart
- `src/transients.py` — broadband impulse detection (the plasma-generated-sound signal class)
- `src/fetch_array.py` — OOI cabled-array hydrophone download via EarthScope FDSN (part 2)
- `src/localize_tdoa.py` — TDOA hyperbolic localization, validated on a M5.0 T-phase (part 2)
- `src/signature_features.py` — propulsion-signature screening: DEMON, kurtosis, spike trains (part 2)
- `src/plasma_transient_search.py` — full-bandwidth plasma-discharge transient template search (part 2)
- `src/batch_screen.py` — archive-scale daily screening campaign, disk-bounded streaming (part 2)
- `src/verify_candidates.py` — pulse-physics adjudication of template matches (part 2)
- `figures/` — generated figures used in the report
- `results/` — per-slice detected-line tables and the feature catalog

## Reproduce

```bash
pip install -r requirements.txt
python3 src/fetch_data.py && python3 src/spectral_analysis.py \
  && python3 src/line_catalog.py && python3 src/transients.py
```
