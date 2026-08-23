# Open-ended anomaly sweep — Channel Islands hydrophones, 2019-07-16 window

"While we have the data, any other anomalies?" A hypothesis-free pass over
the four event-window files (CI01, CI04 ×2, CI05; 48 kHz, ~24 h total):
Long-Term Spectral Average (`fig_ltsa.png`) + per-30 s features
(impulsiveness = crest/kurtosis, band levels, tonality, quantization
richness) via `src/explore_ltsa.py`, then zoom/sonify the standouts via
`src/transient_zoom.py`. Findings, ranked by how interesting they are —
and honestly labelled.

## 1. A dolphin pod at CI05, 05:30 UTC — the clearest biology
`figures/transients/CI05_kurt963_0530.png`, `sound/CI05_kurt963_0530.wav`.
The single most impulsive 30 s in the whole dataset (kurtosis **963**,
crest **180**). It is a marine-mammal encounter: sharp broadband vertical
spikes = **echolocation clicks**, plus wavy 10–20 kHz frequency-modulated
tones = **whistles**. ~30 min before the UAP-event arrival. Unambiguous,
mundane in the good sense — real animals.

## 2. A ship passage at CI05, ~03:15–04:45 UTC
`figures/transients/CI05_ship_0400.png`. The broad hump centred <1 kHz in
the LTSA that rises and falls around closest approach; the waveform is
continuous and Gaussian (crest 4, kurtosis 0). Textbook vessel — the
clearest ship in the data.

## 3. A cross-station impulsive burst cluster, ~05:28–06:53 UTC (CI04↔CI05)
`figures/transients/CI04_kurt193_0605.png`, `sound/CI04_kurt193_0605.wav`.
A run of low-frequency, spindle-shaped bursts (e.g. the CI04 06:05 doublet,
kurtosis 193) that flag on **both** CI04 and CI05 within ~10–80 s of each
other — consistent with a shared low-frequency impulsive source in the
CI04–CI05 area (a fish chorus, a whale, or a moving source). It straddles
the event window, which is why it stands out — but its origin is genuinely
**ambiguous**: it is *not* dolphin (no high-freq clicks/whistles) and *not*
an obvious continuous vessel. Flagged, not over-claimed.

## 4. Persistent narrowband tonal lines
High-resolution averaged spectra (10–200 Hz):
- **CI01: 86.4 Hz (+21 dB above background) with a harmonic at ~173 Hz** —
  the strongest tonal anywhere; a fundamental+2f pair reads as machinery /
  platform self-noise or a continuous distant source.
- A **32.2 Hz line appears on both CI01 and CI04** at nearly the same
  frequency — a shared line across independent stations.
- CI05 shows a 35–56 Hz cluster (42.5 / 49.8 / 35.2 / 55.7 Hz).
These are **constant all night** (horizontal in the LTSA), which argues for
mooring/instrument self-noise or continuous shipping rather than biology
(which would be intermittent). Not resolved further, but catalogued.

## 5. Data-quality characteristics (not anomalies, but worth stating)
- Every channel carries a **DC offset** (CI01 ~6×10⁻³, CI04/CI05 ~3×10⁻³) —
  larger than the AC std at times; must be removed before crest/kurtosis
  (Welch's per-segment detrend hides it from the PSD).
- **CI04 and CI05 are low-gain / quiet** (AC std 3–5×10⁻⁴ vs CI01's
  1.3×10⁻³) and at times ride near the 16-bit **quantization floor**
  (~50–150 distinct levels per 30 s). This is why their snap counts are far
  below CI01's — a site/gain difference, not a fault.

## Note on the event
The CI01 06:01 window has high kurtosis (517) too, but that is just the
snapping-shrimp chorus (energy in 3–20 kHz) — normal for CI01, not a
discrete transient. Nothing in this sweep points to a discrete acoustic
event at the UAP time beyond what the earlier CI01 snap-rate analysis
already documented.
