# Detecting 'The Bloop'

Can we detect the Bloop? Yes — from NOAA's own released recording.
`src/bloop_detect.py`, `figures/fig_bloop.png`.

## What the Bloop is (and isn't)

An ultra-powerful low-frequency sound recorded in **summer 1997** by NOAA's
**Equatorial Pacific autonomous hydrophone array** — loud enough to appear on
sensors >3,000 km apart, source triangulated to ~**50 S 100 W** (remote South
Pacific). Long mythologized as a "sea monster," it was attributed by NOAA
(~2012) to **icequakes**: large icebergs cracking and calving near Antarctica
(non-tectonic cryoseism). The arrival azimuth points to the Bransfield
Strait / Ross Sea / Cape Adare region. Not biological.

## Detection

Downloaded NOAA's public `bloop.wav` (and the `A53` iceberg-calving reference).
The released clip is time-compressed **16x** — that is precisely why a
low-frequency event becomes audible. Correcting the axes back to true units
(freq ÷ 16, time × 16) and running the same spectrogram + spectral-ridge
tooling used elsewhere in this repo:

- The Bloop appears as a **discrete low-frequency transient**, energy
  concentrated **~1–40 Hz (true)**, loudest ~50 s into the clip.
- Its spectral ridge **rises** through the event — the characteristic
  **upsweep** the Bloop is known for.
- The infrasonic-to-low-frequency content is the whole reason NOAA sped the
  clip up 16× to make it hearable.

The iceberg-calving reference clip (played 3×) shows the same cryogenic,
low-frequency, broadband-transient family — consistent with the icequake
attribution.

## Honest notes

- This is NOAA's *released, processed* clip, not the raw 1997 multi-station
  array data. "Detection" here means reproducing and characterizing the Bloop's
  signature from the public recording — not re-triangulating it from raw
  hydrophone streams (that array data would be the next step for a true
  from-scratch detection + localization).
- The Bloop was a 1997 event; it is not something recurring on today's live
  hydrophones. Icequakes themselves are common and *are* recorded routinely
  (e.g. on OOI/PMEL), so the *class* of sound is detectable live — "the Bloop"
  specifically is historical.
