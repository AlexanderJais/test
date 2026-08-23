# 'The Bloop': characterization, re-triangulation feasibility, and a critical assessment

`src/bloop_detect.py`, `figures/fig_bloop.png`.

## What the Bloop is

An ultra-powerful low-frequency sound recorded in **summer 1997** on NOAA's
**Equatorial Pacific autonomous hydrophone array**, source triangulated by
NOAA to ~**50 S 100 W**. Attributed (~2012) to **icequakes** — large icebergs
cracking/calving near Antarctica — not a marine animal. NOAA's public clip is
time-compressed **16x** (a linear speed-up), which is why a low-frequency event
becomes audible.

## Can we RE-TRIANGULATE it? No — and this is a hard result, not a shrug

1. **The public clip is a single mono channel.** Triangulation *in principle*
   needs ≥3 stations' relative arrival times (TDOA) or a beamforming array;
   one channel carries **zero** directional information. The Bloop cannot be
   located from `bloop.wav`, period.
2. **The raw 1997 array data is not reachable here.** The NCEI passive-acoustic
   archive holds only modern projects plus Arctic **Fram Strait** (`fram/` — the
   wrong array and wrong hemisphere); there is no 1997 Equatorial-Pacific
   holding, and IRIS/EarthScope did not serve it. NOAA's "50 S 100 W" is *their*
   array fix; we cannot independently reproduce it without their raw waveforms.

So any honest "re-triangulation" is blocked at the data layer. What we *can* do
rigorously is assess the signal we have.

## Critical assessment of the signal (from the released clip)

Measured, true units (freq ÷16, time ×16):

- **Brief, low-frequency, broadband transient.** Core event ~20 s in this clip
  (fuller NOAA versions ~1 min), energy concentrated **~1.7–34 Hz**, crest
  ~6–7× the background. Consistent with a cryogenic/impulsive source; far too
  low and too loud (see caveat) for any biological source.
- **No robust "upsweep."** The frequency-vs-time trend in the core is **weak and
  method-dependent** (+0.03 to −0.22 Hz/s across methods) — I cannot confirm a
  clean sweep from this clip. *Correction to an earlier claim of mine:* the
  crisp upsweep is the defining feature of NOAA's **separate** sound literally
  named **"Upsweep"** (a different phenomenon, ~54 S 140 W). Popular Bloop
  write-ups conflate the two; I did too, and am correcting it.
- **What the clip CANNOT support:** absolute source level (the file is
  amplitude-normalized), so the famous "loudest sound ever / heard 5,000 km
  apart" claim is *not verifiable from this clip* — it rests on NOAA's raw,
  calibrated multi-station data. And, again, location.
- **Processing note:** the 16× speed-up is linear time-scaling — it preserves
  sweep sign/shape (does not manufacture a sweep) but pushes the true <35 Hz
  content up to <560 Hz in the file, so everything above that in the file is
  empty/noise, not signal.

## Critical assessment of the icequake attribution

**Well supported by:** signature match to the abundant icequake recordings NOAA
later made near Antarctic icebergs (e.g. A53a); a source azimuth pointing south
toward Antarctica; austral-summer (calving-season) timing; and the absence of
any biological source capable of a ~1-min, basin-scale, 1–35 Hz transient.

**Caveats worth stating:** the ~50 S fix is approximate and lies *north* of the
main iceberg belt (the interpretation posits propagation from further south);
and the attribution is by **analogy / pattern-match** to later events, not a
contemporaneous ground-truth of a specific 1997 berg. It is a strong,
well-evidenced inference — not a closed-case single-event proof, and popular
retellings overstate both its mystery and its "loudest ever" superlatives.

## Bottom line

We can characterize the Bloop from NOAA's clip (brief, low-frequency, broadband,
icequake-consistent) and we can state firmly *why* it can't be re-triangulated
from public data. We cannot verify its location or absolute level, and the
"upsweep" descriptor does not survive scrutiny of this clip. Honest limits, not
hand-waving.
