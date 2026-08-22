# M0 Calibration — Entry-Doublet Detector (Discriminant D2)

**Gate rule:** no detector proceeds past M0 without a measured sensitivity curve and
false-alarm rate. The entry-doublet detector (`src/entry_doublet.py`) has now been
injection-calibrated (`src/injection_doublet.py`) against physically modeled water-entry
doublets injected into 12 real MARS slices spanning quiet, biosonar-storm, and
echosounder conditions.

## What the detector looks for
An isolated pair of transients: **slam** (broadband, low-pass, spectral corner
f_c ≈ 240/a Hz) → quiet interlude → **pinch-off pulse** at t_p = 2·√(a/g), optionally
with trailing Minnaert ring. The slam→pinch delay gives body size **a ≈ g·(t_p/2)²**.

## Calibration process caught two real design failures *before any hunt*
1. **False-positive blowout (109 FP in 12 clean runs).** The first gates (`lf_frac`,
   centroid cap) did not separate slams from odontocete clicks — measurement showed real
   MARS clicks carry *high* low-band energy (median lf_frac 0.33). Fixed by switching to
   the **low/high band-energy ratio** (E[0.3–8 kHz]/E[15–60 kHz]): a slam is low-pass
   (ratio ≫ 1), a click is HF (ratio < 2). Empirically **0% of real dolphin clicks
   exceed 2.0** → the correct separator. FP fell 109 → 19.
2. **Recovery decreasing with amplitude (a backwards bug).** A loud pinch-off's own
   Minnaert ring registered as many secondary impulses, tripping the isolation gate;
   and the clean synthetic slam's rising centroid tripped a wrong centroid cap. Fixed by
   counting only *comparable-amplitude* impulses for isolation and by dropping the
   centroid gate. (The synthetic slam was also corrected to be physically low-pass.)

## Measured result
| Condition | Recovery @1.5× thr | @3× | @5× |
|---|---|---|---|
| **Quiet slices** | 30% | **74%** | 65% |
| **Storm slices** (dense biosonar) | 0% | 9% | 20% |

- **False alarms:** 19 doublets in 12 clean 60 s runs = **~1.6 FP/min single-sensor**;
  the grade-A (ring-confirmed) subset does **not** reduce this (12 of 19 were grade-A) —
  low-frequency ocean impulse pairs mimic the full template on one sensor.
- Recovery is body-size dependent (best at a ≈ 0.15 m, t_p ≈ 0.25 s); the smallest body
  (a = 0.05 m, f_c ≈ 4.8 kHz) is hardest because its slam corner approaches the click
  band — physically honest (a 5 cm entry *is* click-like).

## Two conclusions that steer the search
1. **The HF entry-doublet limb is blinded in biosonar storms (0–20%) and works in quiet
   conditions (74%)** — exactly as the signature spec predicted, and the reason the D2
   discriminant is defined as a *cross-band* gate: the OOI LF thump limb survives storms
   that blind the MARS HF limb.
2. **At ~1.6 FP/min, a single-sensor entry-doublet detector is not usable standalone.**
   Its value is only as **one leg of the D2 coincidence** — HF doublet *and* LF thump
   (travel-time corrected) *and* external-catalog anti-coincidence. This is not a
   detector weakness to fix by tuning; it is why M2 (cross-band + array) and M3
   (coincidence catalogs) are *required*, not optional. The calibration converts a design
   assumption into a measured requirement.

## Status
Entry-doublet detector: **M0 gate PASSED** (measured curve exists). Next M0 detectors:
quench-tail and Doppler-drift tracker, then M1 pilot. The FP rate mandates that the M1
pilot treat single-sensor doublets as *candidates for coincidence*, never as detections.
