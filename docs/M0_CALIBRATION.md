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

---

# M0 Calibration — Quench-Tail Detector (Discriminant D4)

**What it looks for:** after an entry doublet, a sustained IRREGULAR DECAYING
impulsive tail (film-boiling "chugging", 5–50 kHz pulses) — the hot-body flag. A cold
body is silent after pinch-off. Gates: elevation over local baseline, kurtosis > 4,
ISI CV > 0.5 (rejects the metronomic echosounder), ≥ 4 sustained windows, decaying trend.

**Calibration caught one bug:** the first version measured level in the full 100 Hz–50 kHz
band, dominated by low-frequency ocean noise, so even 8× injections showed ~0 dB elevation
(and injection was mis-scaled by peak not RMS). Fixed by detecting in a 5–50 kHz "chug band"
where the impulses beat ambient, and scaling injections to unit RMS.

**Measured result:**

| Injected level (× chug-band baseline) | t_q = 15 s | t_q = 35 s |
|---|---|---|
| 1× | 17% | 83% |
| 2× | 92% | 100% |
| ≥ 4× | 100% | 100% |

- **False triggers: 0/60 random triggers on unmodified ocean, and 0/49 on each of the
  heaviest biosonar/echosounder slices.** The onset-then-decay-over-local-baseline
  requirement rejects continuous biosonar (no step-up at a random trigger). Operationally
  the tail detector fires only *after* a doublet, so real FP is lower still. This is the
  cleanest of the three limbs and the sharpest hot-vs-cold discriminator.

---

# M0 Calibration — Doppler-Drift Tracker (Discriminants D1/D3)

**What it looks for:** a narrowband tone/PRF line sweeping the CPA S-curve
f_obs(t) = f0·c/(c+v_r); fits (f0, v, R_cpa, t_c) and reads out **speed** and range.
Gates: good S-curve fit (residual < 0.2% of f0), significant swing (> 0.3%), physical
v < 500 kn.

**Measured result** (speed recovered within 25% of truth):

| Tone SNR (28–44 kHz band) | 50 kn | 150 kn | 400 kn |
|---|---|---|---|
| 0 dB | 67% | 67% | 0% |
| +6 dB | 67% | 67% | 67% |
| +12 dB | 67% | 67% | 67% |
| +20 dB | 100% | 100% | 67% |

- **False transits: 0/12 slices** — ocean tones (ship lines, whistles) do not fit the
  monotonic plateau→inflection→plateau S-curve with a physical speed. Cleanest possible
  FP result; this is the standalone-strongest limb.
- **Two honest limits:** (i) the 400 kn / 0 dB cell is 0% because a 400 kn transit past
  R_cpa = 700 m has a ~3 s CPA core — a fast transit is intrinsically brief (the coverage
  lesson, now on the kinematic side). (ii) The ~67% plateau traces to the 38 kHz
  echosounder falling inside the 28–44 kHz tracking band on some days; a production run
  chooses tracking bands around known lines or tracks multiple ridges.

---

# M0 summary — gate status

| Detector | Recovery | False alarms | Standalone? | Role |
|---|---|---|---|---|
| Entry-doublet (D2) | 74% quiet / storm-blind @3× | ~1.6 FP/min | **No** | coincidence leg (needs LF + catalog) |
| Quench-tail (D4) | 92–100% @≥2× | 0/60, 0/49 in storms | Yes (gated on doublet) | hot-vs-cold flag |
| Doppler tracker (D1/D3) | 67–100% | 0/12 | **Yes** | kinematic core (speed readout) |

**All three M0 detectors PASS the gate (each has a measured sensitivity curve + FP rate).**
The calibrations turned two design assumptions into measured facts: the entry-doublet HF
limb is storm-blind and high-FP, so it is only a coincidence leg (M2/M3 mandatory); the
Doppler tracker and quench-tail are clean enough to carry weight on their own. **M0 complete
→ proceed to M1 (single-site MARS full-coverage month pilot), treating single-sensor
doublets as coincidence candidates, not detections.**
