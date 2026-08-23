# CI01 snap-rate rise: what could the shrimp be reacting to?

Follow-up to `CI01_TEMPERATURE_CORRECTION.md`. The ~8% snap-rate rise at CI01
coincident with the 2019-07-16 (USS Omaha) event is real and
temperature-independent. This tests the mechanistic question: **are the
snapping shrimp reacting to an energy input, and if so in what channel?**

## What a hydrophone can and cannot test

A hydrophone measures acoustic **pressure**. It records an electromagnetic
field only if that field couples *electrically* into the cable/preamp (mains
hum). It cannot measure light, a static or low-frequency magnetic field,
water chemistry, or DC pressure. So this scan can only test:

- acoustic bands (infrasound 1–20 Hz through the 2–20 kHz snap band);
- mains-line pickup at 50/60 Hz, as a proxy for EM coupling on the cable.

Anything else is outside the instrument. Absence of a signal here is **not**
evidence of absence in an unrecorded channel.

## Method (`src/mechanism_scan.py`)

Event file `SanctSound_CI01_02_..._190716033643.flac` (48 kHz), over the
event arrival t0 ±90 min (359 × 30 s bins). Per bin: band energy (dB re the
window median) for each candidate band, plus the snap rate (2–20 kHz
envelope peaks > median + 8·MAD, per minute). Then:

- **lead-lag** cross-correlation of each band vs the snap rate
  (positive lag = band *leads* snaps = stimulus → response);
- **energy step** across the event (mean over t0..t0+20 min minus
  t0−20..t0 min).

A "shrimp reacted to energy X" signature would be: band X shows a coincident
energy **increase** that **leads** the snap-rate rise.

## Result (reproducible run)

Energy step at the event:

| band | Δ energy | 
|---|---|
| infrasound 1–20 Hz | −0.02 dB |
| low 20–100 Hz | **−1.36 dB** |
| ship 100–500 Hz | **−1.66 dB** |
| mid 0.5–2 kHz | −0.02 dB |
| SNAP 2–20 kHz | +0.19 dB (the shrimp) |
| EM 60 Hz | +0.73 dB |
| EM 50 Hz | +0.33 dB |

Lead-lag: the acoustic bands that "lead" are the ones **decreasing** (ship
r = −0.47, low r = −0.28) — anti-correlation, not a driving stimulus. The
two mains lines (50/60 Hz) are the only bands that rose, but both have the
**snaps leading them** (lag ≈ −10 min), so neither precedes the response,
and the rise is only ~0.3–0.7 dB.

## Interpretation

- **Not recorded acoustic energy.** A sound the shrimp heard would appear
  here and lead the response. Instead the acoustic field went *quieter*
  (low/ship bands down 1.4–1.7 dB) while the shrimp snapped *more*.
- **Not electrical/EM pickup on the cable.** No coincident-and-leading
  mains-line change; the 60 Hz wiggle lags the snaps and is sub-dB.

Two honest possibilities remain, not collapsed into one:

1. the trigger was in a channel this sensor cannot record (light, magnetic
   transient, particle motion / low-frequency pressure below the high-pass,
   chemistry); or
2. the rise was not a direct reaction to any external stimulus (internal /
   behavioural, or coincident-by-timing but not caused by a detectable
   input).

What is ruled out cleanly: the shrimp did **not** snap more because the water
got acoustically louder.

## Replication at CI04 and CI05 (`src/mechanism_replicate.py`)

The same scan was run at the two other Channel Islands stations (independent
cables/sites, ~15–30 km away), with a concatenated reader so CI04's event
window could be scanned across its 6 h file boundary. Result: **the CI01
signature does not replicate.**

| station | snap-rate step at event | low 20–100 Hz | ship 100–500 Hz |
|---|---|---|---|
| CI01 | **+8%** (rise) | −1.4 dB | −1.7 dB |
| CI04 | **−37%** (fall) | +1.0 dB | +0.6 dB |
| CI05 | **−24%** (fall) | −1.0 dB | −0.3 dB |

- Neither CI04 nor CI05 shows the snap-rate **rise**. Both snap rates
  **decline** across the event window — consistent with a normal night-time
  diel decrease (the event is ~23:00 local). CI01 rose *against* that
  expectation, which is what made it notable; CI04/CI05 look ordinary.
- The specific "quiet field / louder shrimp" divergence is absent: at CI04
  the bands rose *with* a broadband transient while snaps fell; at CI05 the
  bands and snaps both fell. Each station did its own thing — no coherent
  cross-station acoustic signature.
- Both CI04 and CI05 do show a strong broadband / 60 Hz transient near (but
  not exactly at) the event time; at different clock offsets, most likely
  local passing vessels rather than a common source.

**What this means for the CI01 claim.** The event-coincident, temperature-
independent snap-rate rise at CI01 stays a **single-station** observation.
That neither proves it is noise (it is statistically robust and change-point-
aligned at CI01) nor supports a *regional* physical stimulus — a genuinely
propagating environmental/acoustic driver would more likely leave a coherent
trace across the three stations, and it does not. The honest bound: real and
local to CI01, not corroborated elsewhere.
