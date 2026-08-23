# Live-MARS boom search: 2026-03-23 central-California fireball

First end-to-end run of the "pull live hydrophone data, run boom acoustic"
pipeline on a real, recently-dated event. `src/mars_boom.py`.

**Event:** green fireball over central California, ~2026-03-22 20:19 PDT =
**2026-03-23 03:19:00 UTC**, hundreds of witnesses (NASA/AMS). Many coastal
"UFO sightings" are exactly these.

**Sensor:** MARS (Monterey Bay, 36.713 N 122.186 W, ~890 m, 256 kHz, live).
The 10-min files `MARS_20260323_031000`…`034000` cover 03:10–03:50 UTC.

**Method:** stream each file, decimate 256 kHz → 250 Hz (bolide infrasound is
<~20 Hz), delete the raw file to save disk, then scan the window *after* the
flash for an impulsive low-frequency arrival (air path ~0.34 km/s ⇒ delay =
ground-range/0.34). See `figures/fig_mars_fireball_boom.png`.

## Result — clean negative (no bolide boom), transients characterized

- The 1–20 Hz detector flagged **one** transient, at **+1504 s** (03:44:03
  UTC) → implied ~511 km if it were the bolide. Too far and 25 min too late
  for a central-CA fireball; unrelated ambient infrasound.
- Inside the plausible 0–900 s window there is a **cluster of broadband
  impulsive transients (~+50 to +350 s)** — but broadband, energy to ~120 Hz.
  A distant bolide arrival is **low-frequency and dispersed**, not a sharp
  spike to 120 Hz. The broadband, clustered character = a **local source**
  (vessel cavitation / echosounder), not an airburst 100+ km away.

**Conclusion:** the fireball left no detectable hydroacoustic signature at MARS
in this window — the expected outcome for a high-altitude fireball whose track
was likely far from Monterey Bay, recorded by a deep hydrophone that couples
atmospheric infrasound weakly. A clean negative, not a dismissal.

## What would make this a real detection test

This run used a generic "central California" location. To turn it into a
targeted test:
1. Get the fireball's **actual ground track** (AMS trajectory: start/end
   lat/lon/altitude) → compute the closest-approach ground range to MARS and
   the **exact expected infrasound arrival time**, then examine that precise
   moment (seconds, not a 40-min sweep).
2. Prefer an **energetic, low-altitude bolide near the coast** (CNEOS-logged,
   >~0.05 kt) — those actually put infrasound into the water. A visual green
   fireball is usually too high/weak.
3. Corroborate any candidate against the OOI Oregon cabled hydrophones or a
   land infrasound array (same arrival, different bearing) before believing it.

The pipeline is proven; the limiting factor is finding an event energetic and
close enough to leave a trace on a deep hydrophone.
