# What the "snap rate" actually is, and whether the rise is real

## What is being measured (not a physiological rate)
The "snap rate" is a **detector count**: envelope peaks per minute above a threshold in the
2–20 kHz band, at CI01 (18 m, Channel Islands). At ~2700/min ≈ **45 snaps/second** it is the
**colony aggregate** over the hydrophone's whole detection range — thousands of snapping
shrimp snapping asynchronously — not one animal's rate. Individual shrimp snap at order
1 per few seconds up to a few/second in territorial bursts; the aggregate scales with
abundance × per-animal activity, and is modulated by **temperature (strongly), light
level (crepuscular), tidal current, and social behaviour**. There is no single
"physiological" value: 2600 and 2800/min are both just the colony aggregate fluctuating
~±10% with conditions.

## Cherry-picking check (control absolute levels)
Median snap rate per CI01 file:
```
control days: 2639, 2658, 2766, 2766, 2790, 2826, 2868, 2890, 3090 /min
incident days: 2602 (07-15), 2716 (07-26), 2714 (07-31)
EVENT day 07-16: 2772 /min   <- mid-pack, not special
```
Controls span **2639–3090/min** (windows 2004–3438); the event day is in the middle. The
anomaly was never the absolute level — it was the **within-day after/before ratio**, which
cancels absolute level. So no cherry-picking of ~2700 controls; the null uses the full range.

## Is the +7.7% rise real, or a detection-threshold artifact?
Re-measured the event-day 06:00 UTC before/after change three independent ways:

| metric | before | after | change |
|---|---|---|---|
| MAD-adaptive-threshold count | 2661/min | 2867/min | **+7.7%** |
| FIXED-absolute-threshold count | 2432/min | 2729/min | **+12.2%** |
| band RMS ENERGY (threshold-free) | −58.74 dB | −58.58 dB | **+0.16 dB (+3.7%)** |

All three move the same direction, so the rise is **not** a pure thresholding artifact — the
snap-band energy rose too. But the energy change is small (**0.16 dB**), and the count rose
more than the energy, meaning the extra snaps are individually smaller. So: a **real but
subtle increase in snap density** (a few % more, mostly small snaps), not a loud event.

## Why would shrimp change snap rate here?
Ordinary drivers — a warm-water parcel advecting past, a tidal-current shift, or the
crepuscular light cycle — routinely move colony snap rate by several % over tens of minutes.
A ~4–8% change over ~40 min is within that normal envelope. No known mechanism connects a
185 km-distant event to shrimp at CI01. The specific local driver of this particular rise
is not identified; it is a real, small, environmentally-plausible snap-density fluctuation
whose coincidence with the event time is what makes it notable, not its size.
