# Event-Directed Investigation — 2019 USS Omaha "Splash" Incident

## The incident
Pentagon-confirmed Navy video (USS Omaha, CIC footage) shows a spherical ~2 m object
descending to the sea surface with a "splash" call. Reported time/place: **2019-07-15
~23:00 PDT = 2019-07-16 ~06:00 UTC**, warning area off Southern California, **32.4894 N,
119.3647 W**.

## Why this event, and the honest limits going in
It is the only well-known transmedium report that falls inside the public-hydrophone era
in a region with archived coverage. But the nearest public sensors — the **SanctSound
Channel Islands** array (CI01/CI04/CI05) — are **167–185 km north** of the reported
position. At that range:
- only a strong low-frequency component survives (we search 10–500 Hz);
- a small (~2 m) object's water-entry may be **below detectability regardless of whether
  the event was real** — so a null bounds source level, it does not disprove the event;
- islands/bathymetry can block paths; SoundTrap recorder clocks drift minutes (we use a
  generous ±20 s TDOA tolerance).

## Method
Read 3 h of archived audio per station bracketing the predicted acoustic arrival
(event time + range/1.487 km s⁻¹ ≈ 112–124 s). Scan a 10–500 Hz "entry-thump" band for
impulsive transients; rank by SNR. Then the decisive test: **TDOA geometry**. A source at
the reported position must produce specific inter-station arrival delays —
**CI01−CI04 = +12.5 s, CI04−CI05 = −2.3 s** — independent of the (uncertain) absolute
event time. Test every triple of thumps against that geometry.

## Result — no localizable event at the reported position
- Each station shows ~12 LF thumps over the 3 h window (ships, distant impulses, biologics).
  CI04 has one at 06:01:23 UTC, only −29 s from the naive predicted arrival (SNR 16).
- **But no triple of thumps satisfies the required TDOA geometry within ±20 s.** The lone
  near-arrival CI04 thump has no CI01/CI05 partners at the offsets a co-located source
  demands.
- Chance expectation of a spurious TDOA-consistent triple in this window is **0.02**, so a
  genuine coincident source at the event position would have stood out clearly. It is absent.

![Omaha event window](figures/fig_omaha_tdoa.png)

## Interpretation (stated conservatively)
At the three public hydrophones 167–185 km away, **there is no acoustically-localizable
water-entry event consistent with the reported time and position** of the USS Omaha
incident. This is the *expected* outcome for a ~2 m object at ~170 km even under ordinary
water-entry physics — the array is simply too distant and too far north — so the finding is
a **source-level upper bound at the reported position, not evidence against the event**.

Calibration caveat: SanctSound absolute calibration (ST500 sensitivity) was **not** applied
here; detection and TDOA use relative SNR, which is valid for coincidence but means the
upper bound is order-of-magnitude, not a calibrated pascal figure.

## What would actually test such an event
A hydrophone (or, better, a GPS-timed array) **within ~30 km** of the warning area — e.g.
a future asset in the SOCAL range complex — at ≥ tens of kHz bandwidth. The value of this
exercise is the **template + method**: given a time and place near a capable sensor, the
entry-doublet + quench + TDOA pipeline turns an anecdote into a falsifiable acoustic test.
