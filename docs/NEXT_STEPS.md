# Next Steps — Executing the Transmedium Plasma-Sheath Search

Derived from [SEARCH_PLAN.md](SEARCH_PLAN.md). Ordered by value ÷ cost. Each milestone
is independently shippable and ends in a committed, testable artifact. "Have" = already
built this session; "Build" = new.

## Guiding rule
Never hunt with an unmeasured detector. Every new detector is injection-tested against
synthetic truth (recovery matrix + false-positive rate) **before** it touches archive
data — same discipline as the PRI miner (floor = 1.5× local threshold, 67–92% recovery).

---

## M0 — Build & calibrate the three missing detectors  *(highest priority, ~1 day compute)*
The signature spec's top discriminants need detectors we don't have yet.

- **D2 Entry-doublet detector** — slam shock + cavity pinch-off collapse pulse separated
  by `t_p` (0.1–0.55 s); matched-filter/cepstral pair-finder, outputs `t_p` → implied
  body size, plus anti-coincidence with airgun (fixed 10–15 s) and echosounder PRIs.
  *Build `src/entry_doublet.py`.*
- **D4 Quench-tail detector** — after any doublet, test a long window (up to ~1000 s) for
  an impulsive high-kurtosis tail with irregular ISI (hot-body flag; cold entries are
  silent after pinch-off). *Extend the transient chain with ISI stats.*
- **D1/D3 Doppler-drift tracker** — reassigned-spectrogram tracker for tones/PRF combs
  chirping 4–300 Hz/s through CPA; require identical fractional Doppler across all lines.
  *Build `src/doppler_track.py`.*
- **Calibrate all three** with `injection_test.py`-style synthetic injection into real
  MARS slices → per-detector recovery matrix + FP rate. **Gate: no detector proceeds
  past M0 without a measured sensitivity curve.**

## M1 — Single-site full-coverage pilot on MARS  *(proves the pipeline; ~40 core-h)*
- Run the tiered LF/HF ladder **continuously over one bounded month** (all 144×30 files),
  reusing `deep_dive.py`'s robust streaming, now emitting doublet / quench / Doppler /
  PRI candidates, not just the low-freq band.
- Deliver a candidate stream + measured false-alarm rate over a full month → sets the
  detection threshold and the per-day event budget for a full-year run.
- **Decision point:** month pilot false-alarm rate tolerable → greenlight the 1-year run
  (24.2 TB egress, ~1,300 core-h, ~39 h on 32 cores).

## M2 — Cross-band coincidence + array speed (the two smoking-gun discriminants)
- **D2 cross-band gate:** co-register MARS HF (0.1–30 kHz doublet limb) with OOI LF
  (3–100 kHz thump) within ≤0.1 s travel-time-corrected → rejects single-band mimics.
- **D1 quiet-fast track:** for any surviving candidate, pull the OOI 5-station window and
  run the **validated TDOA pipeline** (`localize_tdoa.py`) to get a speed track; flag the
  impossible region (fast + low broadband level). This is the strongest discriminant and
  reuses code already proven to 25 km at 240 km.

## M3 — External anti-coincidence  *(cheap, high rejection power)*
Cross-match candidate timestamps against public catalogs to kill mimics / support a hit:
- CNEOS fireball/bolide database, WWLLN/GLD360 lightning, space-debris reentry TLEs
  (Space-Track), USGS seismic (T-phase), and — where public — Navy range schedules.
- *Build `src/coincidence.py`* (timestamp → catalog lookups, ±window).

## M4 — Proximity sites for the report-region search
- **SanctSound Channel Islands CI01–CI05** fetcher (NCEI/AWS) — search the 2004 Nimitz /
  2019 Omaha-Russell region directly. *Build `src/fetch_sanctsound.py`.*
- Add ONC 512 kHz (above-Nyquist confirmation tier) and Orcasound (near-surface boom/
  breach limb) as confirmation-only feeds.

## M5 — The science output, detection or not
Regardless of whether a candidate survives: convert full-coverage coverage + measured
sensitivity into a **population upper bound** on the hypothesis — an exclusion curve in
(radiated power × duty cycle × events-per-year) space, per site. A rigorous null is the
primary deliverable; a surviving candidate escalates to M2/M3 confirmation.

---

## Sequencing & cost summary
| Milestone | New code | Compute | Gate to next |
|---|---|---|---|
| M0 detectors + calibration | 2 scripts + ext | ~1 day | measured sensitivity per detector |
| M1 MARS month pilot | reuse deep_dive | ~40 core-h | tolerable false-alarm rate |
| M2 cross-band + TDOA | 1 gate + reuse | per-candidate | — |
| M3 coincidence | 1 script | minutes | — |
| M4 proximity sites | 1–2 fetchers | ~40 core-h/mo/site | — |
| M5 population bound | analysis | minutes | (final) |

**Recommended immediate action:** start M0 — build and injection-calibrate the
entry-doublet detector first (it is the singleton-event discriminant a transit most
depends on, and the one most exposed to the coverage problem).
