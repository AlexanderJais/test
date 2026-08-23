# 'Good' hydrophone spots (and the times gap)

Cross-referencing a geocoded UFO event dataset against the public hydrophone
network to pick where — and when — a passive-acoustic UAP search is worth
running. `src/ufo_spots.py`.

**Dataset:** Larry Hatch's *U database, ~18,116 geocoded records, from
[`richgel999/ufo_data`](https://github.com/richgel999/ufo_data)
(`bin/hatch_udb.json`). Every record has `LatLong` + date; 1,162 carry the
attribute **`SND: UFO sounds heard or recorded`** — the acoustically relevant
subset.

## The hard constraint: spots yes, times no

Hatch's data spans 989–**2003**. Every accessible public hydrophone starts
**2007+** (most 2014+). **Zero temporal overlap** — so this dataset tells us
*where* activity historically clustered near a sensor, but cannot supply an
event *time* that falls inside any sensor's recording window. For times, a
current database (e.g. NUFORC, which runs to the present) must be intersected
with a continuous sensor's coverage.

## Good spots (events within 50 km of a public hydrophone)

| Hydrophone | <25 km | <50 km | <100 km | SND<50 | window |
|---|---|---|---|---|---|
| SanctSound Stellwagen (Mass Bay) | 0 | **49** | 248 | 5 | 2018-2021 |
| SanctSound Monterey Bay | 8 | 22 | 84 | 0 | 2018-2021 |
| **MARS Monterey Bay** | 0 | **15** | 74 | 0 | **2015-present, 256 kHz** |
| SanctSound Florida Keys | 11 | 12 | 13 | 0 | 2018-2021 |
| NRS05 SoCal | 1 | 9 | 78 | 1 | 2014-2015 |
| SanctSound Channel Is. | 0 | 4 | 38 | 0 | 2018-2021 |
| HARP SOCAL / OOI Oregon / others | 0–1 | 0–1 | 0–101 | 0 | — |

(Counts are population-driven — coastal cities generate more sightings — so
read them as "sensor sits in an active area," not "over-water contacts.")

## Recommendation

- **Best actionable spot: Monterey Bay.** It is the only place a *top-tier
  live* sensor (**MARS** — 256 kHz, continuous, 2015-present, public) coincides
  with real historical activity (8 Hatch events within ~40 km — Santa
  Cruz/Monterey, 1950–1973, credibility up to 10). If we monitor anywhere
  forward, it's here.
- **Richest — and most acoustically on-point — cluster: Stellwagen Bank /
  Massachusetts Bay.** Most events overall (49 within 50 km) and a distinct
  set of **SOUND-reported** cases (humming / roar / whine, 1954–1990, e.g.
  Salem cr11, Lynn, Boston). Directly relevant to an acoustic search — but the
  SanctSound hydrophone there was temporary (2018-2021); no continuous public
  sensor remains.
- **Times:** not derivable here (data ends 2003). Pull NUFORC (present-day,
  geocoded) and intersect with MARS's 2015-present coverage to get candidate
  event *times* at a live sensor.

## The acoustically relevant subset

1,162 Hatch events flag reported SOUND. That subset — historical reports of a
UFO *making noise* near what is now a hydrophone site — is the most on-target
seed for a passive-acoustic study, even though it predates the sensors. It is
also a reminder that "UAP + sound" has a long anecdotal record the acoustic
method could, going forward, actually test.

## The 'times' half — resolved (and why it's nearly empty)

Attempting the sensor∩report *time* intersection ran into two walls:
1. **NUFORC lock.** Modern geolocated UAP data is essentially all NUFORC-
   derived, and NUFORC's ToS forbids scraping/redistribution (the one public
   mirror, `timothyrenner/nuforc_sightings_data`, stores its CSV in a private
   GCS bucket → 403). We do **not** scrape NUFORC. So the systematic
   intersection is unavailable through clean data.
2. **Geography.** The documented 2015+ transmedium cases (USS Omaha sphere
   2019-07-15; USS Jackson 2023-02) are **Southern California**, not Monterey.
   MARS — the best *live* sensor — has no documented over-water case to anchor.

Honest ranked candidate list (documented incident × a hydrophone that was
recording):

| # | Incident | Date | Where | Recording sensor | Dist | Status |
|---|---|---|---|---|---|---|
| 1 | USS Omaha sphere→water | 2019-07-15 | off San Diego | SanctSound CI (2018–21) | ~170 km | already analyzed |
| 2 | USS Jackson video | 2023-02 | SoCal | none (CI retired 2021) | — | no coverage |

**Conclusion:** the best {documented incident × recording public hydrophone}
overlap in clean data is the July-2019 SoCal case at the Channel Islands
hydrophones — the exact case analyzed this session. It was not an arbitrary
choice; it is the closest a documented UAP event and a public hydrophone have
come in space *and* time in the accessible record.

Clean ways forward: (a) request NUFORC data legitimately (ToS permits email
request), then run `ufo_spots.py`-style correlation against sensor coverage;
(b) a forward watch on MARS (live, 256 kHz) for new Monterey-area over-water
sightings — catch the future rather than mine the locked past.
