# Event-Directed Test on a Ground-Truth Transmedium Object — 2020 NE Pacific Bolide

## Why a bolide
Documented UAP water-entry reports near a capable public sensor are scarce, and none has
a precise, catalog-grade time+position. A **bolide** does: NASA's CNEOS fireball database
gives date, time, latitude, longitude, radiated energy, velocity and airburst altitude.
A bolide is a genuine high-speed transmedium object (space -> atmosphere -> [ocean]) whose
energy couples acoustically into the sea — an ideal ground-truth test of whether a cabled
hydrophone array can detect and localize such an event, and a real validation of the
entry/transit detection pipeline.

## Candidate selection (CNEOS, NE Pacific box 30-56 N, 118-140 W)
Seven located fireballs; ranked by radiated energy and nearest cabled sensor:

| date (UTC) | lat | lon | E (kt) | v (km/s) | alt (km) | nearest | note |
|---|---|---|---|---|---|---|---|
| 2020-05-09 02:56:11 | 44.8 | -131.0 | 3.53 | 14.5 | 31.2 | **OOI** | **~150 km to Axial cluster; OOI operational; GPS-timed** |
| 2022-06-07 22:53:17 | 40.8 | -127.1 | 2.82 | 24.9 | 41.1 | OOI | ~450 km |
| 2000-08-13 03:00:32 | 36.7 | -127.8 | 6.50 | — | — | MARS | pre-archive (2000) |
| 2004-06-03 09:40:12 | 48.9 | -120.4 | 5.40 | — | — | ONC | pre-ONC / inland |
| 2011-01-08 19:38:44 | 33.3 | -125.8 | 3.28 | — | 48.0 | MARS | pre-archive |

**Chosen: the 2020-05-09 event.** Offshore, in the OOI operational era, and the OOI Axial
hydrophone cluster (AXBA1/AXCC1/AXEC2, 45.9 N 130.0 W) is only ~150 km away while the
Hydrate Ridge pair (HYS14/HYSB1, 44.5 N 125.3 W) is ~460 km — an excellent TDOA baseline —
and OOI is GPS-disciplined, so the arrival-time test is rigorous (unlike the free-running
SanctSound SoundTraps in the Omaha case).

## Method
Fetch 30 min of all five OO/HDH channels (200 Hz, response-removed to Pa) from 3 min before
the burst. The airburst at 31 km altitude couples to the ocean via its infrasound airwave;
at celerity ~0.30 km/s over 150-460 km the ocean arrival is expected ~8-27 min after the
burst. Scan the 2-20 Hz ocean-infrasound band for coherent transients, then cross-correlate
station pairs and grid-search the source, comparing to the CNEOS position.

## Honest expectations
Bolide-to-hydrophone coupling is usually studied with dedicated infrasound arrays; a 3.5 kt
airburst at ~150 km is plausibly detectable in the ocean-infrasound band but not guaranteed,
and the celerity/ducting makes timing uncertain by minutes. A clean coincident transient
localizing back to 44.8 N 131.0 W would be a strong positive; a null is an honest sensitivity
statement, not a failure of the method. Result and figure recorded alongside this file.
