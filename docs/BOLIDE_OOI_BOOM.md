# Meteor into the ocean: 2020-05-09 bolide on OOI Axial — honest negative

The real transmedium target (earthquakes are not meteors): an object from the
sky depositing energy into the sea. `src/bolide_ooi.py`,
`figures/fig_bolide_ooi.png`.

## Target (ground truth)

CNEOS: bolide airburst over open ocean, **2020-05-09 02:56:11 UTC, 44.8 N
131.0 W, radiated 3.5e10 J (~0.12 kt)**. The three OOI **Axial Seamount**
cabled hydrophones (AXBA1/AXCC1/AXEC2, HDH, 200 Hz) sit **~150 km** away and
were recording — the closest an accessible public hydrophone comes to a
2015+ ocean bolide (measured against the full CNEOS catalog). Cross-station
coherence on three sensors is the detection test.

## Result: below detection

Searched airburst+60…+540 s (bracketing both the in-water-coupled ~150–200 s
and atmospheric ~440 s arrival bands). All three stations show only
fluctuating background; the strongest **coherent** (all-station) transient is
**min-z = 2.0** (threshold 6). No meteor boom.

## Why — honest reasons, not excuses

1. **Small energy.** 0.12 kt is a modest bolide. Hydroacoustic bolide
   detection in the literature favors larger events or very close range.
2. **Weak air→water coupling.** An airburst is atmospheric (~30 km up); only a
   small fraction of its energy couples through the sea surface into the ocean.
3. **Noisy site.** Axial Seamount is an **active submarine volcano** — its
   hydrophones sit in a hydrothermally and seismically noisy caldera, raising
   the detection floor. It is the closest accessible sensor but not a quiet one.

The energy/range trade-off confirms this was the best available shot: the more
energetic 2015+ ocean bolides (0.5–0.7 kt) are 1700+ km away, which is *worse*
(amplitude falls with range faster than it rises with energy).

## What it would actually take

- A **large bolide** (≳1 kt — e.g. the 2018-12-18 Bering Sea event, ~173 kt),
  and/or
- a **low-noise, SOFAR-optimized IMS hydroacoustic station** (CTBTO HA-series:
  Wake, Ascension, Diego Garcia, Juan Fernández…). Those routinely detect
  ocean bolides — but their data is access-controlled (IDC/vDEC), not freely
  pullable like MARS/OOI.

## The honest bottom line

Meteors *do* put sound into the water — the physics is sound — but with freely
pullable public hydrophones and the small ocean bolides available during their
coverage, a meteor-into-ocean signal sits **below detection**. That is the real
gap: not that the method fails, but that the accessible-sensor × available-event
combination doesn't cross threshold. The M6.0 T-phase crossed it easily
(`MARS_TPHASE_BOOM.md`) because an M6.0 delivers ~10⁶× the acoustic energy of a
0.12 kt airburst and couples directly in-water. Same detector, honest contrast.
