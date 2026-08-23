# Aguadilla 2013 thermal-video UAP — wind-vs-motion analysis

A standalone investigation (separate from the hydrophone/acoustic UAP work).
This is **not** an acoustic case: the 2013 Aguadilla, Puerto Rico footage is a
thermal (mid-wave IR) video from a U.S. Customs and Border Protection aircraft,
with no hydrophone anywhere near. It is analyzed here on its own because the
lantern-vs-anomaly question has one decisive, checkable test.

## The decisive test

A wind-borne object (sky lantern / balloon) must drift **at wind speed** and
**downwind**. A self-propelled craft need not. So the test is: does the
object's documented motion match the actual wind?

**Wind, pulled independently** from the Rafael Hernández Airport (TJBQ) METAR
archive (Iowa Environmental Mesonet) for the incident hour
(2013-04-26 ~01:22 UTC = 21:22 AST, 25 Apr): ~8–14 mph from the E/ENE through
the evening, **≈9.8 mph from the ENE at the incident** — matching the value
AARO cites.

**Object motion** (AARO STK reconstruction): **8 mph toward the SW**, over
land, ~200 m altitude.

→ The object's motion matches the measured wind vector in **both speed and
direction** = wind-driven drift, not propulsion.

The three "anomalous" impressions each reduce to a known effect:
- **~120 mph** (SCU claim) → motion parallax from the orbiting, zooming
  aircraft (plane–object range tripled during the clip);
- **splitting in two** → two objects the whole time, seen from a changing
  (side-on → top-down) view angle;
- **entering the water** → thermal crossover near sunset (19:48 local; video
  21:22 local), fading the cooling lanterns into the ocean's temperature.
  AARO: the objects did not enter the water.

## Honest residual uncertainty

The lantern *identity* is moderate-confidence (poor video; AARO confirmed local
resorts routinely release sky lanterns). The parallax/range model is AARO's
reconstruction; the raw frames + aircraft track would be needed to redo the
photogrammetry independently. But the **wind-vector match holds regardless of
whose range model you trust** — it only requires the object to head downwind,
which the geometry gives.

## Reproduce

```
python3 src/aguadilla_wind.py    # regenerates figures/fig_aguadilla_wind.png
```

The script is self-contained (the METAR observations are embedded). Wind data
source: Iowa Environmental Mesonet ASOS archive, station TJBQ.
