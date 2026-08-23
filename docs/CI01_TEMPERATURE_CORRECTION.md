# Correction: temperature does NOT explain the CI01 snap-rate rise

## What I got wrong
In the previous step I regressed snap rate on a **pooled, months-long** temperature
correlation (r=-0.555) and used it to claim the event's 0.6 °C cooling "explained" +5.6% of
the +7.8% snap rise, leaving a non-significant residual (p=0.099). That was a
**wrong-timescale over-control**: the pooled correlation is dominated by seasonal/tidal
covariation, not the instantaneous shrimp response over 40 minutes. Using it to erase a
fast event signal was not valid. (The raw control-null p=0.006 already accounts for normal
temperature variability, since control days also have temperature swings — so the residual
step double-counted on a bad basis.)

## The correct test — matched-temperature-change control
No regression assumption: take control-day windows that had the **same ~0.6 °C cooling** as
the event, and ask how often their snap rate rises as much as the event's.

- Event: ΔTemp = -0.575 °C, ΔSnap = **+8.2%**.
- Short-timescale (±20 min) ΔTemp↔ΔSnap correlation is only **-0.18** (not -0.55) — the
  pooled figure massively overstated the instantaneous temperature effect.
- Control windows with ΔTemp in -0.575 ± 0.3 °C: **n=87, median ΔSnap +2.3%, and 0 of 87
  (0.0%) reach the event's +8.2%.** Matched empirical p ≈ 0.000.
- Unconditional: 0.4% of all control windows reach +8.2% (p=0.004).

![matched-temperature test](figures/fig_matched_temperature.png)

## Corrected conclusion
**Temperature does not explain the CI01 snap-rate rise at the Omaha event.** Conditioning on
matched cooling makes the event *more* anomalous, not less: no control window with the same
temperature drop produced a comparable snap rise. My earlier "temperature accounts for most
of it" statement is **retracted**.

## Honest current status of the lead (unchanged facts, no fabricated confounds)
The CI01 rise is a genuine, temperature-independent, statistically robust anomaly
(p≈0.004 unconditional; 0/87 matched-cooling controls reach it; change-point within 24 s of
the predicted arrival; not a diel cycle). The remaining real limitations are factual, not
invented: it is single-station (CI04/CI05 did not replicate for this event), and it appeared
at 1 of the 4 July-2019 incident nights. Those are legitimate caveats. The temperature
"explanation" was not — and it is withdrawn.
