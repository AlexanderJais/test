# Unbiased Case-vs-Control Analysis — USS Omaha Event, Closest Public Mic

Rather than hunting a predefined signature, this asks the data an open question: at the
closest public hydrophone to the 2019 USS Omaha event, does the event window differ from
ordinary recordings of the *same* mic, in **any** band or feature? Case = the 60 s window
at each station's predicted acoustic arrival; control = every window >30 min away (that
mic's own natural variability); features = 16 log PSD bands (10 Hz–24 kHz) + broadband
level, spectral centroid, envelope kurtosis, click count; scored by robust z vs controls
and by multivariate anomaly **rank** among all windows.

## Results

| station | range | event-window anomaly rank | features \|z\|>3 | verdict |
|---|---|---|---|---|
| **CI04** (closest) | 167 km | **#20 / 144** | 4 (LF bands + kurtosis) | mild, explained, uncorroborated |
| CI01 | 185 km | #65 / 72 | 0 | unremarkable |
| CI05 | 170 km | #59 / 72 | 0 | unremarkable |

**Where the event window differs (CI04 only):** the **10–43 Hz** bands are elevated
(z ≈ +3.8 to +4.7; ~10 dB above the control median: −77 vs −87 dB) and envelope kurtosis
is high (the raw z=1163 is a MAD-collapse artifact — control kurtosis is very stable — driven
by a **single** broadband click at t≈10 s, visible in the spectrogram).

**Any clicks / anomalies?** No excess. The event window has **25 clicks — fewer than the
control window's 27**, and far below the most-anomalous window's 68. The one broadband click
is ordinary. Mid/high bands are normal.

**Is the event special?** No. It ranks #20/144 — inside the top 14%, but **19 unrelated
windows are more anomalous**, several driven by whale/click bouts and ship passages. Direct
comparison (see figure): the event window looks like a **distant vessel** (broadband
low-frequency rise) plus one click; the most-anomalous window (08:46 UTC, 7 features >3σ) is
a biosonar click bout with *lower* LF energy. And the two independent stations (CI01, CI05)
show **nothing** at the event time.

![CI04 comparison](figures/fig_case_control_ci04.png)
![Per-band z-scores](figures/fig_case_control.png)

## Verdict
At the closest public mic (167 km), the USS Omaha event window is **acoustically
unremarkable**: a modest 10–43 Hz rise consistent with a distant ship, one ordinary click,
normal click rate, no mid/high-band anomaly, statistically comparable to routine windows and
less anomalous than several unrelated ones, with zero corroboration at two independent
stations. The unbiased search finds **no distinctive signature** at the event.

**Caveat (unchanged):** 167 km is far; a genuine localized water-entry signal from a ~2 m
object would likely be below the noise here regardless. This bounds — it does not disprove.
The value is the method: an assumption-free, interpretable, cross-station case-control test
that says exactly which bands differ and by how much, ready to run the instant data from a
*close* sensor becomes available.
