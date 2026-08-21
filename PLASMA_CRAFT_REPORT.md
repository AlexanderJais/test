# Can public underwater arrays find and locate a plasma-propelled underwater craft?

**Investigation report, part 2 — August 2026** · builds on [REPORT.md](REPORT.md)

## TL;DR

**Detection: yes, in principle — plasma propulsion would make a craft *easier* to hear, not harder.** Every physically credible plasma/electromagnetic propulsion concept is acoustically loud or statistically distinctive (§1). The discriminant is a craft-like source that is *sustained and moving* but shows **no propeller blade-rate modulation** while carrying anomalous broadband or impulse-train content.

**Location: yes, demonstrated with real public data.** Using the public OOI cabled hydrophone array (4 usable stations, 3–294 km baselines), this investigation localized a ground-truth impulsive underwater source — the T-phase of a USGS-reviewed M5.0 earthquake — to **25 km of its true position at ~240 km range**, using classic passive-sonar TDOA processing (§2). The same pipeline applies to any sufficiently loud craft.

**Did we find one? No.** A signature screen of 519 array windows classified everything as ambient, geophysical, biologic (verified fin whale song), vessel, or instrument glitch — zero anomalous-propulsion candidates (§3). And the honest caveat: public scientific arrays are sparse and band-limited compared to the purpose-built naval systems (SOSUS/IUSS) that do this for a living (§4).

---

## 1. Signature model: what would plasma propulsion sound like?

| Concept | Real-world anchor | Acoustic observables |
|---|---|---|
| **MHD seawater drive** ("caterpillar" drive: J×B force on conducted seawater, no moving parts) | Yamato-1, 1992 — functioned, ~8 kn | **No blade-rate lines** (nothing rotates in the water); electrolysis at the electrodes sheds gas bubbles → sustained broadband hiss; power converters at MW scale → harmonic comb tonals at switching frequencies |
| **Supercavitating vehicle** (gas/vapor cavity envelops the hull) | Shkval-class torpedoes (~200 kn class) | Extreme sustained broadband cavitation roar + propulsor noise — among the loudest man-made sources in the ocean; trivially detectable at long range |
| **Plasma sheath** (ionized layer around the hull — the user's scenario) | Lab scale only: electric discharges in water, plasma drag-reduction actuators (in air) | Liquid water quenches plasma, so a sheath must be sustained by **continuous electrical discharge inside a vapor layer**: repetitive micro-discharge shock pulses (a kHz-rate impulse train → heavy-tailed, high-kurtosis amplitude statistics), cavitation-like broadband, MW–GW power draw → boiling bubble wake. Also strong non-acoustic emissions (EM, thermal). |
| Conventional propeller craft (the thing to tell it apart from) | everything afloat | Cavitation broadband whose envelope is **modulated at shaft/blade rate** (DEMON analysis shows discrete modulation lines) + stable machinery tonals |

**Discriminant logic used in this work:** flag a source that is (a) sustained for minutes (not a single packet), (b) broadband-elevated above ambient, (c) impulsive (high kurtosis, many spikes/s), (d) **without** DEMON blade-rate lines, and (e) moving (drifting TDOAs across the array). §3 shows why each condition is needed — dropping any one of them lets earthquakes, whales, or glitches through.

## 2. Localization with a real public array — validated

**The array.** OOI Regional Cabled Array low-frequency hydrophones (network `OO`, channel `HDH`): 200 Hz sampling, response-calibrated to pascals, common GPS timebase — public via EarthScope FDSN web services (`src/fetch_array.py` downloads them with no credentials). Stations: AXBA1/AXCC1/AXEC2 (Axial Seamount cluster, 3–26 km spacing) + HYS14 (Hydrate Ridge, ~380 km east; its pair HYSB1 was offline on the validation day).

**Ground truth.** USGS-reviewed **M5.0 earthquake, 2023-01-11 10:17:18 UTC, 43.9646° N 128.7389° W** (Blanco region), 221–294 km from the stations. Its water-borne T-phase is a loud impulsive source with independently known position — exactly what a "craft transient" would look like to the array.

**Method** (`src/localize_tdoa.py`): band-pass 3–40 Hz → smoothed Hilbert envelopes → cross-correlate all 6 station pairs → TDOAs → hyperbolic grid search over position and water-borne speed.

**Results:**

- Pair correlations 0.93–0.99; TDOAs from −1.5 s (AXCC1↔AXEC2, 3 km apart) to +52.7 s (AXBA1↔HYS14).
- Best fit at group speed **c = 1.47 km/s** (textbook T-phase speed — recovered from the data, not assumed).
- **Acoustic fix 43.78° N, 128.92° W → 25.1 km from the USGS epicenter** (≈10% of range with a one-sided array — the misfit valley in `figures/fig_localization.png` shows the along-range elongation the geometry dictates).
- Independent check: envelope-peak times track predicted range/c moveout across all four stations to within −2.9…+8.0 s over 221–294 km paths (`fig_tphase.png`).

![Record section](figures/fig_tphase.png)
![Localization](figures/fig_localization.png)

**Implication:** any craft radiating comparably strong low-frequency energy inside or near the array's footprint is localizable with this exact pipeline; a moving craft yields TDOA tracks (a sequence of fixes) rather than one point.

## 3. Signature screening of the array data

`src/signature_features.py` computed, for every 30 s window of every station (519 windows across the event day and two comparison windows): broadband elevation over each station's own ambient baseline, amplitude kurtosis, spike-train count, DEMON modulation-line SNR, narrowband tonal count, and biologic pulse-train regularity — then applied the §1 discriminants.

| Class | Windows | Notes |
|---|---:|---|
| ambient | 450 | baseline ocean |
| geophysical energy burst (T-phase class) | 40 | the M5.0 packet + coda across stations |
| biologic pulse train | 12 | **verified fin whale 20 Hz song** — regular ~10 s inter-pulse intervals (CV as low as 0.006), 15–30 Hz downsweeps, January = peak fin whale season at Axial |
| vessel machinery (tonal comb) | 7 | shipping |
| isolated spike (glitch / single close transient) | 5 | single-sample events, kurtosis ~800 |
| broadband elevation, unclassified | 5 | HYS14 around the event window — consistent with a distant vessel near Hydrate Ridge; flagged honestly, not force-classified |
| **anomalous-propulsion candidate** | **0** | |

![Feature space](figures/fig_signatures.png)

**Lessons the screening taught (visible in the figure):** the naive criterion "loud + no blade lines" is *not* sufficient — T-phase windows land in exactly that region. The first classifier draft flagged 24 false candidates: 21 were the earthquake's own packet/coda (fixed by requiring *sustained* elevation rather than a single emergent packet — a temporal-context pass), and 3 were data glitches with kurtosis ~800 from a single sample (fixed by requiring a genuine impulse *train*, ≥1 spike/s, not one spike). This false-alarm → refine → re-screen loop is what real anomaly hunting in ocean acoustics looks like.

## 4. Honest capability assessment

What this demonstrates a **public** array can already do:
- hear strong low-frequency sources at basin scale (the M5 T-phase arrived with ~30 dB of headroom at 294 km),
- localize impulsive sources to tens of km at hundreds of km range, better inside the array footprint,
- classify sources by propulsion physics (blade lines vs none, sustained vs packet, pulse trains).

What it cannot do, and what finding a real plasma-propelled craft would take:
- **Band**: 200 Hz sampling caps analysis at 100 Hz. Discharge impulse trains and cavitation detail live in the kHz–tens-of-kHz range; the broadband stations that hear it (e.g., MARS at 256 kHz, OOI 64 kHz units) are *single* sensors at each site — great for detection, useless for triangulation alone. A purpose-built search wants dense broadband arrays.
- **Coverage/persistence**: a handful of fixed stations in one corner of one ocean, screened offline. Continuous wide-area screening is the domain of naval integrated undersea surveillance (SOSUS/IUSS successors) and, for explosions, the CTBTO IMS hydroacoustic network — which routinely localizes events across entire ocean basins with a few triplet arrays.
- **Quiet targets**: everything above assumes the craft radiates. The acoustic paradox of this question: every plasma-propulsion concept is loud — plasma discharge in water *is* essentially controlled cavitation. A craft that were acoustically silent would, by that very fact, not be plasma-propelled in any known sense; hunting it shifts to non-acoustic channels (magnetic anomaly, wake turbulence/thermal signature, bioluminescent wake, EM emissions — a plasma sheath would be a radio beacon).

**Bottom line:** yes — public array data plus standard passive-sonar processing (demonstrated end-to-end here, validated at 25 km accuracy against ground truth) would find and localize a plasma-propelled craft if one operated near the array at ordinary source levels; the discriminant signature (sustained, moving, broadband/impulsive, no blade lines) is well-defined and none of the surveyed data contains it.

## Reproduce

```bash
pip install -r requirements.txt
python3 src/fetch_array.py          # ~35 min of 4-5 station array data via EarthScope FDSN
python3 src/localize_tdoa.py        # TDOA fix vs USGS truth, fig_tphase, fig_localization
python3 src/signature_features.py   # 519-window screen, fig_signatures
```

Data credits: OOI Regional Cabled Array via EarthScope (IRIS) FDSN services; USGS earthquake catalog.
