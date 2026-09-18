# "Egg" craft — hieratic numerals mapped onto the hull

Recreation of the annotated image, with the 13 claimed hieratic numerals highlighted
*on* the ship instead of being listed in a side chart with arrows.

| file | what it is |
|---|---|
| `01_ship_numerals_highlighted.png` | main panel — same framing as the source annotation, surface enhanced, 13 ringed sites each labelled with its sign and value |
| `02_site_details.png` | each site zoomed, next to the sign it is claimed to be |
| `03_context_whole_craft.png` | the whole craft, showing where the 13 sites sit and which region the main panel covers |
| `positions.csv` | the 13 positions in four coordinate systems + hull-normalised polar coordinates |
| `analysis_meta.json` | registration transform, hull ellipse, provenance |
| `scripts/pipeline.py` | reproduces everything from `source/` |
| `source/` | the two input images |

## What was done

**1. Read the annotation.** The 13 red arrows were located by their heads (morphological
opening isolates the solid head from the thin shaft). Each shaft was then traced back to
its origin: the ray from the head collecting the most red pixels within 2.2 px, walked to
the far end of its contiguous red run. The origin was matched against the numeral chart's
grid, recovered from the table rules at `x = 42/155/165/277/287/429/439/581` and
`y = 565 … 1360`.

This matters — a naive step-along-the-ray walk mis-read three of the arrows where they
cross other arrows. 12 of the 13 traced origins land squarely on a glyph; the exception
is #2, whose stroke runs over the `10` glyph and stops ~55 px past it on the cell corner
(`10` is the only sign it touches, and the two neighbouring cells lie in the opposite
direction from the stroke).

**2. Registered the annotated image to the photo.** The annotated image is *not* a crop of
the posted screenshot — it is a ~3.9× view of the same photograph from a better source,
and it carries real detail the screenshot has lost. Silhouette fitting was useless here
(the visible arc is only part of a very large ellipse, so the fit is under-constrained and
disagreed with itself by 20°). What worked was multi-scale template matching on the cave
rock, which has real texture: **scale 0.25889, rotation −2.318°, translation (519.29,
389.93)**, from 17 correspondences spanning a 1400 px baseline, 1.8 px rms. The annotated
frame covers photo panel `x 519…920, y 390…892`.

**3. Rebuilt the image.** The red arrows were erased using an unbiased local mean
(normalised convolution over unmasked pixels) plus high-frequency texture donated from
26–80 px away, so the repaired strokes keep the surrounding surface grain rather than
leaving flat bands. The hull was then enhanced with a band-pass (σ 2.6 / 34) and mild
local contrast, applied strongly inside the hull and gently outside.

## The 13 readings

| # | value | # | value | # | value |
|---|---|---|---|---|---|
| 1 | 3 | 6 | 3 | 11 | 100 |
| 2 | 10 | 7 | 5 | 12 | 4 |
| 3 | 50 | 8 | 900 | 13 | 80 |
| 4 | 800 | 9 | 10 | | |
| 5 | 500 | 10 | 70 | | |

Sum 2535. `3` and `10` each appear twice; no thousands are used.

## For the position analysis

`positions.csv` carries, besides pixel coordinates, hull-normalised polar coordinates
against the fitted hull ellipse (centre 575.9, 705.2 px; semi-axes 300.3 × 246.8 px;
54.8°):

* `u_major`, `v_minor` — position along the major / minor axis, in units of the semi-axis
* `r_norm` — 1.0 is the fitted rim
* `theta_deg` — angle from the major axis, image y down

Two things are already visible and worth testing properly:

* all 13 sit on the **upper-right quadrant** — the region the annotator zoomed into, so
  this says more about the annotation than the craft;
* they lie in a **narrow band near the rim**: `r_norm` mean 0.82, sd 0.10, and ordering
  them by `theta_deg` walks them around the hull in sequence. Whether that is structure on
  the hull or just an artefact of marks only being legible at grazing incidence near the
  limb is the thing to settle next.

## Caveats

* The rings mark the *annotated* positions. The surface there does carry real, repeatable
  streak-like marks, but at the resolution available most of them do not resolve into
  anything that unambiguously reads as the claimed sign — `02_site_details.png` puts each
  one next to its claimed glyph so the resemblance can be judged directly.
* The glyph images in the labels are cropped from the annotator's own chart, not from the
  hull.
* Per-site contrast stretching was deliberately avoided: at this resolution it amplifies
  the screen half-tone of the photographed display into convincing-looking shapes. A
  single uniform enhancement is used throughout.
