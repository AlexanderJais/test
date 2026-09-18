# "Egg" craft — engravings on the hull

Two passes over the same two source images:

* **`/` (main)** — the engravings themselves, detected and traced, with no script assumed.
* **`annotation/`** — a record of what the original annotator claimed (13 hieratic numerals).
  Kept for reference; the traced shapes do **not** match those glyphs.

| file | what it is |
|---|---|
| `01_engravings_highlighted.png` | the hull with every detected mark group highlighted in place, colour-coded by confidence |
| `02_engraving_plate.png` | each group: enhanced photo beside the extracted shape, with its real size on the source photograph |
| `04_sign_forms.png` | the traced forms alone, normalised to equal height — the sheet to hold against candidate scripts |
| `03_engravings_facsimile.png` | tracing only, in hull position and scale |
| `engravings.csv` | per group: confidence, stroke count, position, size, hull-normalised polar coordinates |
| `scripts/engravings.py` | reproduces all four panels from `source/` |
| `scripts/pipeline.py` | the annotation pass (arrow tracing, registration, arrow removal, enhancement) |

## How the marks were found

The engravings are **dark grooves under raking light**, 4–10 px wide in the enhanced frame.
The competing signal is the photographed screen's scan-line striping, 1–2 px wide — so the
separation is done by *scale*, not by contrast:

1. arrows erased, hull enhanced (band-pass σ 2.6/34 + local contrast);
2. **Hessian ridge filter** (Frangi, dark-line polarity) at σ 2.5 / 3.5 / 4.5 — responds to
   grooves, largely ignores the thinner striping;
3. thresholds taken **relative to a local ridge-response floor**, because the striping
   raises that floor unevenly across the hull. This is what demotes the bottom-left corner,
   where the striping is worst, to `weak`;
4. hysteresis → geodesic reconstruction → oriented closing, so each groove emerges as one
   connected stroke instead of a dashed line;
5. strokes grouped, each group ranked by its contrast against the local floor
   (`rel_contrast`): **strong** ≥ 8, **moderate** 5–8, **weak** < 5.

## What came out

26 mark groups: 10 strong, 13 moderate, 3 weak. They are **not spread over the hull** — they
sit in a band around the visible limb, with the camera-facing centre of the craft empty.
That is the signature of relief being legible only where the light rakes across it, so it is
a lighting artefact of *visibility*, not evidence about where the marks actually are. A large
patch is also simply unreadable: the annotator's chart was pasted over it, and the underlying
pixels are gone from this source.

Forms that recur: single curved strokes, pairs and triplets of near-parallel strokes, and a
few compound figures — `#9` is a clear three-stroke form joined at the base, `#17` is a
comb-like group of parallel strokes, `#18` a set of parallel diagonals with a connector.

Not everything the detector returns is an engraving. `#6`, `#16` and `#22` are long, smooth,
featureless bands: they pass the elongation test but read more like a shading edge or a
surface crease than a cut groove. They are left in rather than quietly dropped — the
confidence column and the photo panel in `02_engraving_plate.png` are there so each one can
be judged directly.

## Before comparing against a script

These marks are **small**. Excluding the weak group, the median size on the source
photograph is about 14 × 16 px and the largest is 26 × 48 px. The annotated image resolves them better than the posted
screenshot because it comes from a better copy of the same photograph, but the underlying
information is still only tens of pixels per mark. Stroke *counts*, *relative angles* and
*grouping* are the properties that survive at that scale; fine terminal shapes and
curvature do not, and should not be used to argue for or against a particular script.

`engravings.csv` gives `photo_w`/`photo_h` for each group so this can be checked per mark
rather than taken on trust.
