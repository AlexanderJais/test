# "Egg" craft — engravings on the hull

Three passes over the same two source images:

* **`/` (main)** — the engravings traced as **relief**: each mark's shadow side *and* its lit
  side, since raking light splits every mark into both.
* **`annotation/`** — what the original annotator claimed (13 hieratic numerals). Kept for
  reference; the traced shapes do not match those glyphs.
* **`source/user_outlines_plate.webp`** — hand sketches supplied as a check.

| file | what it is |
|---|---|
| `01_engravings_highlighted.png` | the hull with every mark highlighted, shadow side and lit side in different colours |
| `02_engraving_plate.png` | per group: enhanced photo, relief traced on the photo, two-tone facsimile, with real size on the source photograph |
| `04_sign_forms.png` | the forms alone at equal height — the sheet to hold against candidate scripts |
| `03_engravings_facsimile.png` | two-tone tracing in hull position and scale |
| `05_check_vs_sketches.png` | each detection beside the hand sketch of the same mark |
| `engravings.csv` | per group: confidence, part count, position, size, raised/sunken reading, hull-normalised polar coordinates |
| `scripts/engravings.py` | reproduces panels 1–4 from `source/` |
| `scripts/validate_user.py` | reproduces panel 5 and the agreement numbers |

## The correction that produced this version

The first attempt used a dark-line ridge filter only. That was wrong: these marks are
**three-dimensional**, so raking light gives each one a shadow on one side and a lit face on
the other, about 4 px apart in this frame. Tracing only the dark ridge follows the shadow and
discards the lit half — roughly half of every mark, and offset from where the mark actually is.

Measured on the hand-sketched marks, the two halves carry the **same** signal: mean local
contrast 3.26 inside the sketches for the shadow channel, 3.05 for the highlight channel. The
highlight was never the weaker cue; it was simply not being looked at.

So the filter is now run at both polarities, σ 3/5/7/9, and the two channels are kept
**separate** through to the output. They are deliberately *not* merged into one silhouette:
morphological closing wide enough to join a shadow to its own highlight also fuses
neighbouring strokes, and the form of a multi-stroke sign is destroyed. Dark and pale in the
facsimile are the two faces of the same relief.

A component is only kept if it shows **both** polarities — that is what a real relief feature
looks like, and it is a strong filter against striping artefacts.

## Checking against the hand sketches

The supplied outlines are loose indications of the 3-D shape, not exact boundaries, and only
cover the marks that could be identified confidently. They are therefore used to test
*detection*, never to tune the shape — fitting them directly produces blobs, which is what a
first attempt at this did.

* **15 / 15** sketched marks are detected.
* **80 %** of what the detector finds inside the sketched zone falls within a sketch.
* Area recall against the sketches is 0.40, and is *expected* to be well below 1: the sketches
  are envelopes drawn around a mark, the tracing follows the relief itself.

Against the old shadow-only run, recall of the sketched area rose from 0.23 to 0.72 when the
detector was still producing filled envelopes — that number is the honest measure of how much
was being missed before.

## What came out

23 mark groups: 4 strong, 9 moderate, 10 weak. They sit in a band around the visible limb
(`r_norm` mean 0.86, sd 0.07), with the camera-facing centre of the craft empty — the
signature of relief being legible only where light rakes across it. This says where marks are
*readable*, not where they are. A large patch is simply unrecoverable: the annotator's chart
was pasted over it in the only good copy of this photograph.

**Raised or sunken** is not settled. Taking the hull's own shading, which puts the light to
the left, 14 groups read as sunken and 8 as raised, and the shadow→highlight direction is only
weakly concentrated (circular R = 0.34). If the light is in fact from the right, every reading
flips — the standard crater illusion, which one image cannot resolve. The `relief` column
records the reading under the left-light assumption. **The traced shape is the same either
way**, so this does not affect script comparison.

## Before comparing against a script

Excluding the weak groups, the median mark is about **20 × 26 px** on the source photograph
and the largest is 39 × 56. The annotated image resolves them better than the posted
screenshot because it comes from a better copy of the same photograph, but this is still tens
of pixels per mark. Part counts, relative angles and grouping survive at that scale; fine
terminal shapes and curvature do not. `engravings.csv` gives `photo_w`/`photo_h` per group so
this can be checked mark by mark.
