# "Egg" craft — engravings on the hull

The marks are **sunken channels**. Each one has two raised rims, and under raking light one
rim is lit while the other is shadowed. The engraved figure is the channel *between* the
rims — not the rims themselves.

| file | what it is |
|---|---|
| `01_engravings_highlighted.png` | the hull with every incised figure highlighted, rims shown faintly behind |
| `02_engraving_plate.png` | per group: enhanced photo, the two rims, the figure between them, with real size on the source photograph |
| `04_sign_forms.png` | the figures alone at equal height — the sheet to hold against candidate scripts |
| `03_engravings_facsimile.png` | the figures in hull position and scale |
| `05_check_vs_sketches.png` | each figure beside the hand sketch of the same mark |
| `engravings.csv` | per group: confidence, part count, position, size, hull-normalised polar coordinates |
| `scripts/engravings.py` | reproduces panels 1–4 |
| `scripts/channel.py` | the rim-pairing operator |
| `scripts/validate_user.py` | reproduces panel 5 and the agreement numbers |
| `annotation/` | the original annotator's claim of 13 hieratic numerals, kept for reference |

## Three passes, two corrections

**1. Shadow only.** A dark-line ridge filter. Wrong: it followed the shadowed rim and threw
away the lit one — half of every mark, offset from where the mark is.

**2. Both rims.** Ridge filter at both polarities. Better, but still wrong in kind: it traced
the *edges* of each channel and called them the mark. The two traces looked alike because
they are the same feature lit from two sides.

**3. The channel between the rims.** Rims are detected at both polarities and
then **paired using the light direction**: a pixel belongs to the figure when a shadowed rim
lies behind it and a lit rim ahead of it along the illumination axis, closer together than
22 px. Pairing runs over a fan of directions (±45°, 15° steps) so channels of any orientation
are recovered, not only those square to the light — a single scan direction misses any stroke
running along the light axis, which is what cut the base off the antler in testing.

The ordering matters and is what makes the pairing unambiguous. Light is from the left, so
crossing a sunken channel left-to-right you meet the shadowed rim first and the lit rim
second. On the flat *between* two channels you meet them in the opposite order. Keying on the
order rejects the gaps and keeps the channels.

**4. Boundaries respected** — this version. A deep dark stroke is where the surface drops: it
is the *edge* of the incision, and nothing may be filled across it. Pass 3 ignored that and
ran its fill straight over the dark strokes, producing blobs. Now:

* **cut edges** are thresholded at the stroke's own scale (σ 3 smoothing first — a
  pixel-level threshold picks up nothing but screen speckle, because the darkest individual
  pixels are noise). A dark stroke is kept only where a lit rim lies nearby, which is what
  separates a cut edge from a stain.
* **the interior** is the pass-3 channel, but cut at the edges and required to touch one, so
  a fill can never cross an edge into the neighbouring mark.

Recall of the hand-sketched areas across the four passes: **0.23 → 0.40 → 0.64 → 0.64**, with
15/15 marks found throughout. Pass 4 does not add area over pass 3 — it puts the same area on
the right side of the boundaries.

## What came out

18 figures: 6 strong, 10 moderate, 2 weak.

**Where this stops.** The trace now follows the cut edges faithfully and respects them. It
does **not** reliably join those edges into the figure a reader sees: at the antler
(group 12) the edges come out as separate strokes rather than one connected three-tined form.
Joining them is a gestalt judgement, and at a median 25 × 31 px per mark there is not enough
in the pixels to make it a measurement. The `cut edges alone` pane of the plate is the most
neutral artefact for that work — it is what the surface actually shows, with the grouping
left open.

They sit in a band around the visible limb (`r_norm` mean 0.88, sd 0.08), the camera-facing
centre of the craft empty. That is where raking light makes relief legible, so it says where
marks are *readable*, not where they are. One large patch is unrecoverable: the annotator's
chart was pasted over it in the only good copy of this photograph.

## Before comparing against a script

Excluding weak groups, the median figure is about **25 × 31 px** on the source photograph. Part counts, relative angles and topology survive at that scale; fine
terminal shapes and curvature do not. `engravings.csv` gives `photo_w`/`photo_h` per group.

The hand sketches were used to test *detection* — did we find the mark, and roughly its extent
— never to tune the shape. Fitting them directly produces blobs, which an earlier pass here
did before that was caught.
