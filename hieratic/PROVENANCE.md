# Provenance of the source image

## What it is

The photograph analysed here comes from **4chan /x/ post No.39683660, "Another Egg UFO
photo", 21 January 2025** — the post number and timestamp are legible in the screenshot
(`source/post_screenshot.webp`).

It belongs to a set posted across three threads on 20–21 January 2025 (39677079, 39683084,
39683660) claiming to show an egg-shaped craft recovered from an Antarctic ice cave in 2022.
The set rode a wave of attention following a NewsNation segment with Jacob Barber, who
described an "egg-shaped" object recovered under a UAP retrieval programme.

## What is known about it

* **The poster described the material as a LARP.** An image in this series carried the
  caption *"This is a LARP. Maybe I will answer some questions. Everything is going to come
  to light soon anyways."*
* **Contemporary assessment was that the images are synthetic.** They were widely identified
  as AI-generated or game-engine renders; the most-cited analysis argued the cave surface
  textures were computer-generated — *"created in Unreal Engine with a filter over it"* — and
  that each image had a different background, which does not fit a set of photographs taken
  at one location.
* **An independent look at the carvings reached the same place this analysis did.** Sculptor
  Peter Osborne, examining enhanced versions, judged the marks *"rushed, almost amateurish"*,
  lacking the symmetry and coherence expected of either ancient or advanced work.

None of that is proof of how the images were made — no creator has been identified — but it
is the state of the public record, and it is the context any reading of the marks sits in.

## Attempt to find a better copy

Requested, and attempted. Results:

* **Higher-resolution copies of other images in the set were found** and are kept in
  `source/related/`: the five-frame cave strip at **1848 × 4987** (about 4× the width it has
  in the screenshot) and the "pod" frame at **1920 × 1080**. The two cave frames in that strip
  do show an egg-shaped object, but blown out to flat white — no surface detail at all.
* **No better copy of this particular cave-egg photograph was found.** The 4chan archive
  mirrors (desuarchive, archived.moe, thebarchive, archive.palanq.win) all refuse requests
  from here with HTTP 403.
* **Be wary of "enhanced" and "upscaled" versions in circulation** — at least one is archived
  publicly. An AI upscale invents exactly the kind of fine surface detail this analysis
  measures, so an upscaled copy is worse than useless here, however much sharper it looks.

## Is the annotated image itself an upscale?

Worth checking, since everything traced here depends on it: `source/annotated.webp` carries
detail the screenshot does not, and an AI upscale would be an alternative explanation for
that.

Its radial power spectrum rolls off smoothly from f = 0.05 to f = 0.90 of Nyquist with **no
cliff at f = 0.259**, which is where the screenshot's own Nyquist falls. Interpolating a
low-resolution image up leaves a sharp break there; there isn't one. So the annotated image
is a genuinely better copy of the same photograph rather than a naive upscale of the
screenshot. A *learned* upscaler cannot be excluded by this test — it would fill that band
smoothly too.
