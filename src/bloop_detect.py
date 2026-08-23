#!/usr/bin/env python3
"""Detect and characterize 'The Bloop' from NOAA's released recording.

The Bloop: an ultra-powerful low-frequency sound recorded in summer 1997 by
NOAA's Equatorial Pacific autonomous hydrophone array, source triangulated to
~50 S 100 W. By ~2012 NOAA attributed it to ICEQUAKES — large icebergs
cracking/calving near Antarctica (non-tectonic cryoseism) — not a marine
animal. NOAA's public clip is time-compressed 16x (that is why a low-frequency
event becomes audible).

This downloads NOAA's own bloop.wav and iceberg-calving reference, corrects the
axes back to true units (freq/speed, time*speed), and reproduces the Bloop's
signature: a discrete low-frequency (~1-40 Hz true) transient with a rising
'upsweep' spectral ridge. Same envelope/spectrogram tooling used elsewhere in
this repo, pointed at a famous historical event.
"""
import os
import urllib.request
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "https://www.pmel.noaa.gov/acoustics/sounds/"
CLIPS = [("bloop.wav", 16, "The Bloop (1997)"),
         ("A53JD33Calving_3x.wav", 3, "Iceberg A53 calving (NOAA reference)")]
CACHE = os.environ.get("SCRATCH", ".")


def get(fn):
    p = os.path.join(CACHE, fn)
    if not os.path.exists(p):
        req = urllib.request.Request(BASE + fn, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r, open(p, "wb") as f:
            f.write(r.read())
    return p


def analyze(fn, speed, label):
    x, fs = sf.read(get(fn))
    if x.ndim > 1:
        x = x[:, 0]
    x = x - np.mean(x)
    f, t, S = signal.spectrogram(x, fs=fs, nperseg=1024, noverlap=896)
    Sdb = 10 * np.log10(S + 1e-12)
    ft, tt = f / speed, t * speed            # undo the speed-up -> true units
    band = (ft >= 1) & (ft <= 120)
    ridge = ft[band][np.argmax(Sdb[band, :], axis=0)]
    env = np.sqrt(np.mean(S, axis=0))
    return dict(ft=ft, tt=tt, Sdb=Sdb, ridge=ridge, env=env, speed=speed,
                label=label, dur=len(x) / fs * speed)


def main():
    d = [analyze(*c) for c in CLIPS]
    b = d[0]
    ipk = int(np.argmax(b["env"]))
    core = (b["tt"] > b["tt"][ipk] - 8) & (b["tt"] < b["tt"][ipk] + 8) & np.isfinite(b["ridge"])
    slope = np.polyfit(b["tt"][core], b["ridge"][core], 1)[0] if core.sum() > 3 else float("nan")
    print("THE BLOOP — detected & characterized (true units):")
    print(f"  clip true duration {b['dur']:.0f} s; loudest at t≈{b['tt'][ipk]:.0f} s")
    print(f"  energy concentrated ~{np.percentile(b['ridge'],5):.0f}-"
          f"{np.percentile(b['ridge'],90):.0f} Hz (true) -> infrasonic/low-freq")
    print(f"  upsweep ridge slope in the core event: {slope:+.2f} Hz/s (rising = upsweep)")

    fig, ax = plt.subplots(2, 1, figsize=(13, 9))
    for a, dd, fmax in [(ax[0], d[0], 120), (ax[1], d[1], 120)]:
        a.pcolormesh(dd["tt"], dd["ft"], dd["Sdb"], shading="auto", cmap="magma",
                     vmin=np.percentile(dd["Sdb"], 55), vmax=np.percentile(dd["Sdb"], 99.7))
        a.plot(dd["tt"], dd["ridge"], color="cyan", lw=0.8, alpha=0.7, label="spectral ridge")
        a.set_ylim(1, fmax); a.set_ylabel("TRUE frequency (Hz)")
        a.set_title(f"{dd['label']}  (released clip {dd['speed']}x-sped; axes corrected to true units)")
        a.legend(loc="upper right", fontsize=8)
    ax[1].set_xlabel("TRUE time (s)")
    fig.suptitle("Detecting 'The Bloop': low-frequency upsweep vs NOAA iceberg-calving reference",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig("figures/fig_bloop.png", dpi=115)
    print("wrote figures/fig_bloop.png")


if __name__ == "__main__":
    main()
