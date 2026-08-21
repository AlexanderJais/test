#!/usr/bin/env python3
"""Broadband transient (impulse) detection — the nearest *real* link
between plasma physics and hydrophone data.

Electron plasma oscillations themselves cannot exist in seawater, but two
plasma-*generated* acoustic sources genuinely appear in ocean recordings:

* lightning strikes on the sea surface — the return stroke is a ~30,000 K
  plasma channel whose shock/thunder couples into the water column;
* snapping-shrimp cavitation collapse — the imploding bubble briefly forms
  a sonoluminescent microplasma; the audible "snap" is the collapse.

Both arrive as broadband impulses.  This script runs a simple matched
detector for that signal class on the survey slices: bandpass 5-40 kHz,
Hilbert envelope, robust (median + k*MAD) threshold, minimum-separation
peak picking.  It reports per-slice impulse rates and plots the strongest
event.

Localization note (the "locate it" half of the question): a single fixed
hydrophone gives arrival times but no bearing, so impulses can only be
*counted* here, not positioned.  Real localization uses time-difference-
of-arrival across arrays/networks, or cross-matching event times against
global lightning networks (WWLLN/GLD360).
"""

import argparse
import glob
import json
import os

import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spectral_analysis import load_pressure

BAND = (5000.0, 40000.0)
MAD_K = 12.0            # detection threshold: median + K * MAD of envelope
MIN_SEP_S = 0.005       # refractory time between detections


def detect_impulses(p, fs):
    sos = signal.butter(6, BAND, btype="bandpass", fs=fs, output="sos")
    x = signal.sosfiltfilt(sos, p)
    env = np.abs(signal.hilbert(x))
    med = np.median(env)
    mad = np.median(np.abs(env - med))
    thr = med + MAD_K * mad
    peaks, props = signal.find_peaks(env, height=thr,
                                     distance=int(MIN_SEP_S * fs))
    return x, env, thr, peaks, props["peak_heights"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", default=os.environ.get(
        "PACIFIC_SOUND_DIR", "data/slices"))
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    summary, best = [], None
    for path in sorted(glob.glob(os.path.join(args.slices, "*.wav"))):
        name = os.path.basename(path).split("_first")[0]
        p, fs = load_pressure(path)
        x, env, thr, peaks, heights = detect_impulses(p, fs)
        dur = len(p) / fs
        summary.append({"slice": name, "n_impulses": int(len(peaks)),
                        "rate_per_s": round(len(peaks) / dur, 2),
                        "threshold_upa": round(float(thr), 1)})
        print(f"{name}: {len(peaks)} impulses ({len(peaks)/dur:.1f}/s)")
        if len(peaks) and (best is None or heights.max() > best[0]):
            i = peaks[np.argmax(heights)]
            best = (heights.max(), name, fs, i,
                    x[max(0, i - fs // 100):i + fs // 100].copy())

    with open(os.path.join(args.out, "transients.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4),
                                   gridspec_kw={"width_ratios": [1, 1.4]})
    names = [s["slice"].replace("MARS_", "") for s in summary]
    ax1.bar(range(len(summary)), [s["rate_per_s"] for s in summary],
            color="tab:blue")
    ax1.set_xticks(range(len(summary)))
    ax1.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax1.set_ylabel("impulses per second")
    ax1.set_title(f"Broadband impulse rate ({BAND[0]/1e3:.0f}–"
                  f"{BAND[1]/1e3:.0f} kHz band)", fontsize=10)
    if best is not None:
        _, name, fs, i, w = best
        t = (np.arange(len(w)) - len(w) // 2) / fs * 1000
        ax2.plot(t, w / 1e6, lw=0.7, color="tab:red")
        ax2.set_xlabel("time re. peak [ms]")
        ax2.set_ylabel("band-passed pressure [Pa]")
        ax2.set_title(f"strongest impulse — {name} (t = {i/fs:.2f} s)",
                      fontsize=10)
    fig.suptitle("Impulsive transients: the signal class plasma-generated "
                 "sound (lightning, snapping-shrimp sonoluminescence) "
                 "belongs to", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_transients.png"), dpi=150)
    print(f"wrote {args.figdir}/fig_transients.png")


if __name__ == "__main__":
    main()
