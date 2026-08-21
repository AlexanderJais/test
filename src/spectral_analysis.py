#!/usr/bin/env python3
"""Spectral survey of the MARS 256 kHz hydrophone slices.

For every WAV slice produced by ``fetch_data.py`` this script:

* converts samples to calibrated sound pressure (µPa) using the published
  MARS chain: 24-bit full scale = 3 V peak, icListen HF sensitivity
  ~ -177.9 dB re V/µPa (flat approximation; the measured curve varies
  -177..-181 dB across the band),
* computes a ~1 Hz resolution Welch PSD over the full 10 Hz - 128 kHz band,
* detects narrowband spectral "lines" (the only signal class that could
  even superficially resemble a coherent oscillation) as peaks >= 6 dB
  above a median-smoothed baseline,
* writes per-slice PSDs (.npz) + detected lines (.json) and a 6-panel
  spectrogram figure.

Usage:
    python3 src/spectral_analysis.py [--slices DIR] [--out results]
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
from matplotlib.colors import Normalize

VOLTS_PEAK = 3.0                 # full-scale = 3 V peak (file metadata)
SENS_DB_RE_V_PER_UPA = -177.9    # icListen HF #1689, flat approximation
V_TO_UPA = 10 ** (-SENS_DB_RE_V_PER_UPA / 20)  # ~7.85e8 µPa per volt

NPERSEG_PSD = 2 ** 18            # 0.977 Hz resolution at fs=256 kHz
LINE_SNR_DB = 6.0                # detection threshold above local baseline
FMIN = 10.0                      # ignore lowest bins (window leakage, DC)


def load_pressure(path):
    x, fs = sf.read(path, dtype="float64", always_2d=False)
    return x * VOLTS_PEAK * V_TO_UPA, fs


def welch_psd(p, fs):
    f, pxx = signal.welch(p, fs=fs, window="hann", nperseg=NPERSEG_PSD,
                          noverlap=NPERSEG_PSD // 2, detrend="constant")
    return f, pxx


def detect_lines(f, pxx):
    """Peaks >= LINE_SNR_DB above a 201-bin median baseline."""
    sel = f >= FMIN
    fq, p = f[sel], pxx[sel]
    baseline = signal.medfilt(p, kernel_size=201)
    snr_db = 10 * np.log10(p / np.maximum(baseline, 1e-30))
    peaks, props = signal.find_peaks(snr_db, height=LINE_SNR_DB, distance=5)
    return [
        {
            "freq_hz": round(float(fq[i]), 3),
            "snr_db": round(float(snr_db[i]), 2),
            "level_db_re_upa2hz": round(float(10 * np.log10(p[i])), 2),
        }
        for i in peaks
    ]


def spectrogram_panel(ax, p, fs, title):
    f, t, sxx = signal.spectrogram(p, fs=fs, window="hann", nperseg=8192,
                                   noverlap=4096)
    sxx_db = 10 * np.log10(np.maximum(sxx, 1e-30))
    lo, hi = np.percentile(sxx_db[f >= FMIN], [5, 99.8])
    m = ax.pcolormesh(t, f[1:], sxx_db[1:], shading="auto", cmap="magma",
                      norm=Normalize(lo, hi), rasterized=True)
    ax.set_yscale("log")
    ax.set_ylim(FMIN, fs / 2)
    ax.set_title(title, fontsize=9)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", default=os.environ.get(
        "PACIFIC_SOUND_DIR", "data/slices"))
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    paths = sorted(glob.glob(os.path.join(args.slices, "*.wav")))
    if not paths:
        raise SystemExit(f"no WAV slices found in {args.slices}")

    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharey=True)
    mesh = None
    for ax, path in zip(axes.ravel(), paths):
        name = os.path.basename(path).split("_first")[0]
        print(f"analyzing {name} ...")
        p, fs = load_pressure(path)
        f, pxx = welch_psd(p, fs)
        lines = detect_lines(f, pxx)
        np.savez_compressed(os.path.join(args.out, f"{name}_psd.npz"),
                            f=f.astype(np.float32), pxx=pxx.astype(np.float32))
        with open(os.path.join(args.out, f"{name}_lines.json"), "w") as fh:
            json.dump({"slice": name, "fs": fs, "n_lines": len(lines),
                       "lines": lines}, fh, indent=1)
        print(f"  {len(lines)} narrowband lines >= {LINE_SNR_DB} dB")
        mesh = spectrogram_panel(ax, p, fs, name)

    for ax in axes[:, 0]:
        ax.set_ylabel("frequency [Hz]")
    for ax in axes[1, :]:
        ax.set_xlabel("time [s]")
    fig.suptitle("MARS hydrophone (Monterey Bay, 891 m) — 60 s spectrograms, "
                 "10 Hz – 128 kHz", fontsize=12)
    cb = fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.01)
    cb.set_label("PSD [dB re 1 µPa²/Hz]")
    fig.savefig(os.path.join(args.figdir, "fig_spectrograms.png"), dpi=150)
    print(f"wrote {args.figdir}/fig_spectrograms.png")


if __name__ == "__main__":
    main()
