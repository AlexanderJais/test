#!/usr/bin/env python3
"""Zoom in on the standout transients from the anomaly sweep.

For each (file, UTC, label) target: read a 30 s window, draw waveform +
zoomed spectrogram, auto-locate the single loudest ~2 s and draw a fine view,
and export a normalized wav clip. Prints simple descriptors (peak crest,
kurtosis, dominant band, and whether it looks tonal/vessel-like).
"""
import os
import datetime
import re
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("sound", exist_ok=True)
os.makedirs("figures/transients", exist_ok=True)

CI01E = "data/omaha/SanctSound_CI01_02_671883305_190716033643.flac"
CI04A = "data/omaha/SanctSound_CI04_02_671924265_190715235545.flac"
CI04B = "data/omaha/SanctSound_CI04_02_671924265_190716055543.flac"
CI05E = "data/omaha/SanctSound_CI05_02_671359013_190716014401.flac"

TARGETS = [
    (CI05E, "05:30:31", "CI05_kurt963_0530"),   # highest kurtosis anywhere
    (CI04B, "06:05:43", "CI04_kurt193_0605"),   # most impulsive, on the event
    (CI05E, "02:51:15", "CI05_streak_0251"),    # coincident CI04/CI05 streak
    (CI05E, "04:00:01", "CI05_ship_0400"),      # middle of the ship-passage blob
    (CI01E, "06:01:45", "CI01_atevent_0601"),   # CI01 at the event, for contrast
]
WIN_S = 30.0


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1) + m.group(2), "%y%m%d%H%M%S")


def load(path, hhmmss, win=WIN_S):
    f0 = fstart(path)
    tt = datetime.datetime.strptime(hhmmss, "%H:%M:%S").time()
    target = datetime.datetime.combine(f0.date(), tt)
    if target < f0:
        target += datetime.timedelta(days=1)
    off = (target - f0).total_seconds()
    info = sf.info(path)
    fs = info.samplerate
    with sf.SoundFile(path) as fh:
        fh.seek(int(off * fs))
        x = fh.read(int(win * fs), dtype="float64")
    if x.ndim > 1:
        x = x[:, 0]
    return x - x.mean(), fs, target


def descr(x, fs):
    xh = signal.sosfiltfilt(signal.butter(2, 5, "high", fs=fs, output="sos"), x)
    rms = np.sqrt(np.mean(xh ** 2)) + 1e-20
    crest = np.max(np.abs(xh)) / rms
    from scipy import stats
    kurt = stats.kurtosis(xh)
    ff, pxx = signal.welch(xh, fs=fs, nperseg=4096)
    # band split
    bands = {"low<0.5k": (5, 500), "mid0.5-3k": (500, 3000), "hi3-20k": (3000, 20000)}
    bp = {n: pxx[(ff >= a) & (ff < b)].sum() for n, (a, b) in bands.items()}
    dom = max(bp, key=bp.get)
    # tonality: peak above median-smoothed
    pdb = 10 * np.log10(pxx + 1e-20)
    tonal = np.percentile(pdb - signal.medfilt(pdb, 51), 99.5)
    peakf = ff[np.argmax(pxx)]
    return crest, kurt, dom, tonal, peakf


def main():
    for path, hhmmss, label in TARGETS:
        x, fs, tutc = load(path, hhmmss)
        crest, kurt, dom, tonal, peakf = descr(x, fs)
        st = re.search(r"_(CI\d\d)_", path).group(1)
        print(f"{label:24s} {st} {tutc:%H:%M:%S}  crest={crest:6.1f} kurt={kurt:8.1f} "
              f"dom={dom:9s} tonal={tonal:4.1f}dB peakf={peakf:6.0f}Hz")

        # locate loudest 2 s
        env = np.abs(signal.hilbert(
            signal.sosfiltfilt(signal.butter(2, 200, "high", fs=fs, output="sos"), x)))
        w = int(2 * fs)
        csum = np.convolve(env, np.ones(w), "valid")
        i0 = int(np.argmax(csum))
        xz = x[i0:i0 + w]

        fig, ax = plt.subplots(3, 1, figsize=(12, 9))
        t = np.arange(len(x)) / fs
        ax[0].plot(t, x, lw=0.3, color="k")
        ax[0].axvspan(i0 / fs, (i0 + w) / fs, color="tab:red", alpha=0.15)
        ax[0].set_title(f"{label}  {st} {tutc:%H:%M:%S} UTC  "
                        f"(crest {crest:.0f}, kurt {kurt:.0f}, dom {dom})")
        ax[0].set_ylabel("amp"); ax[0].set_xlabel("s (30 s window)")
        f2, tt2, Sxx = signal.spectrogram(x, fs=fs, nperseg=2048, noverlap=1536)
        ax[1].pcolormesh(tt2, f2, 10 * np.log10(Sxx + 1e-20), shading="auto",
                         cmap="magma",
                         vmin=np.percentile(10*np.log10(Sxx+1e-20), 40),
                         vmax=np.percentile(10*np.log10(Sxx+1e-20), 99.5))
        ax[1].set_yscale("log"); ax[1].set_ylim(20, fs/2); ax[1].set_ylabel("Hz (30 s)")
        f3, tt3, Sz = signal.spectrogram(xz, fs=fs, nperseg=512, noverlap=448)
        ax[2].pcolormesh(tt3, f3, 10 * np.log10(Sz + 1e-20), shading="auto",
                         cmap="magma",
                         vmin=np.percentile(10*np.log10(Sz+1e-20), 40),
                         vmax=np.percentile(10*np.log10(Sz+1e-20), 99.5))
        ax[2].set_yscale("log"); ax[2].set_ylim(20, fs/2)
        ax[2].set_ylabel("Hz (loudest 2 s)"); ax[2].set_xlabel("s")
        fig.tight_layout()
        fig.savefig(f"figures/transients/{label}.png", dpi=110)
        plt.close(fig)

        # wav clip: loudest 6 s, normalized
        w6 = int(6 * fs)
        j0 = max(0, min(i0 - w6 // 3, len(x) - w6))
        clip = x[j0:j0 + w6]
        clip = clip / (np.max(np.abs(clip)) + 1e-12) * 0.95
        sf.write(f"sound/{label}.wav", clip.astype(np.float32), fs)
    print("\nwrote figures/transients/*.png and sound/*.wav")


if __name__ == "__main__":
    main()
