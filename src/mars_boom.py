#!/usr/bin/env python3
"""Boom/infrasound search on live MARS audio around a dated fireball.

Event: green fireball over central California, ~2026-03-22 20:19 PDT =
2026-03-23 03:19:00 UTC, hundreds of witnesses (NASA/AMS). MARS (Monterey
Bay, 36.713N 122.186W, 256 kHz, ~890 m depth) is live and has the audio.

A bolide airburst radiates INFRASOUND (~0.1-20 Hz) that travels the air path
(~0.34 km/s) and couples through the sea surface -> arrives at MARS MINUTES
after the visual flash, delay = ground-range/0.34 km/s. So we scan the window
AFTER 03:19 for an impulsive low-frequency arrival.

Method (per 10-min MARS file, streamed then deleted to save disk):
  * read in blocks, decimate 256 kHz -> 250 Hz (Nyquist 125 Hz),
  * concatenate to a continuous low-rate series across the window,
  * infrasound band (1-20 Hz) envelope -> transient detector (MAD threshold),
  * low band (20-80 Hz) for corroboration,
  * spectrogram + band time series, fireball flash + expected-arrival band
    marked.
Honest: a high-altitude green fireball far from Monterey may leave NO trace.
This is a search; it reports whatever is (or isn't) there.
"""
import os
import sys
import datetime
import urllib.request
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUCKET = "https://pacific-sound-256khz-2026.s3.amazonaws.com"
FILES = ["03/MARS_20260323_031000.wav",   # 03:10-03:20 (pre / flash at 03:19)
         "03/MARS_20260323_032000.wav",   # 03:20-03:30
         "03/MARS_20260323_033000.wav",   # 03:30-03:40
         "03/MARS_20260323_034000.wav"]   # 03:40-03:50
FLASH_UTC = datetime.datetime(2026, 3, 23, 3, 19, 0)
TARGET_FS = 250
SCRATCH = os.environ.get("SCRATCH", ".")


def file_start(key):
    m = key.split("/")[-1]              # MARS_YYYYMMDD_HHMMSS.wav
    return datetime.datetime.strptime(m[5:20], "%Y%m%d_%H%M%S")


def decimate_to(x, fs, target):
    # staged FIR decimation 256000 -> 250 (factors 8*8*16=1024)
    for q in (8, 8, 16):
        x = signal.decimate(x, q, ftype="fir", zero_phase=True)
        fs //= q
    return x, fs


def fetch_and_decimate(key):
    url = f"{BUCKET}/{key}"
    local = os.path.join(SCRATCH, key.split("/")[-1])
    print(f"  downloading {key} ...", flush=True)
    urllib.request.urlretrieve(url, local)
    info = sf.info(local)
    fs = info.samplerate
    blk = fs * 10                       # 10 s blocks
    low = []
    with sf.SoundFile(local) as fh:
        while True:
            x = fh.read(blk, dtype="float64")
            if len(x) == 0:
                break
            if x.ndim > 1:
                x = x[:, 0]
            if len(x) < 400:
                break
            xd, fsd = decimate_to(x, fs, TARGET_FS)
            low.append(xd)
    os.remove(local)                    # free disk immediately
    return np.concatenate(low), TARGET_FS, file_start(key)


def main():
    segs, t0 = [], None
    for k in FILES:
        try:
            xd, fsd, fst = fetch_and_decimate(k)
        except Exception as e:
            print(f"  skip {k}: {e}", flush=True)
            continue
        if t0 is None:
            t0 = fst
        segs.append((fst, xd))
    if not segs:
        print("no data"); return
    fs = TARGET_FS
    # assemble continuous series on an absolute-time grid
    total = int((segs[-1][0] - t0).total_seconds() * fs) + len(segs[-1][1])
    sig = np.full(total, np.nan)
    for fst, xd in segs:
        i0 = int((fst - t0).total_seconds() * fs)
        sig[i0:i0 + len(xd)] = xd
    good = np.isfinite(sig)
    sig[~good] = 0.0
    t_utc = [t0 + datetime.timedelta(seconds=i / fs) for i in range(total)]
    flash_s = (FLASH_UTC - t0).total_seconds()

    def band_env(lo, hi):
        sos = signal.butter(4, (lo, hi), btype="bandpass", fs=fs, output="sos")
        return np.abs(signal.hilbert(signal.sosfiltfilt(sos, sig)))

    infra = band_env(1, 20)
    lowb = band_env(20, 80)
    # transient detector on infrasound envelope (1 s smoothing)
    w = fs
    sm = np.convolve(infra, np.ones(w) / w, "same")
    med = np.median(sm); mad = np.median(np.abs(sm - med)) + 1e-30
    z = 0.6745 * (sm - med) / mad
    thr = 8
    pk, _ = signal.find_peaks(z, height=thr, distance=fs * 5)
    print(f"\nwindow {t0:%H:%M}-{t_utc[-1]:%H:%M} UTC, flash at {FLASH_UTC:%H:%M:%S}")
    print(f"infrasound (1-20Hz) transients over MAD z>={thr}:")
    if len(pk) == 0:
        print("  none")
    for i in pk:
        dt = (t_utc[i] - FLASH_UTC).total_seconds()
        rng = dt * 0.34  # implied air-path range if this were the bolide
        print(f"  {t_utc[i]:%H:%M:%S} UTC  z={z[i]:.1f}  ({dt:+.0f}s after flash"
              f" -> implied ~{rng:.0f} km if bolide)")

    # figure
    tm = (np.arange(total) / fs) - flash_s   # seconds relative to flash
    fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    f2, tt, Sxx = signal.spectrogram(sig, fs=fs, nperseg=1024, noverlap=768)
    ax[0].pcolormesh(tt - flash_s, f2, 10 * np.log10(Sxx + 1e-20), shading="auto",
                     cmap="magma", vmin=np.percentile(10*np.log10(Sxx+1e-20), 30),
                     vmax=np.percentile(10*np.log10(Sxx+1e-20), 99))
    ax[0].set_ylim(0, 125); ax[0].set_ylabel("Hz")
    ax[0].set_title("MARS 2026-03-23 ~03:19 UTC fireball window — decimated to 250 Hz")
    ax[1].plot(tm, infra, lw=0.6, color="tab:blue", label="infrasound 1-20 Hz env")
    ax[1].plot(tm, lowb, lw=0.6, color="tab:orange", alpha=0.7, label="low 20-80 Hz env")
    for a in ax:
        a.axvline(0, color="red", ls="--", lw=2)
        a.axvspan(0, 900, color="yellow", alpha=0.08)  # plausible arrival (0-300km)
    ax[1].set_xlabel("seconds relative to fireball flash (red); yellow = 0-300 km arrival band")
    ax[1].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig("figures/fig_mars_fireball_boom.png", dpi=115)
    print("\nwrote figures/fig_mars_fireball_boom.png")


if __name__ == "__main__":
    main()
