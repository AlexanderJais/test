#!/usr/bin/env python3
"""Find a real 'boom in water': an earthquake T-phase on live MARS audio.

Ground truth (USGS): M6.0, 2026-01-16 03:25:53 UTC, epicenter 43.69 N
128.06 W, "off the coast of Oregon" — a large OCEANIC event (efficient
T-phase source). Range to MARS (36.713 N, 122.186 W) = 921 km. A T-phase
travels the SOFAR channel at ~1.48 km/s, so it should arrive at MARS about
921/1.48 = 622 s after origin -> ~03:36:15 UTC, as an emergent, spindle-
shaped broadband (~2-90 Hz) arrival tens of seconds long. This predicts the
time from independent seismology; the hydrophone either shows a boom there or
it does not.

Method: stream MARS 256 kHz files over 03:10-03:50, decimate to 250 Hz
(delete raw to save disk), build a broadband (2-90 Hz) envelope, find the
strongest transient, and compare its arrival to the seismic prediction.
"""
import os
import datetime
import urllib.request
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUCKET = "https://pacific-sound-256khz-2026.s3.amazonaws.com"
FILES = ["01/MARS_20260116_031000.wav",   # 03:10-03:20 baseline
         "01/MARS_20260116_032000.wav",   # 03:20-03:30
         "01/MARS_20260116_033000.wav",   # 03:30-03:40 (predicted arrival)
         "01/MARS_20260116_034000.wav"]   # 03:40-03:50
ORIGIN = datetime.datetime(2026, 1, 16, 3, 25, 53)
RANGE_KM = 921.0
C_TPHASE = 1.48
PRED_ARR = ORIGIN + datetime.timedelta(seconds=RANGE_KM / C_TPHASE)
TARGET_FS = 250
SCRATCH = os.environ.get("SCRATCH", ".")


def fstart(key):
    m = key.split("/")[-1]
    return datetime.datetime.strptime(m[5:20], "%Y%m%d_%H%M%S")


def decimate_to(x, fs, target):
    for q in (8, 8, 16):
        x = signal.decimate(x, q, ftype="fir", zero_phase=True)
        fs //= q
    return x, fs


def fetch_dec(key):
    local = os.path.join(SCRATCH, key.split("/")[-1])
    print(f"  downloading {key} ...", flush=True)
    urllib.request.urlretrieve(f"{BUCKET}/{key}", local)
    fs = sf.info(local).samplerate
    out = []
    with sf.SoundFile(local) as fh:
        while True:
            x = fh.read(fs * 10, dtype="float64")
            if len(x) < 400:
                break
            if x.ndim > 1:
                x = x[:, 0]
            xd, _ = decimate_to(x, fs, TARGET_FS)
            out.append(xd)
    os.remove(local)
    return np.concatenate(out), fstart(key)


def main():
    segs, t0 = [], None
    for k in FILES:
        try:
            xd, fst = fetch_dec(k)
        except Exception as e:
            print(f"  skip {k}: {e}", flush=True); continue
        t0 = fst if t0 is None else t0
        segs.append((fst, xd))
    if not segs:
        print("no data"); return
    fs = TARGET_FS
    total = int((segs[-1][0] - t0).total_seconds() * fs) + len(segs[-1][1])
    sig = np.zeros(total)
    for fst, xd in segs:
        i0 = int((fst - t0).total_seconds() * fs)
        sig[i0:i0 + len(xd)] = xd
    t_of = lambda i: t0 + datetime.timedelta(seconds=i / fs)

    # broadband T-phase envelope (2-90 Hz), 2 s smoothing
    sos = signal.butter(4, (2, 90), btype="bandpass", fs=fs, output="sos")
    env = np.abs(signal.hilbert(signal.sosfiltfilt(sos, sig)))
    sm = np.convolve(env, np.ones(2 * fs) / (2 * fs), "same")
    med = np.median(sm); mad = np.median(np.abs(sm - med)) + 1e-30
    z = 0.6745 * (sm - med) / mad
    ipk = int(np.argmax(sm))
    tpk = t_of(ipk)
    dt_pred = (tpk - PRED_ARR).total_seconds()
    print(f"\nM6.0 origin {ORIGIN:%H:%M:%S} UTC | predicted T-arrival "
          f"{PRED_ARR:%H:%M:%S} UTC (range {RANGE_KM:.0f} km @ {C_TPHASE} km/s)")
    print(f"strongest broadband transient: {tpk:%H:%M:%S} UTC  peak z={z[ipk]:.0f}")
    print(f"  offset from seismic prediction: {dt_pred:+.0f} s "
          f"(implied speed {RANGE_KM/((tpk-ORIGIN).total_seconds()):.3f} km/s)")
    verdict = ("MATCH — this is the T-phase boom" if abs(dt_pred) < 90
               else "does not match prediction")
    print(f"  -> {verdict}")

    # figure
    tm = (np.arange(total) / fs)
    tm_rel = tm - (PRED_ARR - t0).total_seconds()
    fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    f2, tt, Sxx = signal.spectrogram(sig, fs=fs, nperseg=1024, noverlap=768)
    ax[0].pcolormesh(tt - (PRED_ARR - t0).total_seconds(), f2,
                     10 * np.log10(Sxx + 1e-20), shading="auto", cmap="magma",
                     vmin=np.percentile(10*np.log10(Sxx+1e-20), 40),
                     vmax=np.percentile(10*np.log10(Sxx+1e-20), 99.5))
    ax[0].set_ylim(0, 125); ax[0].set_ylabel("Hz")
    ax[0].set_title(f"MARS 2026-01-16 — M6.0 off Oregon T-phase "
                    f"(predicted arrival = 0 s, {PRED_ARR:%H:%M:%S} UTC)")
    ax[1].plot(tm_rel, sm, lw=0.8, color="tab:blue", label="2-90 Hz envelope (2 s)")
    ax[1].plot(tm_rel[ipk], sm[ipk], "rv", label=f"peak (z={z[ipk]:.0f})")
    for a in ax:
        a.axvline(0, color="lime", ls="--", lw=2)
    ax[1].set_xlabel("seconds relative to predicted T-phase arrival (green dashed)")
    ax[1].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig("figures/fig_mars_tphase.png", dpi=115)
    print("wrote figures/fig_mars_tphase.png")


if __name__ == "__main__":
    main()
