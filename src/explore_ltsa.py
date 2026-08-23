#!/usr/bin/env python3
"""Open-ended anomaly sweep of the Channel Islands hydrophone data.

Two modes:
  compute FILE  -> stream the whole file, build a Long-Term Spectral Average
                   (LTSA) plus per-column features; cache to FILE.ltsa.npz
  report        -> load all event-window caches, draw a multi-panel LTSA and
                   print an objective anomaly catalog + cross-station
                   coincidence check.

Per 30 s column we store:
  * PSD reduced to 160 log-spaced frequency bins (dB) -> the LTSA,
  * band RMS (dB) in infrasound/low/ship/mid/snap bands,
  * crest factor (max|x| / RMS) and kurtosis -> impulsiveness,
  * a tonal score: how peaky the spectrum is vs its local smooth (persistent
    narrow lines = machinery/ships).
No detector is tuned to a hypothesis here; this just surfaces what stands out.
"""
import sys
import glob
import datetime
import re
import numpy as np
import soundfile as sf
from scipy import signal, stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL_S = 30.0
NPERSEG = 4096
NF_LOG = 160
BANDS = [("infra", 1, 20), ("low", 20, 100), ("ship", 100, 500),
         ("mid", 500, 2000), ("snap", 2000, 20000)]
EVENT_FILES = [
    "data/omaha/SanctSound_CI01_02_671883305_190716033643.flac",
    "data/omaha/SanctSound_CI04_02_671924265_190715235545.flac",
    "data/omaha/SanctSound_CI04_02_671924265_190716055543.flac",
    "data/omaha/SanctSound_CI05_02_671359013_190716014401.flac",
]


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1) + m.group(2), "%y%m%d%H%M%S")


def compute(path):
    cache = path + ".ltsa.npz"
    info = sf.info(path)
    fs = info.samplerate
    n = int(COL_S * fs)
    ncol = info.frames // n
    f = np.fft.rfftfreq(NPERSEG, 1 / fs)
    # log-spaced frequency edges 10 Hz .. Nyquist
    fedges = np.logspace(np.log10(10), np.log10(fs / 2 * 0.99), NF_LOG + 1)
    fidx = np.digitize(f, fedges) - 1
    fcenters = np.sqrt(fedges[:-1] * fedges[1:])

    ltsa = np.full((NF_LOG, ncol), np.nan)
    band_db = {b[0]: np.full(ncol, np.nan) for b in BANDS}
    crest = np.full(ncol, np.nan)
    kurt = np.full(ncol, np.nan)
    tonal = np.full(ncol, np.nan)
    nuniq = np.full(ncol, np.nan)       # quantization richness (data quality)
    hp = signal.butter(2, 5.0, btype="high", fs=fs, output="sos")

    with sf.SoundFile(path) as fh:
        for c in range(ncol):
            fh.seek(c * n)
            x = fh.read(n, dtype="float64")
            if x.ndim > 1:
                x = x[:, 0]
            ff, pxx = signal.welch(x, fs=fs, nperseg=NPERSEG)   # detrends per seg
            pxx_db = 10 * np.log10(pxx + 1e-20)
            for k in range(NF_LOG):
                sel = fidx == k
                if sel.any():
                    ltsa[k, c] = pxx_db[sel].mean()
            df = ff[1] - ff[0]
            for name, lo, hi in BANDS:
                sel = (ff >= lo) & (ff < hi)
                band_db[name][c] = 10 * np.log10(pxx[sel].sum() * df + 1e-20)
            nuniq[c] = len(np.unique(x[::10]))
            xh = signal.sosfiltfilt(hp, x)                       # kill DC + swell
            rms = np.sqrt(np.mean(xh ** 2)) + 1e-20
            crest[c] = np.max(np.abs(xh)) / rms
            kurt[c] = stats.kurtosis(xh)
            sm = signal.medfilt(pxx_db, 51)
            tonal[c] = np.percentile(pxx_db - sm, 99)

    np.savez(cache, f0=fstart(path).timestamp(), fcenters=fcenters,
             ltsa=ltsa, crest=crest, kurt=kurt, tonal=tonal, nuniq=nuniq,
             **{f"band_{k}": v for k, v in band_db.items()})
    print(f"cached {cache}: {ncol} cols x {NF_LOG} freqs")


def station_of(p):
    m = re.search(r"_(CI\d\d)_", p)
    return m.group(1) if m else p


def report():
    caches = [p + ".ltsa.npz" for p in EVENT_FILES]
    Z = {p: np.load(p + ".ltsa.npz", allow_pickle=True) for p in EVENT_FILES}

    # ---- figure: LTSA per file ----
    fig, axes = plt.subplots(len(EVENT_FILES), 1, figsize=(15, 12), sharex=False)
    for ax, p in zip(axes, EVENT_FILES):
        z = Z[p]
        ltsa = z["ltsa"]; fc = z["fcenters"]
        t0 = datetime.datetime.fromtimestamp(z["f0"].item())
        tt = np.arange(ltsa.shape[1]) * COL_S / 3600.0
        vmin, vmax = np.nanpercentile(ltsa, [5, 99])
        im = ax.pcolormesh(tt, fc, ltsa, shading="auto", vmin=vmin, vmax=vmax,
                           cmap="magma")
        ax.set_yscale("log"); ax.set_ylim(10, fc.max())
        ax.set_ylabel(f"{station_of(p)}\nHz")
        ax.set_title(f"{station_of(p)}  start {t0:%m-%d %H:%M} UTC", fontsize=9,
                     loc="left")
        fig.colorbar(im, ax=ax, pad=0.01, label="dB")
    axes[-1].set_xlabel("hours since file start")
    fig.suptitle("Channel Islands LTSA (event-window files)")
    fig.tight_layout()
    fig.savefig("figures/fig_ltsa.png", dpi=110)
    print("wrote figures/fig_ltsa.png\n")

    # ---- objective anomaly catalog ----
    print("=== per-file anomaly flags (robust z >= 5 vs that file's own median) ===")
    events = []   # (utc, station, kind, value)
    for p in EVENT_FILES:
        z = Z[p]
        t0 = z["f0"].item()
        st = station_of(p)
        feats = {"crest": z["crest"], "kurt": z["kurt"], "tonal": z["tonal"],
                 "ship_dB": z["band_ship"], "low_dB": z["band_low"]}
        for name, v in feats.items():
            med = np.nanmedian(v)
            mad = np.nanmedian(np.abs(v - med)) + 1e-9
            zc = 0.6745 * (v - med) / mad
            for c in np.where(zc >= 5)[0]:
                utc = datetime.datetime.fromtimestamp(t0) + datetime.timedelta(seconds=c * COL_S)
                events.append((utc, st, name, float(v[c]), float(zc[c])))
    events.sort()
    # print top 25 by z
    top = sorted(events, key=lambda e: -e[4])[:25]
    for utc, st, name, val, zc in top:
        print(f"  {utc:%m-%d %H:%M:%S} {st:5s} {name:8s} val={val:8.2f}  z={zc:5.1f}")

    # ---- cross-station coincidences (same feature within +/-90 s on >=2 st) ----
    print("\n=== cross-station coincident flags (<=90 s apart, different stations) ===")
    found = 0
    for i in range(len(events)):
        u0, s0, n0, v0, z0 = events[i]
        for j in range(i + 1, len(events)):
            u1, s1, n1, v1, z1 = events[j]
            dt = abs((u1 - u0).total_seconds())
            if dt > 90:
                break
            if s0 != s1 and n0 == n1:
                print(f"  {n0:8s}: {s0} {u0:%H:%M:%S} <-> {s1} {u1:%H:%M:%S} "
                      f"(dt={dt:.0f}s, z={z0:.0f}/{z1:.0f})")
                found += 1
    if not found:
        print("  none -- flagged transients are local to single stations")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "compute":
        compute(sys.argv[2])
    elif len(sys.argv) >= 2 and sys.argv[1] == "report":
        report()
    else:
        print("usage: explore_ltsa.py compute FILE | report")
