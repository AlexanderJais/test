#!/usr/bin/env python3
"""Decisive test of the CI01 click-rate anomaly: is the dip-then-rise at
~06:00 UTC unique to the event day, or a recurring diel/tidal cycle?

Overlays the CI01 click-rate profile vs LOCAL time-of-day for the event
day and 3 control days (no reported event), all covering 06:00 UTC. If
the inflection recurs on control days at the same clock time, it is the
shallow-station snapping-shrimp / diel cycle coinciding with t0 (anomaly
explained). If it is unique to the event day, the p=0.003 finding stands
as genuinely event-associated.
"""
import datetime, glob, re, json
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)
BIN = 30
CLICK_BAND = (2000.0, 20000.0)


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def clickrate_window(path, utc_lo, utc_hi):
    f0 = fstart(path); info = sf.info(path); fs = info.samplerate
    f1 = f0 + datetime.timedelta(seconds=info.frames/fs)
    a, b = max(utc_lo, f0), min(utc_hi, f1)
    if a >= b:
        return [], []
    sos = signal.butter(4, CLICK_BAND, btype="bandpass", fs=fs, output="sos")
    n = BIN*fs
    ts, rs = [], []
    with sf.SoundFile(path) as fh:
        i0 = int((a-f0).total_seconds()*fs)
        i1 = int((b-f0).total_seconds()*fs)
        for off in range(i0, i1-n, n):
            fh.seek(off); x = fh.read(n, dtype="float64")
            if x.ndim > 1: x = x[:, 0]
            xb = signal.sosfiltfilt(sos, x); env = np.abs(signal.hilbert(xb))
            med = np.median(env); mad = np.median(np.abs(env-med))+1e-30
            pk, _ = signal.find_peaks(env, height=med+8*mad, distance=int(0.01*fs))
            wt = f0 + datetime.timedelta(seconds=off/fs+BIN/2)
            ts.append(wt); rs.append(len(pk)*(60/BIN))
    return ts, rs


def main():
    # UTC hour-of-day relative to 06:00 for overlay
    fig, ax = plt.subplots(figsize=(13, 6))
    days = [("EVENT 2019-07-16", datetime.date(2019,7,16),
             sorted(glob.glob("data/omaha/SanctSound_CI01*.flac")), "red", 2.2),
            ("control 2019-07-09", datetime.date(2019,7,9),
             sorted(glob.glob("data/omaha/ctrl_*190709*.flac")), "0.4", 1.0),
            ("control 2019-07-23", datetime.date(2019,7,23),
             sorted(glob.glob("data/omaha/ctrl_*190723*.flac")), "tab:blue", 1.0),
            ("control 2019-07-02", datetime.date(2019,7,2),
             sorted(glob.glob("data/omaha/ctrl_*190702*.flac")), "tab:green", 1.0)]
    summary = {}
    for label, d, paths, color, lw in days:
        if not paths:
            print(f"{label}: no file"); continue
        lo = datetime.datetime(d.year,d.month,d.day,6,0,0)-datetime.timedelta(minutes=45)
        hi = datetime.datetime(d.year,d.month,d.day,6,0,0)+datetime.timedelta(minutes=45)
        ts, rs = [], []
        for p in paths:
            a, b = clickrate_window(p, lo, hi)
            ts += a; rs += b
        if not ts:
            print(f"{label}: no coverage of 06:00 UTC"); continue
        t0 = datetime.datetime(d.year,d.month,d.day,6,0,0)
        xrel = [(t-t0).total_seconds()/60 for t in ts]
        rs = np.array(rs, float)
        # normalize each day to its own median for shape comparison
        rsn = rs/np.median(rs)
        ax.plot(xrel, rsn, color=color, lw=lw, label=label, alpha=0.85)
        # before/after ratio around 06:00
        xrel = np.array(xrel)
        bef = rsn[(xrel>=-20)&(xrel<0)]; aft = rsn[(xrel>0)&(xrel<=20)]
        ratio = float(np.mean(aft)/np.mean(bef)) if len(bef) and len(aft) else float('nan')
        summary[label] = round(ratio, 3)
        print(f"{label}: after/before click ratio at 06:00 UTC = {ratio:.3f} "
              f"({'RISE' if ratio>1.05 else 'drop' if ratio<0.95 else 'flat'})")
    ax.axvline(0, color="red", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel("minutes relative to 06:00 UTC (event arrival time)")
    ax.set_ylabel("CI01 click rate / median (normalized)")
    ax.set_title("CI01 click-rate shape at 06:00 UTC: event day vs 3 control days\n"
                 "(if controls show the same dip->rise, it is a diel/tidal cycle, not the event)")
    ax.legend()
    fig.tight_layout(); fig.savefig("figures/fig_ci01_control_days.png", dpi=140)
    json.dump(summary, open("results/ci01_control_days.json","w"), indent=1)
    print("\nafter/before ratios:", summary)
    print("wrote figures/fig_ci01_control_days.png")


if __name__ == "__main__":
    main()
