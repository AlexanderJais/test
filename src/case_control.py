#!/usr/bin/env python3
"""Unbiased case-vs-control analysis at the closest mic to a UAP event.

Instead of matching a predefined signature, this asks the data an open
question: does the event window differ from ordinary recordings of the
SAME hydrophone, in ANY frequency band or acoustic feature?

Method (per SanctSound Channel Islands station, 48 kHz):
  * slice all available hours into 60 s windows every 5 min;
  * for each window compute an interpretable feature vector -- PSD in 16
    log-spaced bands (10 Hz-24 kHz), broadband RMS, spectral centroid,
    envelope kurtosis (impulsiveness), and an impulsive-click count;
  * CASE = windows within +/-3 min of the event's predicted acoustic
    arrival; CONTROL = all windows >30 min from the event (the mic's own
    natural variability -- diel, tide, shipping, weather, biology);
  * robust z-score (median/MAD over controls) of every feature for the
    case window -> which bands/features are anomalous;
  * multivariate anomaly score (count of |z|>3 and max|z|) for EVERY
    window, then report where the event window RANKS among all windows.
    If it is an ordinary window, that is a clean null; if it is a top
    outlier, report exactly which bands/clicks drive it.
Cross-station: an anomaly at the event time on INDEPENDENT stations is
meaningful; a one-station blip is local.
"""

import argparse
import datetime
import glob
import json
import os
import re

import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)
EVENT_LAT, EVENT_LON = 32.4894, -119.3647
C_KM_S = 1.487
STATIONS = {
    "CI01": (34.0438, -120.0811), "CI04": (33.849, -120.118),
    "CI05": (34.0178, -119.3172),
}
NBANDS = 16
WIN_S = 60
STEP_S = 300
CASE_TOL_S = 180
CTRL_EXCLUDE_S = 1800


def hav(a, b, c, e):
    la1, lo1, la2, lo2 = map(np.radians, [a, b, c, e])
    return 2*6371*np.arcsin(np.sqrt(np.sin((la2-la1)/2)**2 +
                                    np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2))


def fstart(path):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", path)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def features(x, fs, band_edges):
    f, pxx = signal.welch(x, fs=fs, nperseg=8192)
    feat = {}
    for i in range(len(band_edges)-1):
        m = (f >= band_edges[i]) & (f < band_edges[i+1])
        feat[f"b{i:02d}_{int(band_edges[i])}-{int(band_edges[i+1])}Hz"] = \
            float(10*np.log10(np.mean(pxx[m]) + 1e-20)) if m.any() else np.nan
    feat["broadband_db"] = float(10*np.log10(np.mean(pxx)+1e-20))
    feat["centroid_hz"] = float(np.sum(f*pxx)/np.sum(pxx))
    env = np.abs(signal.hilbert(signal.sosfiltfilt(
        signal.butter(4, (1000, min(20000, fs/2*0.95)), btype="bandpass",
                      fs=fs, output="sos"), x)))
    m, s = env.mean(), env.std()
    feat["kurtosis"] = float(((env-m)**4).mean()/(s**4+1e-30))
    med = np.median(env)
    mad = np.median(np.abs(env-med))+1e-30
    pk, _ = signal.find_peaks(env, height=med+8*mad, distance=int(0.01*fs))
    feat["click_count"] = int(len(pk))
    return feat


def scan_station(code, paths, band_edges):
    st = STATIONS[code]
    rng = hav(EVENT_LAT, EVENT_LON, st[0], st[1])
    t_arr = EVENT_UTC + datetime.timedelta(seconds=rng/C_KM_S)
    rows = []
    for path in sorted(paths):
        f0 = fstart(path)
        info = sf.info(path)
        fs = info.samplerate
        nfile = info.frames
        with sf.SoundFile(path) as fh:
            for off in range(0, int(nfile - WIN_S*fs), STEP_S*fs):
                fh.seek(off)
                x = fh.read(WIN_S*fs, dtype="float64", always_2d=False)
                if x.ndim > 1:
                    x = x[:, 0]
                wt = f0 + datetime.timedelta(seconds=off/fs + WIN_S/2)
                feat = features(x, fs, band_edges)
                feat["utc"] = wt.isoformat()
                feat["dt_event_s"] = (wt - t_arr).total_seconds()
                rows.append(feat)
    return rows, t_arr, rng


def analyze(code, rows, band_edges):
    feat_keys = [k for k in rows[0] if k not in ("utc", "dt_event_s")]
    case = [r for r in rows if abs(r["dt_event_s"]) <= CASE_TOL_S]
    ctrl = [r for r in rows if abs(r["dt_event_s"]) > CTRL_EXCLUDE_S]
    if not case or len(ctrl) < 10:
        return None
    # robust stats over controls
    med = {k: np.median([r[k] for r in ctrl]) for k in feat_keys}
    mad = {k: np.median(np.abs(np.array([r[k] for r in ctrl])-med[k]))+1e-9
           for k in feat_keys}
    def z(r, k):
        return (r[k]-med[k])/(1.4826*mad[k])
    # event = the single case window closest to predicted arrival
    ev = min(case, key=lambda r: abs(r["dt_event_s"]))
    zev = {k: float(z(ev, k)) for k in feat_keys}
    # multivariate anomaly score for ALL windows: count |z|>3, and max|z|
    def score(r):
        zs = [abs(z(r, k)) for k in feat_keys]
        return sum(v > 3 for v in zs), max(zs)
    scored = [(r["utc"], float(r["dt_event_s"]), int(score(r)[0]), float(score(r)[1])) for r in rows]
    scored_sorted = sorted(scored, key=lambda s: (-s[2], -s[3]))
    ev_rank = next(i for i, s in enumerate(scored_sorted)
                   if s[0] == ev["utc"])
    return {"station": code, "n_windows": int(len(rows)), "n_control": int(len(ctrl)),
            "event_utc": ev["utc"], "event_dt_s": float(ev["dt_event_s"]),
            "z": {k: float(v) for k, v in zev.items()},
            "n_bands_over3": int(score(ev)[0]), "max_abs_z": float(score(ev)[1]),
            "event_rank": int(ev_rank),
            "top5_outlier_windows": [
                {"utc": s[0], "dt_event_s": round(s[1], 0),
                 "n_over3": s[2], "max_z": round(s[3], 1)}
                for s in scored_sorted[:5]]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="data/omaha")
    args = ap.parse_args()
    band_edges = np.logspace(np.log10(10), np.log10(24000), NBANDS+1)
    results = {}
    allrows = {}
    for code in STATIONS:
        paths = glob.glob(os.path.join(args.datadir, f"*{code}*.flac"))
        if not paths:
            continue
        rows, t_arr, rng = scan_station(code, paths, band_edges)
        allrows[code] = rows
        res = analyze(code, rows, band_edges)
        if res:
            res["range_km"] = round(float(rng), 1)
            results[code] = res
            print(f"{code} ({rng:.0f} km): {res['n_windows']} windows, "
                  f"event window ranks #{res['event_rank']+1}/"
                  f"{res['n_windows']} by anomaly; "
                  f"{res['n_bands_over3']} features |z|>3, "
                  f"max|z|={res['max_abs_z']:.1f}")
    json.dump(results, open("results/case_control.json", "w"), indent=1)

    # figure: per-band z of the event window at each station
    feat_keys = [k for k in allrows[list(allrows)[0]][0]
                 if k.startswith("b")]
    fig, ax = plt.subplots(figsize=(14, 6))
    xlab = [k.split("_", 1)[1] for k in feat_keys]
    x = np.arange(len(feat_keys))
    for code in results:
        zv = [results[code]["z"][k] for k in feat_keys]
        ax.plot(x, zv, "o-", label=f"{code} ({results[code]['range_km']} km)")
    for extra in ["broadband_db", "centroid_hz", "kurtosis", "click_count"]:
        pass
    ax.axhspan(-3, 3, color="0.9", zorder=0, label="±3σ control band")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(xlab, rotation=60, ha="right",
                                         fontsize=7)
    ax.set_ylabel("robust z-score of EVENT window vs control distribution")
    ax.set_title("USS Omaha event window vs same-station controls — per-band "
                 "anomaly (closest mic CI04 = 167 km). Inside grey = normal.")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/fig_case_control.png", dpi=150)
    print("wrote figures/fig_case_control.png and results/case_control.json")
    # also print the non-band features for the event
    for code in results:
        r = results[code]
        print(f"\n{code} event-window z: broadband={r['z']['broadband_db']:+.1f} "
              f"centroid={r['z']['centroid_hz']:+.1f} "
              f"kurtosis={r['z']['kurtosis']:+.1f} "
              f"clicks={r['z']['click_count']:+.1f}")
        print(f"  top outlier windows (any time): "
              + "; ".join(f"{o['utc'][11:19]}Z dt={o['dt_event_s']:+.0f}s "
                          f"n>3σ={o['n_over3']}" for o in r["top5_outlier_windows"]))


if __name__ == "__main__":
    main()
