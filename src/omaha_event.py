#!/usr/bin/env python3
"""Event-directed forensics: USS Omaha 'splash' incident.

Documented event: 2019-07-15 ~23:00 PDT = 2019-07-16 ~06:00 UTC at
32.4894 N, 119.3647 W (Pentagon-confirmed video; spherical ~2 m object
descends to the sea surface; 'splash' called; submarine search found
nothing). Three SanctSound Channel Islands stations were recording:

  CI01 34.0438 -120.0811  18 m   range 185 km  travel ~124 s
  CI04 33.849  -120.118  153 m   range 166 km  travel ~112 s
  CI05 34.0178 -119.3172 136 m   range 170 km  travel ~114 s

Approach (transparent forensic scan, not black-box classification):
for each station, read the window around the event (+travel time,
+/- generous uncertainty for report time and SoundTrap clock drift),
scan for impulsive transients in a low-frequency 'entry thump' band
(10-500 Hz, where a water-entry event carries energy that survives
~170 km) and broadband, produce spectrograms and a ranked transient
table, then cross-match stations within a +/-3 min drift-tolerant
window. Any coincident LF transient gets waveform-level examination
for slam/pinch doublet structure. A null yields a source-level upper
bound at the event position.

Usage: python3 src/omaha_event.py [--window-min 90]
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
    "CI01": {"lat": 34.0438, "lon": -120.0811, "depth": 18},
    "CI04": {"lat": 33.849, "lon": -120.118, "depth": 153},
    "CI05": {"lat": 34.0178, "lon": -119.3172, "depth": 136},
}
THUMP_BAND = (10.0, 500.0)
BB_BAND = (300.0, 22000.0)
MAD_K = 7.0


def hav_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*6371*np.arcsin(np.sqrt(a))


def fname_start(path):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", path)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def bp(x, fs, band):
    ny = fs / 2
    hi = min(band[1], ny * 0.95)
    sos = signal.butter(4, (band[0], hi), btype="bandpass", fs=fs,
                        output="sos")
    return signal.sosfiltfilt(sos, x)


def read_window_decimated(paths, w0, w1, target_fs=2000):
    """Memory-safe read: stream file chunks, decimate to target_fs for the
    LF analysis; also return per-chunk broadband transient events."""
    lf_parts, bb_events, fs_native = [], [], None
    t_abs0 = None
    for path in sorted(paths):
        f0 = fname_start(path)
        info = sf.info(path)
        fs = info.samplerate
        fs_native = fs
        f1 = f0 + datetime.timedelta(seconds=info.frames / fs)
        a, b = max(w0, f0), min(w1, f1)
        if a >= b:
            continue
        if t_abs0 is None:
            t_abs0 = a
        dec = max(1, int(round(fs / target_fs)))
        i0 = int((a - f0).total_seconds() * fs)
        n_total = int((b - a).total_seconds() * fs)
        chunk = 10 * 60 * fs                      # 10-min chunks
        with sf.SoundFile(path) as fh:
            fh.seek(i0)
            done = 0
            while done < n_total:
                x = fh.read(min(chunk, n_total - done), dtype="float64",
                            always_2d=False)
                if len(x) == 0:
                    break
                if x.ndim > 1:
                    x = x[:, 0]
                # LF: anti-alias + decimate
                lf_parts.append(signal.decimate(x, dec, ftype="fir",
                                                zero_phase=True))
                # broadband transient scan within the chunk
                xb = bp(x, fs, BB_BAND)
                env = np.abs(signal.hilbert(xb))
                med = np.median(env)
                mad = np.median(np.abs(env - med)) + 1e-30
                pk, props = signal.find_peaks(env, height=med + MAD_K * mad,
                                              distance=int(0.5 * fs))
                tc0 = a + datetime.timedelta(seconds=done / fs)
                for i, h in sorted(zip(pk, props["peak_heights"]),
                                   key=lambda z: -z[1])[:4]:
                    bb_events.append({
                        "band": "broadband",
                        "utc": (tc0 + datetime.timedelta(
                            seconds=float(i / fs))).isoformat(),
                        "snr": round(float(h / (med + MAD_K * mad)), 2)})
                done += len(x)
                del x, xb, env
    if not lf_parts:
        return None
    return t_abs0, np.concatenate(lf_parts), target_fs, bb_events, fs_native


def scan_station(code, paths, window_min):
    st = STATIONS[code]
    rng = hav_km(EVENT_LAT, EVENT_LON, st["lat"], st["lon"])
    t_arr = EVENT_UTC + datetime.timedelta(seconds=rng / C_KM_S)
    w0 = t_arr - datetime.timedelta(minutes=window_min)
    w1 = t_arr + datetime.timedelta(minutes=window_min)
    out = {"station": code, "range_km": round(float(rng), 1),
           "predicted_arrival_utc": t_arr.isoformat(),
           "window": [w0.isoformat(), w1.isoformat()],
           "events": [], "files": []}
    rw = read_window_decimated(paths, w0, w1)
    if rw is None:
        out["status"] = "no coverage"
        return out, None
    t_abs0, xlf, fsd, bb_events, fs_native = rw
    out["files"] = [os.path.basename(p) for p in paths]
    out["fs"] = fs_native
    out["coverage_s"] = round(len(xlf) / fsd, 1)

    # LF thump scan on the decimated stream
    xb = bp(xlf, fsd, THUMP_BAND)
    env = np.abs(signal.hilbert(xb))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    thr = med + MAD_K * mad
    pk, props = signal.find_peaks(env, height=thr, distance=int(0.5 * fsd))
    out["thump_rms"] = float(np.sqrt(np.mean(xb ** 2)))
    out["thump_thr"] = float(thr)
    for i, h in sorted(zip(pk, props["peak_heights"]),
                       key=lambda z: -z[1])[:12]:
        tt = t_abs0 + datetime.timedelta(seconds=float(i / fsd))
        a0, b0 = max(0, i - fsd), min(len(env), i + fsd)
        e = env[a0:b0]
        dur = float(np.sum(e > 0.1 * h)) / fsd
        out["events"].append({
            "band": "thump", "utc": tt.isoformat(),
            "snr": round(float(h / thr), 2), "dur_s": round(dur, 3),
            "dt_from_pred_s": round((tt - t_arr).total_seconds(), 1)})
    # broadband events from chunked scan (annotate offset from prediction)
    for e in bb_events:
        e["dt_from_pred_s"] = round((datetime.datetime.fromisoformat(
            e["utc"]) - t_arr).total_seconds(), 1)
    out["events"] += sorted(bb_events, key=lambda e: -e["snr"])[:12]
    out["status"] = "ok"
    return out, (t_abs0, xlf, fsd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window-min", type=int, default=90)
    ap.add_argument("--datadir", default="data/omaha")
    args = ap.parse_args()

    results, waves = [], {}
    for code in STATIONS:
        paths = sorted(glob.glob(os.path.join(
            args.datadir, f"*{code}*.flac")))
        if not paths:
            results.append({"station": code, "status": "no files"})
            continue
        r, w = scan_station(code, paths, args.window_min)
        results.append(r)
        if w:
            waves[code] = w
        print(f"{code}: {r['status']} range={r.get('range_km')}km "
              f"cover={r.get('coverage_s')}s "
              f"events={len(r.get('events', []))}")

    # cross-station coincidence: LF thump events within +/-180 s
    coinc = []
    evs = [(r["station"], e) for r in results for e in r.get("events", [])
           if e["band"] == "thump"]
    for i in range(len(evs)):
        for j in range(i + 1, len(evs)):
            si, ei = evs[i]
            sj, ej = evs[j]
            if si == sj:
                continue
            dt = abs((datetime.datetime.fromisoformat(ei["utc"])
                      - datetime.datetime.fromisoformat(ej["utc"]))
                     .total_seconds())
            if dt <= 180:
                coinc.append({"a": si, "b": sj, "utc_a": ei["utc"],
                              "utc_b": ej["utc"], "dt_s": round(dt, 1),
                              "snr_a": ei["snr"], "snr_b": ej["snr"]})
    print(f"cross-station thump coincidences (+/-180 s): {len(coinc)}")

    os.makedirs("results", exist_ok=True)
    json.dump({"event_utc": EVENT_UTC.isoformat(),
               "event_pos": [EVENT_LAT, EVENT_LON],
               "stations": results, "coincidences": coinc},
              open("results/omaha_event.json", "w"), indent=1)

    # spectrogram figure: LF band per station around predicted arrival
    fig, axes = plt.subplots(len(waves), 1, figsize=(14, 3.2 * len(waves)),
                             sharex=False)
    for ax, (code, (t0, x, fs)) in zip(np.atleast_1d(axes), waves.items()):
        f, t, sxx = signal.spectrogram(x, fs=fs, nperseg=4096,
                                       noverlap=2048)
        db = 10 * np.log10(np.maximum(sxx, 1e-30))
        m = ax.pcolormesh(t / 60, f, db, cmap="magma", shading="auto",
                          vmin=np.percentile(db, 30),
                          vmax=np.percentile(db, 99.5), rasterized=True)
        st = STATIONS[code]
        rngkm = hav_km(EVENT_LAT, EVENT_LON, st["lat"], st["lon"])
        t_pred = (EVENT_UTC + datetime.timedelta(seconds=rngkm / C_KM_S)
                  - t0).total_seconds() / 60
        ax.axvline(t_pred, color="cyan", ls="--", lw=1.2)
        ax.set_yscale("log")
        ax.set_ylim(8, 990)
        ax.set_ylabel(f"{code}\nHz")
        ax.set_title(f"{code} — window start {t0.isoformat()}Z; cyan = "
                     f"predicted arrival for reported event time", fontsize=9)
    np.atleast_1d(axes)[-1].set_xlabel("minutes into window")
    fig.suptitle("USS Omaha event window — SanctSound Channel Islands, "
                 "8 Hz–2 kHz", fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/fig_omaha_event.png", dpi=140)
    print("wrote figures/fig_omaha_event.png and results/omaha_event.json")


if __name__ == "__main__":
    main()
