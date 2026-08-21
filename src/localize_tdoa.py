#!/usr/bin/env python3
"""TDOA localization of an impulsive underwater source with the OOI
low-frequency hydrophone array — validated against a ground-truth event.

Method (classic passive-sonar / hydroacoustic monitoring):

1. band-pass each station to the T-phase band (3-40 Hz),
2. form smoothed Hilbert envelopes (T-phases are emergent wave packets;
   envelope correlation is the robust way to time them),
3. cross-correlate all station pairs -> time-differences-of-arrival,
4. grid-search latitude/longitude (and water-borne group speed c) for the
   source position whose predicted TDOAs best fit the measurements
   (hyperbolic fixing).

Validation target: USGS-reviewed M5.0 earthquake, 2023-01-11 10:17:18 UTC
at 43.9646 N, 128.7389 W. Its T-phase was recorded by 4 stations of the
array at 220-295 km range. The script reports the miss distance between
the acoustic fix and the USGS epicenter.

Usage:
    python3 src/localize_tdoa.py [--data data/array/event.npz]
"""

import argparse
import itertools
import json
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TRUTH = {"lat": 43.9646, "lon": -128.7389, "name": "USGS M5.0 us7000j3ld"}
ORIGIN_OFFSET_S = 60.0        # traces start 60 s before origin time
BAND = (3.0, 40.0)
ENV_LP_HZ = 0.4               # envelope smoothing
ENV_FS = 10.0                 # envelope resample rate
XCORR_WIN = (100.0, 400.0)    # s into trace: T-phase, excludes seismic P
MAX_LAG_S = 120.0
R_EARTH = 6371.0


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_EARTH * np.arcsin(np.sqrt(a))


def load_array(path):
    z = np.load(path)
    stations = sorted({k.split("_")[0] for k in z.files if k.endswith("_x")})
    out = {}
    for s in stations:
        out[s] = {"x": z[f"{s}_x"].astype(float), "fs": float(z[f"{s}_fs"]),
                  "lat": float(z[f"{s}_lat"]), "lon": float(z[f"{s}_lon"])}
    return out


def envelope(x, fs):
    sos = signal.butter(4, BAND, btype="bandpass", fs=fs, output="sos")
    xb = signal.sosfiltfilt(sos, x)
    env = np.abs(signal.hilbert(xb))
    sos_lp = signal.butter(4, ENV_LP_HZ, btype="low", fs=fs, output="sos")
    env = signal.sosfiltfilt(sos_lp, env)
    step = int(round(fs / ENV_FS))
    return env[::step]


def measure_tdoas(arr):
    """Envelope cross-correlation for every station pair.

    Returns list of (sta_i, sta_j, tdoa_seconds, peak_correlation) where
    tdoa > 0 means the signal reached sta_i EARLIER than sta_j.
    """
    i0, i1 = int(XCORR_WIN[0] * ENV_FS), int(XCORR_WIN[1] * ENV_FS)
    envs = {}
    for s, d in arr.items():
        e = envelope(d["x"], d["fs"])[i0:i1]
        envs[s] = (e - e.mean()) / e.std()
    out, max_lag = [], int(MAX_LAG_S * ENV_FS)
    for si, sj in itertools.combinations(sorted(envs), 2):
        a, b = envs[si], envs[sj]
        cc = signal.correlate(b, a, mode="full") / len(a)
        lags = signal.correlation_lags(len(b), len(a), mode="full")
        keep = np.abs(lags) <= max_lag
        cc, lags = cc[keep], lags[keep]
        k = int(np.argmax(cc))
        # parabolic sub-sample refinement
        if 0 < k < len(cc) - 1:
            denom = cc[k - 1] - 2 * cc[k] + cc[k + 1]
            k_frac = 0.5 * (cc[k - 1] - cc[k + 1]) / denom if denom else 0.0
        else:
            k_frac = 0.0
        tdoa = (lags[k] + k_frac) / ENV_FS
        out.append((si, sj, float(tdoa), float(cc[k])))
    return out


def grid_search(arr, tdoas, lat_rng=(42.0, 47.5), lon_rng=(-132.5, -124.5),
                step=0.02, c_range=(1.44, 1.52), c_step=0.005):
    lats = np.arange(*lat_rng, step)
    lons = np.arange(*lon_rng, step)
    glon, glat = np.meshgrid(lons, lats)
    dist = {s: haversine_km(glat, glon, d["lat"], d["lon"])
            for s, d in arr.items()}
    best = None
    for c in np.arange(c_range[0], c_range[1] + 1e-9, c_step):
        misfit = np.zeros_like(glat)
        wsum = 0.0
        for si, sj, tdoa, w in tdoas:
            pred = (dist[sj] - dist[si]) / c   # sta_i earlier => positive
            misfit += w * (pred - tdoa) ** 2
            wsum += w
        rms = np.sqrt(misfit / wsum)
        k = np.unravel_index(np.argmin(rms), rms.shape)
        if best is None or rms[k] < best["rms"]:
            best = {"lat": float(glat[k]), "lon": float(glon[k]),
                    "rms": float(rms[k]), "c": float(c),
                    "surface": rms, "lats": lats, "lons": lons}
    return best


def peak_moveout_check(arr, best_c):
    """Secondary validation: T-phase envelope PEAK times must move out as
    range/c across the array (up to one common offset for source duration).
    Peaks are used rather than onsets because seafloor hydrophones also
    record the earlier seismic (crustal) arrival, which contaminates
    threshold-based onset picks, while the energy peak is T-phase."""
    rows = []
    for s, d in arr.items():
        e = envelope(d["x"], d["fs"])
        i0, i1 = int(XCORR_WIN[0] * ENV_FS), int(XCORR_WIN[1] * ENV_FS)
        t_pk = XCORR_WIN[0] + np.argmax(e[i0:i1]) / ENV_FS - ORIGIN_OFFSET_S
        rng = haversine_km(TRUTH["lat"], TRUTH["lon"], d["lat"], d["lon"])
        rows.append({"station": s, "range_km": round(float(rng), 1),
                     "peak_s_after_origin": round(float(t_pk), 1),
                     "predicted_travel_s": round(float(rng / best_c), 1)})
    off = float(np.median([r["peak_s_after_origin"]
                           - r["predicted_travel_s"] for r in rows]))
    for r in rows:
        r["moveout_residual_s"] = round(
            r["peak_s_after_origin"] - r["predicted_travel_s"] - off, 1)
    return rows, off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/array/event.npz")
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()

    arr = load_array(args.data)
    print(f"stations: {', '.join(arr)}")
    tdoas = measure_tdoas(arr)
    for si, sj, tdoa, w in tdoas:
        print(f"  TDOA {si}->{sj}: {tdoa:+8.2f} s   (peak corr {w:.2f})")

    best = grid_search(arr, tdoas)
    miss = float(haversine_km(best["lat"], best["lon"],
                              TRUTH["lat"], TRUTH["lon"]))
    print(f"\nacoustic fix: {best['lat']:.3f} N, {best['lon']:.3f} E  "
          f"(c={best['c']:.3f} km/s, rms={best['rms']:.2f} s)")
    print(f"USGS truth  : {TRUTH['lat']:.3f} N, {TRUTH['lon']:.3f} E")
    print(f"miss distance: {miss:.1f} km")

    onsets, off = peak_moveout_check(arr, best["c"])
    print(f"envelope-peak moveout check (common offset {off:+.1f} s):")
    for r in onsets:
        print(f"  {r['station']}: peak {r['peak_s_after_origin']} s, "
              f"travel-time {r['predicted_travel_s']} s at "
              f"{r['range_km']} km -> residual {r['moveout_residual_s']} s")

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "localization.json"), "w") as fh:
        json.dump({"tdoas": [{"pair": f"{a}-{b}", "tdoa_s": t, "corr": w}
                             for a, b, t, w in tdoas],
                   "fix": {k: best[k] for k in ("lat", "lon", "rms", "c")},
                   "truth": TRUTH, "miss_km": round(miss, 1),
                   "peak_moveout": onsets, "peak_offset_s": off},
                  fh, indent=1)

    # --- record section -------------------------------------------------
    fig, axes = plt.subplots(len(arr), 1, figsize=(11, 7), sharex=True)
    for ax, (s, d) in zip(np.atleast_1d(axes), sorted(arr.items())):
        sos = signal.butter(4, BAND, btype="bandpass", fs=d["fs"],
                            output="sos")
        xb = signal.sosfiltfilt(sos, d["x"])
        t = np.arange(len(xb)) / d["fs"] - ORIGIN_OFFSET_S
        rng = haversine_km(TRUTH["lat"], TRUTH["lon"], d["lat"], d["lon"])
        ax.plot(t, xb, lw=0.3, color="k")
        ax.axvline(rng / best["c"], color="tab:red", ls="--", lw=1)
        ax.set_ylabel(f"{s}\n{rng:.0f} km", fontsize=8)
        ax.set_xlim(-60, 500)
    np.atleast_1d(axes)[-1].set_xlabel("time after earthquake origin [s]")
    fig.suptitle(f"T-phase record section, {TRUTH['name']} — red dashes: "
                 f"predicted water-borne arrival (range / {best['c']:.2f} "
                 "km/s)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_tphase.png"), dpi=150)

    # --- misfit map -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 7))
    m = ax.pcolormesh(best["lons"], best["lats"], best["surface"],
                      cmap="viridis", shading="auto",
                      vmax=np.percentile(best["surface"], 60))
    fig.colorbar(m, ax=ax, label="weighted TDOA rms misfit [s]")
    for s, d in arr.items():
        ax.plot(d["lon"], d["lat"], "^", color="w", mec="k", ms=9)
        ax.annotate(s, (d["lon"], d["lat"]), xytext=(4, 4),
                    textcoords="offset points", fontsize=8, color="w")
    ax.plot(TRUTH["lon"], TRUTH["lat"], "*", color="tab:red", ms=17,
            mec="k", label="USGS epicenter (truth)")
    ax.plot(best["lon"], best["lat"], "o", color="tab:orange", ms=9,
            mec="k", label=f"acoustic TDOA fix (miss {miss:.0f} km)")
    ax.set_xlabel("longitude [°E]")
    ax.set_ylabel("latitude [°N]")
    ax.set_title("Hyperbolic localization with the public OOI hydrophone "
                 "array", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_localization.png"), dpi=150)
    print(f"wrote {args.figdir}/fig_tphase.png and fig_localization.png")


if __name__ == "__main__":
    main()
