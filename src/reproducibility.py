#!/usr/bin/env python3
"""Reproducibility / placebo test of the CI01 event-coincident snap-rate rise.

A single anomaly found after many looks is weak. Two things test whether it
is more than chance, both using the identical after/before snap-rate ratio
(+/-20 min) at each station:

  (1) PLACEBO NULL per station: slide the test window to every valid center
      across all available data for that station -> distribution of ratios a
      random 'event time' would produce. The event's empirical p is where
      the TRUE Omaha arrival falls on that null.
  (2) CROSS-STATION REPLICATION: a real response to one physical event should
      appear at the independent stations too. Compare the event ratio at
      CI01 / CI04 / CI05.

Verdict logic (stated, not fudged): the rise replicates -> compelling; it is
a lone tail outlier at one station -> consistent with look-elsewhere chance.
"""
import glob, re, datetime, json
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)
EVENT_LAT, EVENT_LON = 32.4894, -119.3647
C_KM_S = 1.487
STA = {"CI01": (34.0438, -120.0811), "CI04": (33.849, -120.118),
       "CI05": (34.0178, -119.3172)}
BIN = 30; HALF_MIN = 20; CLICK_BAND = (2000.0, 20000.0)


def hav(a, b, c, e):
    a, b, c, e = map(np.radians, [a, b, c, e])
    return 2*6371*np.arcsin(np.sqrt(np.sin((c-a)/2)**2+np.cos(a)*np.cos(c)*np.sin((e-b)/2)**2))


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def snap_series(path):
    """(start_dt, np.array of snaps per BIN s)."""
    info = sf.info(path); fs = info.samplerate; n = BIN*fs
    sos = signal.butter(4, CLICK_BAND, btype="bandpass", fs=fs, output="sos")
    rs = []
    with sf.SoundFile(path) as fh:
        for off in range(0, info.frames-n, n):
            fh.seek(off); x = fh.read(n, dtype="float64")
            if x.ndim > 1: x = x[:, 0]
            xb = signal.sosfiltfilt(sos, x); env = np.abs(signal.hilbert(xb))
            med = np.median(env); mad = np.median(np.abs(env-med))+1e-30
            pk, _ = signal.find_peaks(env, height=med+8*mad, distance=int(0.01*fs))
            rs.append(len(pk))
    return fstart(path), np.array(rs, float)


def ratio_null(series_list, half_bins):
    out = []
    for _, rs in series_list:
        for c in range(half_bins, len(rs)-half_bins):
            b = rs[c-half_bins:c].mean()
            if b > 0:
                out.append(rs[c:c+half_bins].mean()/b)
    return np.array(out)


def event_ratio(series_list, t_arr, half_bins):
    """after/before ratio centered on the bin containing t_arr."""
    for f0, rs in series_list:
        end = f0 + datetime.timedelta(seconds=len(rs)*BIN)
        if f0 <= t_arr < end:
            c = int((t_arr-f0).total_seconds()/BIN)
            if half_bins <= c < len(rs)-half_bins:
                return rs[c:c+half_bins].mean()/rs[c-half_bins:c].mean()
    return None


def main():
    half_bins = HALF_MIN*(60//BIN)
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    res = {}
    for ax, code in zip(axes, STA):
        # event-day file(s) = those NOT prefixed ctrl_
        ev_paths = sorted(p for p in glob.glob(f"data/omaha/*{code}*.flac")
                          if "/ctrl_" not in p and not p.split("/")[-1].startswith("ctrl_"))
        ctrl_paths = sorted(glob.glob(f"data/omaha/ctrl_*{code}*.flac"))
        all_paths = ev_paths + ctrl_paths
        series = [snap_series(p) for p in all_paths]
        rng = hav(EVENT_LAT, EVENT_LON, *STA[code])
        t_arr = EVENT_UTC + datetime.timedelta(seconds=rng/C_KM_S)
        er = event_ratio(series, t_arr, half_bins)
        null = ratio_null(series, half_bins)
        p = float(np.mean(null >= er)) if er else None
        res[code] = {"range_km": round(float(rng), 1), "event_ratio": er,
                     "p_empirical": p, "n_null": len(null),
                     "n_files": len(all_paths),
                     "null_p95": float(np.percentile(null, 95)),
                     "null_p99": float(np.percentile(null, 99))}
        ax.hist(null, bins=80, color="0.7",
                label=f"placebo null (n={len(null)}, {len(all_paths)} files)")
        if er:
            ax.axvline(er, color="red", lw=2,
                       label=f"Omaha event = {er:.3f} (p={p:.3f})")
        ax.axvline(np.percentile(null, 95), color="k", ls=":", label="95th pct")
        ax.set_title(f"{code} ({res[code]['range_km']} km): event ratio "
                     f"{er:.3f}, p={p:.3f}" if er else f"{code}: no event coverage",
                     fontsize=10)
        ax.set_ylabel("count"); ax.legend(fontsize=8)
        print(f"{code} ({res[code]['range_km']}km): event after/before={er:.3f}, "
              f"p={p:.3f}, null95={res[code]['null_p95']:.3f}")
    axes[-1].set_xlabel("after/before snap-rate ratio (±20 min)")
    fig.suptitle("Reproducibility: does the CI01 rise replicate across "
                 "independent stations for the SAME event?", fontsize=12)
    fig.tight_layout(); fig.savefig("figures/fig_reproducibility.png", dpi=140)
    json.dump(res, open("results/reproducibility.json", "w"), indent=1)
    # verdict
    rises = [c for c in res if res[c]["event_ratio"] and
             res[c]["p_empirical"] is not None and res[c]["p_empirical"] < 0.05]
    print(f"\nstations with a significant (p<0.05) coincident rise: {rises}")
    print("replicates across independent stations -> compelling; lone outlier "
          "-> consistent with chance")
    print("wrote figures/fig_reproducibility.png")


if __name__ == "__main__":
    main()
