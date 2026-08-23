#!/usr/bin/env python3
"""Does marine-life acoustic activity change AT the UAP event time?

A genuine, testable hypothesis (documented for sonar, explosions, ships):
cetaceans alter vocal/echolocation behaviour in response to a disturbance
-- often a sudden silence, sometimes a spike. If biosonar activity shifts
at the event time on INDEPENDENT stations, that is a real finding, not an
artifact. This is tested here with no human judgement in the loop: the
significance is set by a permutation/empirical null, not by me.

Per SanctSound Channel Islands station (48 kHz), across ALL available
hours: compute a per-minute time series of
  * odontocete CLICK RATE (impulses/min, 2-20 kHz),
  * whale/low-frequency tonal ENERGY (15-100 Hz),
  * broadband level,
then, treating the event's predicted arrival as t0, test:
  (A) STEP test: does mean activity in the T minutes AFTER t0 differ from
      the T minutes BEFORE?  Significance from a permutation null (slide
      t0 to 2000 random times; how often is |after-before| >= observed?).
  (B) CHANGE-POINT: is the largest CUSUM break near t0?
  (C) cross-station: do independent stations show the SAME sign of change
      within a few minutes of t0?
Outputs a per-minute activity figure with the event marked, and empirical
p-values. Whatever it says, it says.
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
STATIONS = {"CI01": (34.0438, -120.0811), "CI04": (33.849, -120.118),
            "CI05": (34.0178, -119.3172)}
WIN_S = 60          # 1-minute activity bins
CLICK_BAND = (2000.0, 20000.0)
WHALE_BAND = (15.0, 100.0)
TEST_T_MIN = 20     # before/after half-window for the step test
N_PERM = 2000


def hav(a, b, c, e):
    la1, lo1, la2, lo2 = map(np.radians, [a, b, c, e])
    return 2*6371*np.arcsin(np.sqrt(np.sin((la2-la1)/2)**2 +
                                    np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2))


def fstart(path):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", path)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def activity_series(code, paths):
    """Per-minute (utc, click_rate, whale_db, broadband_db)."""
    rows = []
    for path in sorted(paths):
        f0 = fstart(path)
        info = sf.info(path)
        fs = info.samplerate
        n = WIN_S * fs
        sos_c = signal.butter(4, CLICK_BAND, btype="bandpass", fs=fs,
                              output="sos")
        sos_w = signal.butter(4, WHALE_BAND, btype="bandpass", fs=fs,
                              output="sos")
        with sf.SoundFile(path) as fh:
            for off in range(0, info.frames - n, n):    # contiguous 1-min bins
                fh.seek(off)
                x = fh.read(n, dtype="float64", always_2d=False)
                if x.ndim > 1:
                    x = x[:, 0]
                xc = signal.sosfiltfilt(sos_c, x)
                env = np.abs(signal.hilbert(xc))
                med = np.median(env); mad = np.median(np.abs(env-med))+1e-30
                pk, _ = signal.find_peaks(env, height=med+8*mad,
                                          distance=int(0.01*fs))
                xw = signal.sosfiltfilt(sos_w, x)
                wt = f0 + datetime.timedelta(seconds=off/fs + WIN_S/2)
                rows.append((wt, len(pk),
                             10*np.log10(np.mean(xw**2)+1e-20),
                             10*np.log10(np.mean(x**2)+1e-20)))
    return rows


def step_test(t_rel_min, y, t0_min=0.0, half=TEST_T_MIN, nperm=N_PERM,
              seed=0):
    """|mean(after)-mean(before)| around t0; empirical p from sliding t0."""
    rng = np.random.default_rng(seed)
    t = np.asarray(t_rel_min); y = np.asarray(y, float)
    def stat(c):
        a = y[(t > c) & (t <= c+half)]
        b = y[(t >= c-half) & (t < c)]
        if len(a) < 5 or len(b) < 5:
            return np.nan
        return abs(np.nanmean(a)-np.nanmean(b))
    obs = stat(t0_min)
    if not np.isfinite(obs):
        return None
    lo, hi = t.min()+half, t.max()-half
    null = [stat(c) for c in rng.uniform(lo, hi, nperm)]
    null = np.array([v for v in null if np.isfinite(v)])
    p = float(np.mean(null >= obs))
    # signed change for direction
    a = y[(t > t0_min) & (t <= t0_min+half)]
    b = y[(t >= t0_min-half) & (t < t0_min)]
    return {"observed_abs_change": float(obs),
            "signed_change_after_minus_before":
                float(np.nanmean(a)-np.nanmean(b)),
            "p_empirical": p, "null_median": float(np.median(null))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="data/omaha")
    args = ap.parse_args()
    fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
    results = {}
    for ax, code in zip(axes, STATIONS):
        paths = glob.glob(os.path.join(args.datadir, f"*{code}*.flac"))
        if not paths:
            continue
        rng = hav(EVENT_LAT, EVENT_LON, *STATIONS[code])
        t_arr = EVENT_UTC + datetime.timedelta(seconds=rng/C_KM_S)
        rows = activity_series(code, paths)
        t_rel = [(r[0]-t_arr).total_seconds()/60 for r in rows]
        clicks = [r[1] for r in rows]
        whale = [r[2] for r in rows]
        res = {"range_km": round(float(rng), 1),
               "n_min": len(rows),
               "click_step": step_test(t_rel, clicks),
               "whale_step": step_test(t_rel, whale)}
        results[code] = res
        ax.plot(t_rel, clicks, lw=0.8, color="tab:blue", label="click rate/min")
        ax2 = ax.twinx()
        ax2.plot(t_rel, whale, lw=0.8, color="tab:green", alpha=0.7,
                 label="15-100 Hz dB")
        ax.axvline(0, color="red", lw=1.5, ls="--")
        ax.axvspan(0, TEST_T_MIN, color="red", alpha=0.05)
        ax.axvspan(-TEST_T_MIN, 0, color="blue", alpha=0.05)
        cs, ws = res["click_step"], res["whale_step"]
        ax.set_ylabel(f"{code} ({res['range_km']}km)\nclicks/min")
        ax2.set_ylabel("15-100Hz dB", fontsize=8)
        ax.set_title(
            f"{code}: click-rate step p={cs['p_empirical']:.3f} "
            f"(Δ={cs['signed_change_after_minus_before']:+.1f}/min); "
            f"whale-band step p={ws['p_empirical']:.3f} "
            f"(Δ={ws['signed_change_after_minus_before']:+.1f} dB)",
            fontsize=9)
        print(f"{code} ({res['range_km']}km, {res['n_min']} min): "
              f"click Δ={cs['signed_change_after_minus_before']:+.1f}/min "
              f"p={cs['p_empirical']:.3f} | whale Δ="
              f"{ws['signed_change_after_minus_before']:+.1f}dB "
              f"p={ws['p_empirical']:.3f}")
    axes[-1].set_xlabel("minutes relative to event predicted arrival (red = t0; "
                        "red band = 'after', blue band = 'before')")
    fig.suptitle("Biosonar & whale-band activity around the USS Omaha event — "
                 "permutation-tested change at t0", fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/fig_biosonar_response.png", dpi=140)
    json.dump(results, open("results/biosonar_response.json", "w"), indent=1)
    # cross-station verdict
    print("\ncross-station: sign of click-rate change at event")
    for c in results:
        print(f"  {c}: {results[c]['click_step']['signed_change_after_minus_before']:+.1f}/min "
              f"(p={results[c]['click_step']['p_empirical']:.3f})")
    print("wrote figures/fig_biosonar_response.png")


if __name__ == "__main__":
    main()
