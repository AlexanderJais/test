#!/usr/bin/env python3
"""Empirical null for the CI01 click-rate rise.

Slides the +/-20 min before/after window across all 9 CI01 control days
(no reported event) to build the natural distribution of the after/before
click-rate ratio at this station, then reports where the EVENT day's ratio
(measured at the predicted arrival) falls -> a real empirical p-value for
"is the event-coincident rise unusual for this hydrophone".
"""
import glob, re, datetime, json
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BIN = 30; HALF = 20
CLICK_BAND = (2000.0, 20000.0)


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def clickrate_series(path):
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
    return np.array(rs, float)


def ratios_over_day(rs, bins_per_min=60//BIN):
    """after/before ratio for every valid center bin in the day."""
    half = HALF*bins_per_min
    out = []
    for c in range(half, len(rs)-half):
        b = rs[c-half:c].mean(); a = rs[c:c+half].mean()
        if b > 0:
            out.append(a/b)
    return out


def main():
    # event ratio at 06:00 UTC on CI01 event file
    evp = sorted(glob.glob("data/omaha/SanctSound_CI01*.flac"))[0]
    ev_rs = clickrate_series(evp)
    f0 = fstart(evp); fs = sf.info(evp).samplerate
    ev_center_s = (datetime.datetime(2019,7,16,6,0,0)-f0).total_seconds()
    c = int(ev_center_s/BIN)
    half = HALF*(60//BIN)
    ev_ratio = ev_rs[c:c+half].mean()/ev_rs[c-half:c].mean()

    # null from all control days
    null = []
    ctrls = sorted(glob.glob("data/omaha/ctrl_*CI01*.flac"))
    for p in ctrls:
        null += ratios_over_day(clickrate_series(p))
    null = np.array(null)
    p_emp = float(np.mean(null >= ev_ratio))
    print(f"control days: {len(ctrls)}; null samples (sliding windows): {len(null)}")
    print(f"EVENT after/before ratio at 06:00 UTC = {ev_ratio:.3f}")
    print(f"null ratio pct[50,90,95,99]={np.percentile(null,[50,90,95,99]).round(3)}")
    print(f"empirical p (fraction of natural windows >= event) = {p_emp:.4f}")

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.hist(null, bins=80, color="0.7", label=f"natural CI01 windows (n={len(null)}, {len(ctrls)} control days)")
    ax.axvline(ev_ratio, color="red", lw=2,
               label=f"EVENT day at 06:00 UTC = {ev_ratio:.3f} (p={p_emp:.4f})")
    ax.axvline(np.percentile(null,95), color="k", ls=":", label="95th pct of null")
    ax.set_xlabel("after/before click-rate ratio (±20 min)")
    ax.set_ylabel("count")
    ax.set_title("CI01: is the event-coincident click-rate rise unusual?\n"
                 "empirical null from sliding the window across control days")
    ax.legend()
    fig.tight_layout(); fig.savefig("figures/fig_ci01_null.png", dpi=140)
    json.dump({"event_ratio": ev_ratio, "p_empirical": p_emp,
               "n_control_days": len(ctrls), "n_null": len(null),
               "null_p95": float(np.percentile(null,95)),
               "null_p99": float(np.percentile(null,99))},
              open("results/ci01_null.json","w"), indent=1)
    print("wrote figures/fig_ci01_null.png")


if __name__ == "__main__":
    main()
