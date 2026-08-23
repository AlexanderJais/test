#!/usr/bin/env python3
"""Scan CI01 snap rate around the July 2019 SoCal UAP incident nights.

Incidents (Navy destroyer 'drone-swarm' series, Channel Islands warning
area; local evening -> next-day ~06:00 UTC):
  2019-07-15  (night of Jul 14: USS Kidd, first detections)
  2019-07-16  (night of Jul 15: Peralta/Russell/John Finn; USS Omaha sphere)
  2019-07-26  (night of Jul 25)
  2019-07-31  (night of Jul 30)

For each incident night, at station CI01 (18 m, Channel Islands), compute
the snapping-shrimp snap rate (impulses/min, 2-20 kHz) and the after/before
ratio (+/-20 min) at 06:00 UTC. Control null = same ratio slid across all
downloaded non-incident CI01 days. Output: per-incident snap-rate elevation
and its percentile against the control null.
"""
import glob, re, datetime, json
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BIN = 30; HALF_MIN = 20; CLICK_BAND = (2000.0, 20000.0)
INCIDENTS = {  # UTC date at ~06:00 : label
    "2019-07-15": "Jul14 night (Kidd)",
    "2019-07-16": "Jul15 night (Omaha sphere)",
    "2019-07-26": "Jul25 night",
    "2019-07-31": "Jul30 night",
}


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1)+m.group(2), "%y%m%d%H%M%S")


def snap_series(path):
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


def ratio_at(series, t_utc, half_bins):
    f0, rs = series
    c = int((t_utc - f0).total_seconds()/BIN)
    if half_bins <= c < len(rs)-half_bins:
        b = rs[c-half_bins:c].mean(); a = rs[c:c+half_bins].mean()
        if b > 0:
            return float(a/b), float(b*60/BIN), float(a*60/BIN)
    return None, None, None


def main():
    half_bins = HALF_MIN*(60//BIN)
    # incident files: event-day (no ctrl_/inc_ prefix) + inc_ prefixed
    inc_files = (sorted(glob.glob("data/omaha/SanctSound_CI01*.flac"))
                 + sorted(glob.glob("data/omaha/inc_*CI01*.flac")))
    ctrl_files = sorted(glob.glob("data/omaha/ctrl_*CI01*.flac"))
    inc_series = {}
    for p in inc_files:
        f0, rs = snap_series(p)
        inc_series[f0.strftime("%Y-%m-%d")] = (f0, rs)
    # control null
    null = []
    for p in ctrl_files:
        f0, rs = snap_series(p)
        for c in range(half_bins, len(rs)-half_bins):
            b = rs[c-half_bins:c].mean()
            if b > 0:
                null.append(rs[c:c+half_bins].mean()/b)
    null = np.array(null)
    p95 = float(np.percentile(null, 95)); p99 = float(np.percentile(null, 99))

    rows = []
    for datestr, label in INCIDENTS.items():
        if datestr not in inc_series:
            rows.append({"date": datestr, "label": label, "status": "no data"})
            continue
        t = datetime.datetime.fromisoformat(datestr+"T06:00:00")
        r, before, after = ratio_at(inc_series[datestr], t, half_bins)
        if r is None:
            rows.append({"date": datestr, "label": label, "status": "no window"})
            continue
        pval = float(np.mean(null >= r)) if len(null) else None
        rows.append({"date": datestr, "label": label, "ratio": round(r, 3),
                     "before_per_min": round(before, 0),
                     "after_per_min": round(after, 0),
                     "pct_change": round((r-1)*100, 1),
                     "p_vs_control": pval, "status": "ok"})
    json.dump({"null_p95": p95, "null_p99": p99, "n_null": len(null),
               "incidents": rows}, open("results/incident_scan.json", "w"),
              indent=1)
    print(f"CI01 snap-rate at UAP incident nights (control null p95={p95:.3f}, "
          f"p99={p99:.3f}, n={len(null)}):")
    print(f"{'date':12s} {'label':26s} {'before':>8s} {'after':>8s} "
          f"{'%chg':>6s} {'p':>7s}")
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['date']:12s} {r['label']:26s}  {r['status']}")
            continue
        print(f"{r['date']:12s} {r['label']:26s} {r['before_per_min']:8.0f} "
              f"{r['after_per_min']:8.0f} {r['pct_change']:+6.1f} "
              f"{r['p_vs_control']:7.3f}")

    # figure: per-incident click-rate curve around 06:00 UTC + elevation bars
    ok = [r for r in rows if r["status"] == "ok"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 5))
    for datestr, label in INCIDENTS.items():
        if datestr not in inc_series: continue
        f0, rs = inc_series[datestr]
        t = datetime.datetime.fromisoformat(datestr+"T06:00:00")
        c = int((t-f0).total_seconds()/BIN)
        w = 90*(60//BIN)
        seg = rs[max(0, c-w):c+w]
        tt = (np.arange(len(seg))-min(c, w))*BIN/60
        a1.plot(tt, seg/np.median(seg), lw=1,
                label=f"{datestr} {label}")
    a1.axvline(0, color="red", ls="--"); a1.set_xlabel("min from 06:00 UTC")
    a1.set_ylabel("CI01 snap rate / median"); a1.legend(fontsize=7)
    a1.set_title("CI01 snap-rate around each UAP incident night")
    labels = [r["label"] for r in ok]; vals = [r["pct_change"] for r in ok]
    ps = [r["p_vs_control"] for r in ok]
    colors = ["red" if p < 0.05 else "0.5" for p in ps]
    a2.bar(range(len(ok)), vals, color=colors)
    a2.axhline((p95-1)*100, color="k", ls=":", label="control 95th pct")
    a2.set_xticks(range(len(ok)))
    a2.set_xticklabels([f"{r['date'][5:]}\n{r['label'][:10]}" for r in ok],
                       fontsize=7)
    a2.set_ylabel("snap-rate change after/before (%)")
    a2.set_title("Elevation per incident (red = p<0.05 vs control null)")
    a2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/fig_incident_scan.png", dpi=140)
    print("wrote figures/fig_incident_scan.png")


if __name__ == "__main__":
    main()
