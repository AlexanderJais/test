#!/usr/bin/env python3
"""Replicate the CI01 mechanism scan at CI04 and CI05.

Tests whether the "acoustic field goes quieter while the shrimp snap more"
divergence seen at CI01 around the 2019-07-16 event also appears on the two
other Channel Islands stations. Independent cables/sites: if the same
divergence recurs it is much harder to dismiss as local; if it does not, the
CI01 event stays single-station (as its snap-rate rise already was).

Same method and bands as src/mechanism_scan.py, but with a concatenated
reader so CI04's event window (which straddles a 6 h file boundary) can be
scanned continuously.
"""
import datetime
import re
import numpy as np
import soundfile as sf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mechanism_scan import (band_rms, snap_count, haversine_km, BANDS,
                            WIN_S, SPAN_MIN, MAX_LAG_MIN, C_KM_S)

EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)
EVENT_LATLON = (32.4894, -119.3647)

STATIONS = {
    "CI04": {"latlon": (33.849, -120.118),
             "files": ["data/omaha/SanctSound_CI04_02_671924265_190715235545.flac",
                       "data/omaha/SanctSound_CI04_02_671924265_190716055543.flac"]},
    "CI05": {"latlon": (34.0178, -119.3172),
             "files": ["data/omaha/SanctSound_CI05_02_671359013_190716014401.flac"]},
}


def fstart(p):
    m = re.search(r"_(\d{6})(\d{6})\.flac$", p)
    return datetime.datetime.strptime(m.group(1) + m.group(2), "%y%m%d%H%M%S")


class ConcatReader:
    """Read n frames starting at an absolute UTC time across contiguous files."""
    def __init__(self, files):
        self.meta = []           # (path, start_dt, frames, fs, handle)
        for p in files:
            info = sf.info(p)
            self.meta.append([p, fstart(p), info.frames, info.samplerate,
                              sf.SoundFile(p)])
        self.fs = self.meta[0][3]
        self.t0 = self.meta[0][1]
        end = self.meta[-1][1] + datetime.timedelta(
            seconds=self.meta[-1][2] / self.meta[-1][3])
        self.duration_s = (end - self.t0).total_seconds()

    def read_at(self, abs_dt, n):
        out = np.empty(0)
        need = n
        cur = abs_dt
        for p, sdt, frames, fs, fh in self.meta:
            edt = sdt + datetime.timedelta(seconds=frames / fs)
            if cur >= edt or need <= 0:
                continue
            if cur < sdt:                    # gap before this file -> stop
                break
            loc = int(round((cur - sdt).total_seconds() * fs))
            fh.seek(loc)
            x = fh.read(min(need, frames - loc), dtype="float64")
            if x.ndim > 1:
                x = x[:, 0]
            out = np.concatenate([out, x])
            need -= len(x)
            cur = sdt + datetime.timedelta(seconds=(loc + len(x)) / fs)
        return out


def scan(station):
    st = STATIONS[station]
    rdr = ConcatReader(st["files"])
    fs = rdr.fs
    n = int(WIN_S * fs)
    rng = haversine_km(st["latlon"], EVENT_LATLON)
    t0_s = (EVENT_UTC - rdr.t0).total_seconds() + rng / C_KM_S
    lo_s = max(0.0, t0_s - SPAN_MIN * 60)
    hi_s = min(rdr.duration_s, t0_s + SPAN_MIN * 60)

    starts_s = np.arange(lo_s, hi_s - WIN_S, WIN_S)
    tmin = starts_s / 60.0
    t0_min = t0_s / 60.0

    energy = {lbl: np.full(len(starts_s), np.nan) for lbl, _, _ in BANDS}
    snaprate = np.full(len(starts_s), np.nan)
    for i, s in enumerate(starts_s):
        x = rdr.read_at(rdr.t0 + datetime.timedelta(seconds=s), n)
        if len(x) < n * 0.9:
            continue
        for lbl, kind, p in BANDS:
            energy[lbl][i] = band_rms(x, fs, kind, p)
        snaprate[i] = snap_count(x, fs) * (60.0 / WIN_S)

    edb = {}
    for lbl in energy:
        e = energy[lbl]
        med = np.nanmedian(e)
        edb[lbl] = 20 * np.log10(e / med) if med > 0 else np.full_like(e, np.nan)

    # lead-lag
    lags = np.arange(-int(MAX_LAG_MIN * 60 / WIN_S),
                     int(MAX_LAG_MIN * 60 / WIN_S) + 1)
    sr = snaprate - np.nanmean(snaprate)
    leadlag = {}
    for lbl, _, _ in BANDS:
        if lbl.startswith("SNAP"):
            continue
        b = edb[lbl] - np.nanmean(edb[lbl])
        best_c, best_l = 0.0, 0
        for L in lags:
            if L >= 0:
                a, c = b[:len(b) - L], sr[L:]
            else:
                a, c = b[-L:], sr[:len(sr) + L]
            m = np.isfinite(a) & np.isfinite(c)
            if m.sum() < 10 or np.std(a[m]) == 0 or np.std(c[m]) == 0:
                continue
            cc = np.corrcoef(a[m], c[m])[0, 1]
            if abs(cc) > abs(best_c):
                best_c, best_l = cc, L
        leadlag[lbl] = (best_c, best_l * WIN_S / 60)

    before = (tmin >= t0_min - 20) & (tmin < t0_min)
    after = (tmin >= t0_min) & (tmin < t0_min + 20)
    step = {lbl: np.nanmean(edb[lbl][after]) - np.nanmean(edb[lbl][before])
            for lbl, _, _ in BANDS}

    snap_step = (np.nanmean(snaprate[after]) - np.nanmean(snaprate[before]))
    snap_pct = 100 * snap_step / np.nanmean(snaprate[before])

    print(f"\n=== {station} (range {rng:.0f} km, arrival t0={t0_min:.1f} min "
          f"since {rdr.t0:%H:%M} UTC, {int(np.isfinite(snaprate).sum())} bins) ===")
    print(f"snap rate step at event: {snap_step:+.0f}/min ({snap_pct:+.1f}%)")
    print("energy step at event (after-before, +/-20 min):")
    for lbl, _, _ in BANDS:
        c, l = leadlag.get(lbl, (float('nan'), float('nan')))
        tag = "" if lbl.startswith("SNAP") else f"   lead-lag |r|={c:+.2f}@{l:+.1f}min"
        print(f"  {lbl:18s}: {step[lbl]:+.2f} dB{tag}")
    return dict(tmin=tmin, t0_min=t0_min, edb=edb, snaprate=snaprate,
                snap_pct=snap_pct, step=step)


def main():
    res = {s: scan(s) for s in STATIONS}
    fig, axes = plt.subplots(1, 2, figsize=(17, 6))
    colors = {"infrasound 1-20Hz": "tab:blue", "low 20-100Hz": "tab:orange",
              "ship 100-500Hz": "tab:green", "EM 60Hz": "tab:red"}
    for ax, st in zip(axes, STATIONS):
        r = res[st]
        for lbl in colors:
            ax.plot(r["tmin"], r["edb"][lbl], color=colors[lbl], lw=0.7, label=lbl)
        ax.axvline(r["t0_min"], color="red", ls="--", lw=2)
        ax.set_xlabel("min since file start (UTC)")
        ax.set_ylabel("band energy (dB re median)")
        axr = ax.twinx()
        axr.plot(r["tmin"], r["snaprate"], color="black", lw=1.5, label="SNAP rate")
        axr.set_ylabel("snap rate/min")
        ax.set_title(f"{st}: snap step {r['snap_pct']:+.1f}%  "
                     f"(low {r['step']['low 20-100Hz']:+.1f} dB, "
                     f"ship {r['step']['ship 100-500Hz']:+.1f} dB)")
        if st == list(STATIONS)[0]:
            ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("Replication of the CI01 mechanism scan at CI04 & CI05 "
                 "(2019-07-16 Omaha event)")
    fig.tight_layout()
    fig.savefig("figures/fig_mechanism_replicate.png", dpi=110)
    print("\nwrote figures/fig_mechanism_replicate.png")


if __name__ == "__main__":
    main()
