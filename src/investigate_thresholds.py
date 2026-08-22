#!/usr/bin/env python3
"""Threshold investigation: are the M1 clean-gate thresholds too strict?

The pilot logs only passes, so this re-runs the Doppler and quench detectors
on a sample of files and records RAW scores regardless of pass/fail, then
shows where the real population sits relative to each gate. If signals pile
up just under a threshold, it is too strict; if they sit far away, the 0/0
is genuine.
"""
import argparse, datetime, json, os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from deep_dive import list_day_keys, file_start_dt, robust_range
from fetch_data import BUCKET_FMT, http_get, parse_riff_header
from spectral_analysis import VOLTS_PEAK, V_TO_UPA
import entry_doublet as ed
import quench_tail as qt
import doppler_track as dp

BANDS = [(55000.0, 95000.0), (28000.0, 44000.0)]


def _decode(blob):
    raw = np.frombuffer(blob, dtype=np.uint8); n = len(raw)//3
    raw = raw[:n*3].reshape(-1,3).astype(np.int32)
    s = raw[:,0]|(raw[:,1]<<8)|(raw[:,2]<<16)
    s = np.where(s & 0x800000, s-0x1000000, s)
    return s.astype(np.float64)/0x7FFFFF*VOLTS_PEAK*V_TO_UPA*1e-6


def process(args):
    key, year = args
    url = f"{BUCKET_FMT.format(year=year)}/{key}"
    t0 = file_start_dt(key)
    dop, qu = [], []
    try:
        head = http_get(url, {"Range": "bytes=0-4095"})
        off, fs, ch, bits = parse_riff_header(head); wb = 60*fs*ch*(bits//8)
        for wi in range(10):
            a = off + wi*wb
            try: p = _decode(robust_range(url, a, a+wb-1))
            except Exception: continue
            for band in BANDS:
                t, f, snr, _ton = dp.track_ridge(p, fs, band)
                fit = dp.fit_scurve(t, f, snr) if len(t) else None
                if fit:
                    swing = (f.max()-f.min())/fit["f0"]
                    dop.append({"band": band[0], "v_kn": round(fit["v_ms"]/dp.KN,1),
                                "swing_frac": float(swing),
                                "resid_frac": fit["resid_hz"]/fit["f0"],
                                "n_ridge": len(t)})
            for d in ed.detect(p, fs):
                tail = qt.detect_tail(p, fs, d["t_slam_s"], max_tail_s=min(40.0,60-d["t_slam_s"]))
                if tail:
                    qu.append({"tp": d["tp_s"], "tail_windows": tail["tail_windows"],
                               "decay": tail["decay_db_per_s"],
                               "mean_level_up": tail["mean_level_up_db"],
                               "isi_cv": tail["median_isi_cv"]})
    except Exception as e:
        return {"key": key, "err": str(e), "dop": [], "qu": []}
    return {"key": key, "utc": t0.isoformat(), "dop": dop, "qu": qu}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=12)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    keys = []
    for d in [datetime.date(2022,7,15), datetime.date(2022,7,16), datetime.date(2022,7,17)]:
        keys += [(k, d.year) for k in list_day_keys(d.year, d.month, d.day)][::a.every]
    print(f"sampling {len(keys)} files")
    dop, qu = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for r in ex.map(process, keys):
            dop += r.get("dop", []); qu += r.get("qu", [])
    json.dump({"dop": dop, "qu": qu}, open("results/threshold_probe.json","w"))
    # Doppler distributions vs gates (swing>=3e-3, resid<=2e-3)
    sw = np.array([d["swing_frac"] for d in dop]); rf = np.array([d["resid_frac"] for d in dop])
    print(f"\nDOPPLER fits computed: {len(dop)}")
    print(f"  swing_frac pct[50,90,99,max]={np.percentile(sw,[50,90,99,100]).round(5)} (gate>=0.003)")
    print(f"  resid_frac pct[10,50,90]={np.percentile(rf,[10,50,90]).round(5)} (gate<=0.002)")
    print(f"  pass swing gate: {np.mean(sw>=0.003):.3f}; pass resid gate: {np.mean(rf<=0.002):.3f}; "
          f"pass BOTH: {np.mean((sw>=0.003)&(rf<=0.002)):.3f}")
    print("  top-8 by swing (closest to / over gate):")
    for d in sorted(dop,key=lambda x:-x["swing_frac"])[:8]:
        ok = d["swing_frac"]>=0.003 and d["resid_frac"]<=0.002
        print(f"    v={d['v_kn']}kn swing={d['swing_frac']:.4f} resid={d['resid_frac']:.4f} "
              f"n={d['n_ridge']} band={d['band']/1e3:.0f}k {'** PASSES' if ok else ''}")
    # quench decay distribution vs gate (decay < -0.02)
    dc = np.array([q["decay"] for q in qu]); tw = np.array([q["tail_windows"] for q in qu])
    print(f"\nQUENCH tails computed (on doublets): {len(qu)}")
    if len(qu):
        print(f"  decay_db_per_s pct[1,10,50]={np.percentile(dc,[1,10,50]).round(3)} (gate < -0.02 to flag hot)")
        print(f"  tail_windows pct[50,90,max]={np.percentile(tw,[50,90,100]).round(0)}")
        print(f"  would flag as decaying+sustained: {np.mean((dc<-0.02)&(tw>=4)):.3f}")
        print("  top-6 by (most negative decay, sustained):")
        for q in sorted([q for q in qu if q['tail_windows']>=4],key=lambda x:x['decay'])[:6]:
            print(f"    decay={q['decay']}dB/s windows={q['tail_windows']} level_up={q['mean_level_up']}dB isi_cv={q['isi_cv']}")


if __name__ == "__main__":
    main()
