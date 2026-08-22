#!/usr/bin/env python3
"""M1 pilot: full-coverage multi-detector hunt on MARS.

Streams every 10-minute MARS file over a date range (robust per-60s-window
range fetch, stream-and-delete, restartable) and runs the calibrated M0
detector suite on each 60 s window:

  * entry-doublet (D2)  -> slam+pinch pairs, body size from t_p
  * quench-tail  (D4)   -> hot-body flag, GATED on each doublet (coincidence)
  * Doppler tracker (D1/D3) in two bands: a clean anomaly band 55-95 kHz and
    the calibrated 28-44 kHz band -> moving narrowband sources + speed

Because M0 measured the entry-doublet as storm-blind and high-FP, ranking
is by COINCIDENCE strength, not any single limb:

  rank 1  HOT ENTRY   : doublet AND quench tail at the same time
  rank 2  KINEMATIC   : Doppler transit with a physical speed (clean, 0 FP
                        in calibration) -- the real 'bogies' (ships/movers)
  rank 3  DOUBLET+???  : doublet alone (mostly false; logged for coincidence
                        with the OOI LF limb + external catalogs in M2/M3)

One JSONL record per file. Usage:
  python3 src/m1_pilot.py --start 2022-07-15 --days 3 [--workers 2]
  python3 src/m1_pilot.py --aggregate
"""

import argparse
import datetime
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from deep_dive import list_day_keys, file_start_dt, robust_range
from fetch_data import BUCKET_FMT, http_get, parse_riff_header
from spectral_analysis import VOLTS_PEAK, V_TO_UPA
import entry_doublet as ed
import quench_tail as qt
import doppler_track as dp

DOPPLER_BANDS = [(55000.0, 95000.0), (28000.0, 44000.0)]
WIN_S = 60


def _decode(blob):
    raw = np.frombuffer(blob, dtype=np.uint8)
    n = len(raw) // 3
    raw = raw[:n * 3].reshape(-1, 3).astype(np.int32)
    s = raw[:, 0] | (raw[:, 1] << 8) | (raw[:, 2] << 16)
    s = np.where(s & 0x800000, s - 0x1000000, s)
    return s.astype(np.float64) / 0x7FFFFF * VOLTS_PEAK * V_TO_UPA * 1e-6


def process_file(args):
    key, year = args
    base = BUCKET_FMT.format(year=year)
    url = f"{base}/{key}"
    t0 = file_start_dt(key)
    rec = {"key": key, "utc": t0.isoformat(),
           "hot_entry": [], "kinematic": [], "doublet_only": [], "n_win": 0}
    try:
        head = http_get(url, {"Range": "bytes=0-4095"})
        data_off, fs, ch, bits = parse_riff_header(head)
        wb = WIN_S * fs * ch * (bits // 8)
        for wi in range(10):
            a = data_off + wi * wb
            try:
                p = _decode(robust_range(url, a, a + wb - 1))
            except Exception:
                continue
            rec["n_win"] += 1
            wt = t0 + datetime.timedelta(seconds=wi * WIN_S)
            # entry doublets + quench coincidence
            for d in ed.detect(p, fs):
                tail = qt.detect_tail(p, fs, d["t_slam_s"],
                                      max_tail_s=min(40.0,
                                                     WIN_S - d["t_slam_s"]))
                item = {"utc": wt.isoformat(), **d}
                if tail and tail["is_decaying"]:
                    item["quench"] = tail
                    rec["hot_entry"].append(item)
                else:
                    rec["doublet_only"].append(item)
            # Doppler transits
            for band in DOPPLER_BANDS:
                fit = dp.detect(p, fs, band)
                if fit and fit["is_transit"]:
                    rec["kinematic"].append({"utc": wt.isoformat(),
                                             "band": band, **fit})
        rec["status"] = "ok" if rec["n_win"] >= 8 else f"partial:{rec['n_win']}"
    except Exception as e:
        rec["status"] = f"error: {type(e).__name__}: {e}"
    return rec


def run(start, days, workers, out):
    d0 = datetime.date.fromisoformat(start)
    keys = []
    for i in range(days):
        d = d0 + datetime.timedelta(days=i)
        keys += [(k, d.year) for k in list_day_keys(d.year, d.month, d.day)]
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["key"] for l in open(out) if l.strip()}
    keys = [k for k in keys if k[0] not in done]
    print(f"{len(keys)} files to screen ({days} d from {start})")
    tot = {"hot_entry": 0, "kinematic": 0, "doublet_only": 0}
    with open(out, "a") as fh, ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(process_file, k): k[0] for k in keys}
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            fh.write(json.dumps(r) + "\n"); fh.flush()
            for k in tot:
                tot[k] += len(r.get(k, []))
            if r.get("hot_entry") or r.get("kinematic"):
                print(f"  *** {r['key']}: {len(r['hot_entry'])} HOT-ENTRY, "
                      f"{len(r['kinematic'])} KINEMATIC", flush=True)
            if i % 25 == 0:
                print(f"  {i}/{len(keys)} files | hot={tot['hot_entry']} "
                      f"kin={tot['kinematic']} doublet-only="
                      f"{tot['doublet_only']}", flush=True)
    print(f"done: {tot}")


def aggregate(out_glob="results/m1_pilot*.jsonl"):
    recs = []
    for p in sorted(glob.glob(out_glob)):
        recs += [json.loads(l) for l in open(p) if l.strip()]
    ok = [r for r in recs if r.get("status", "").startswith(("ok", "partial"))]
    hot = [(r["utc"], h) for r in ok for h in r.get("hot_entry", [])]
    kin = [(r["utc"], k) for r in ok for k in r.get("kinematic", [])]
    ndo = sum(len(r.get("doublet_only", [])) for r in ok)
    mins = len(ok) * 10
    print(f"files screened: {len(ok)} (~{mins} min coverage)")
    print(f"RANK 1 hot-entry coincidences: {len(hot)}")
    for u, h in hot[:20]:
        print(f"   {u}  t_p={h['tp_s']}s a~{h['implied_body_a_m']}m "
              f"quench {h['quench']['tail_duration_s']}s "
              f"decay {h['quench']['decay_db_per_s']}dB/s")
    print(f"RANK 2 kinematic (Doppler) transits: {len(kin)}")
    for u, k in sorted(kin, key=lambda x: -x[1]["swing_frac"])[:20]:
        print(f"   {u}  v={k['v_kn']}kn R_cpa={k['r_cpa_m']:.0f}m "
              f"f0={k['f0']:.0f}Hz swing={k['swing_frac']} band={k['band']}")
    print(f"RANK 3 doublet-only (high-FP, for M2/M3 coincidence): {ndo} "
          f"= {ndo/max(mins,1):.2f}/min")
    json.dump({"files": len(ok), "coverage_min": mins,
               "hot_entry": len(hot), "kinematic": len(kin),
               "doublet_only": ndo}, open("results/m1_summary.json", "w"),
              indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start")
    ap.add_argument("--days", type=int, default=3)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default="results/m1_pilot.jsonl")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()
    if args.aggregate:
        aggregate()
    else:
        run(args.start, args.days, args.workers, args.out)


if __name__ == "__main__":
    main()
