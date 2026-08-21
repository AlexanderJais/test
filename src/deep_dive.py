#!/usr/bin/env python3
"""Full-coverage deep dive around a candidate signature.

The daily screen samples only the first 60 s of one 10-minute file per day
(~0.07% of the record). That is fine for a first pass but cannot tell a
brief transient from a persistent source, and judges candidates on 60 s.
This tool removes that limit for a chosen date range: it streams EVERY
10-minute file (144/day = continuous 24 h), splits each into 60 s windows
so memory stays bounded, mines discharge-like chains in every window, and
assembles a continuous timeline of the candidate signature.

What the timeline reveals that a snippet cannot:
* persistence  -- a fixed installation is present hour after hour at a
  stable rate/centroid/echo; a transient craft appears then leaves;
* motion       -- a transiting source drifts smoothly in received level
  and (Doppler) repetition rate; a fixed source does not;
* diel pattern -- biologic sources track dawn/dusk and foraging bouts.

Files are deleted immediately after processing (disk-bounded). Restartable:
windows already in the output JSONL are skipped.

Usage:
    python3 src/deep_dive.py --start 2021-12-27 --days 3 \
        [--cen-lo 12000 --cen-hi 22000] [--workers 3] [--every 1]
"""

import argparse
import datetime
import json
import os
import re
import struct
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fetch_data import BUCKET_FMT, http_get, parse_riff_header
from plasma_transient_search import (detect, load_pressure_pa, screen_events,
                                     DET_BAND)
from spectral_analysis import VOLTS_PEAK, V_TO_UPA

WIN_S = 60


def list_day_keys(year, month, day):
    base = BUCKET_FMT.format(year=year)
    prefix = f"{month:02d}/MARS_{year}{month:02d}{day:02d}"
    keys, token = [], None
    while True:
        url = f"{base}/?list-type=2&max-keys=1000&prefix={prefix}"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token)}"
        xml = http_get(url).decode()
        keys += re.findall(r"<Key>([^<]+)</Key>", xml)
        m = re.search(r"<NextContinuationToken>([^<]+)</", xml)
        if not m:
            break
        token = m.group(1)
    return sorted(keys)


def file_start_dt(key):
    m = re.search(r"MARS_(\d{8})_(\d{6})", key)
    return datetime.datetime.strptime(m.group(1) + m.group(2),
                                      "%Y%m%d%H%M%S")


def process_file(args):
    key, year, cen_lo, cen_hi = args
    base = BUCKET_FMT.format(year=year)
    url = f"{base}/{key}"
    t0 = file_start_dt(key)
    try:
        head = http_get(url, {"Range": "bytes=0-4095"})
        data_off, fs, ch, bits = parse_riff_header(head)
        bps = fs * ch * (bits // 8)
        blob = http_get(url, {"Range":
                              f"bytes=0-{data_off + 600 * bps - 1}"})
        raw = np.frombuffer(blob[data_off:], dtype=np.uint8)
        n = (len(raw) // 3)
        raw = raw[:n * 3].reshape(-1, 3).astype(np.int32)
        samp = (raw[:, 0] | (raw[:, 1] << 8) | (raw[:, 2] << 16))
        samp = np.where(samp & 0x800000, samp - 0x1000000, samp)
        p = samp.astype(np.float64) / 0x7FFFFF * VOLTS_PEAK * V_TO_UPA * 1e-6
        out = []
        for w in range(0, len(p) - WIN_S * fs + 1, WIN_S * fs):
            seg = p[w:w + WIN_S * fs]
            x, env, peaks, heights = detect(seg, fs)
            _, chains = screen_events(x, env, fs, peaks, heights,
                                      feature_cap=400)
            wt = t0 + datetime.timedelta(seconds=w / fs)
            for c in chains:
                if not c["class"].startswith("DISCHARGE"):
                    continue
                if cen_lo <= c["med_centroid_hz"] <= cen_hi:
                    out.append({"utc": wt.isoformat(), "n": c["n"],
                                "rate": c["rate_hz"], "cv": c["ipi_cv"],
                                "cen": c["med_centroid_hz"],
                                "echo": c["med_echo_lag_ms"],
                                "hf": c.get("med_hf_frac", 0.0),
                                "n_events": int(len(peaks))})
        return {"key": key, "utc": t0.isoformat(), "status": "ok",
                "hits": out}
    except Exception as e:
        return {"key": key, "utc": t0.isoformat(),
                "status": f"error: {type(e).__name__}: {e}", "hits": []}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--days", type=int, default=1)
    ap.add_argument("--cen-lo", type=float, default=0.0)
    ap.add_argument("--cen-hi", type=float, default=128000.0)
    ap.add_argument("--every", type=int, default=1,
                    help="process every Nth 10-min file (1 = all)")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--out", default="results/deep_dive.jsonl")
    args = ap.parse_args()

    d0 = datetime.date.fromisoformat(args.start)
    keys = []
    for i in range(args.days):
        d = d0 + datetime.timedelta(days=i)
        dk = list_day_keys(d.year, d.month, d.day)[::args.every]
        keys += [(k, d.year, args.cen_lo, args.cen_hi) for k in dk]
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)["key"] for l in open(args.out) if l.strip()}
    keys = [k for k in keys if k[0] not in done]
    print(f"{len(keys)} ten-minute files to process "
          f"({args.days} day(s) from {args.start}, every {args.every})")

    n_hits = 0
    with open(args.out, "a") as fh, \
            ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_file, k): k[0] for k in keys}
        for i, fu in enumerate(as_completed(futs), 1):
            rec = fu.result()
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n_hits += len(rec["hits"])
            if i % 20 == 0 or rec["hits"]:
                print(f"  {i}/{len(keys)} files, "
                      f"{n_hits} in-band discharge chains so far"
                      + (f"  <-- {rec['key']}: "
                         f"{len(rec['hits'])} hits" if rec["hits"] else ""),
                      flush=True)
    print(f"done: {n_hits} in-band discharge-chain detections")


if __name__ == "__main__":
    main()
