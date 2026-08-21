#!/usr/bin/env python3
"""Archive-scale batch screening for plasma-discharge transient trains.

Streams a stratified sample of the MBARI MARS 256 kHz public archive —
60 s from EVERY DAY of a year, with the sampled hour rotating through
all 24 diel hours (hour = 7*day_of_year mod 24, 7 being coprime with 24
so the full diel cycle is covered every ~3.4 weeks) — and runs the
discharge-template transient screen from plasma_transient_search.py on
each slice. Slices are deleted after processing, so disk use stays
bounded no matter how many days are screened.

Each day yields one JSONL record: transient count, train statistics and
classifications, and any DISCHARGE-LIKE template matches (metronomic +
broadband + impulsive). Restartable: days already in the output file are
skipped.

Usage:
    python3 src/batch_screen.py --year 2022 [--workers 3]
    python3 src/batch_screen.py --aggregate   # summary + figure from all JSONL
"""

import argparse
import datetime
import glob
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_data import fetch_slice, list_keys
from plasma_transient_search import (detect, load_pressure_pa,
                                     screen_events)

FEATURE_CAP = 3000     # max per-event waveform-feature extractions per slice


def day_list(year):
    d0 = datetime.date(year, 1, 1)
    days = []
    d = d0
    while d.year == year:
        doy = (d - d0).days + 1
        days.append((d, (7 * doy) % 24))
        d += datetime.timedelta(days=1)
    return days


def screen_day(args):
    date, hour, seconds, tmpdir = args
    ymd = date.strftime("%Y%m%d")
    rec = {"date": str(date), "hour": hour}
    path = None
    try:
        os.makedirs(tmpdir, exist_ok=True)
        for h in (hour, (hour + 1) % 24):    # tolerate single-hour gaps
            prefix = f"{ymd}_{h:02d}"
            path = fetch_slice(str(date.year), f"{date.month:02d}", prefix,
                               seconds, tmpdir)
            if path:
                rec["hour"] = h
                break
        if not path:
            # bigger outage: fall back to the day's nearest available hour
            keys = list_keys(str(date.year), f"{date.month:02d}", f"{ymd}_")
            if keys:
                hours = sorted({int(k.split("_")[-1][:2]) for k in keys},
                               key=lambda h: min(abs(h - hour),
                                                 24 - abs(h - hour)))
                path = fetch_slice(str(date.year), f"{date.month:02d}",
                                   f"{ymd}_{hours[0]:02d}", seconds, tmpdir)
                rec["hour"] = hours[0]
        if not path:
            rec["status"] = "no_data"
            return rec
        p, fs = load_pressure_pa(path)
        x, env, peaks, heights = detect(p, fs)
        gap, chains = screen_events(x, env, fs, peaks, heights,
                                    feature_cap=FEATURE_CAP)
        rec["status"] = "ok"
        rec["n_events"] = int(len(peaks))
        rec["trains"] = gap
        rec["chains"] = chains
        rec["n_candidates"] = sum(
            1 for t in gap + chains if t["class"].startswith("DISCHARGE"))
    except Exception as e:
        rec["status"] = f"error: {type(e).__name__}: {e}"
    finally:
        if path and os.path.exists(path):
            os.remove(path)
    return rec


def run_year(year, seconds, workers, out_path, tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            done = {json.loads(l)["date"] for l in fh if l.strip()}
    todo = [(d, h) for d, h in day_list(year) if str(d) not in done]
    print(f"{year}: {len(todo)} days to screen ({len(done)} already done)")
    t0 = time.time()
    n_done, n_cand = 0, 0
    with open(out_path, "a") as out, \
            ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(screen_day, (d, h, seconds, tmpdir)): d
                for d, h in todo}
        for fu in as_completed(futs):
            rec = fu.result()
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n_done += 1
            n_cand += rec.get("n_candidates", 0)
            if rec.get("n_candidates"):
                print(f"  *** CANDIDATE on {rec['date']}: {rec['trains']}")
            if n_done % 20 == 0:
                rate = n_done / (time.time() - t0)
                eta = (len(todo) - n_done) / rate / 60
                print(f"  {n_done}/{len(todo)} days "
                      f"({rate*60:.1f}/min, ETA {eta:.0f} min), "
                      f"candidates so far: {n_cand}", flush=True)
    print(f"{year} complete: {n_done} days, {n_cand} candidates")


TAG = ""


def aggregate(figdir="figures", out="results"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    recs = []
    pat = f"batch_screen_*{TAG}.jsonl" if TAG else "batch_screen_2???.jsonl"
    for path in sorted(glob.glob(os.path.join(out, pat))):
        with open(path) as fh:
            recs += [json.loads(l) for l in fh if l.strip()]
    ok = [r for r in recs if r.get("status") == "ok"]
    n_cand = sum(r.get("n_candidates", 0) for r in ok)
    cls_days = {}
    for r in ok:
        for t in r["trains"] + r.get("chains", []):
            key = t["class"].split(" (")[0]
            cls_days.setdefault(key, set()).add(r["date"])
    summary = {
        "days_screened": len(ok),
        "days_no_data": sum(1 for r in recs if r.get("status") == "no_data"),
        "days_error": sum(1 for r in recs
                          if str(r.get("status", "")).startswith("error")),
        "audio_minutes": len(ok),
        "total_transients": sum(r["n_events"] for r in ok),
        "total_trains": sum(len(r["trains"]) for r in ok),
        "discharge_candidates": n_cand,
        "days_with_class": {k: len(v) for k, v in sorted(cls_days.items())},
    }
    with open(os.path.join(out, "batch_screen_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))

    # figure: daily transient counts (log) + class presence by month
    months = sorted({r["date"][:7] for r in ok})
    classes = sorted(cls_days)
    colors = {"DISCHARGE-LIKE CANDIDATE": "tab:red",
              "odontocete click activity": "tab:blue",
              "echosounder / engineered narrowband ping activity":
                  "tab:green",
              "irregular impulse cluster": "tab:gray",
              "metronomic train, unclassified": "tab:purple"}
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8),
                                   gridspec_kw={"height_ratios": [2, 1.4]})
    dates = [datetime.date.fromisoformat(r["date"]) for r in ok]
    ax1.semilogy([d for d, r in zip(dates, ok)],
                 [max(r["n_events"], 0.5) for r in ok],
                 ".", ms=4, color="0.4", label="transients per 60 s sample")
    for r, d in zip(ok, dates):
        for t in r["trains"]:
            key = t["class"].split(" (")[0]
            ax1.plot([d], [max(r["n_events"], 0.5)], "o", ms=7,
                     mfc="none", mec=colors.get(key, "k"), mew=1.4)
    ax1.set_ylabel("transients per daily 60 s sample (log)")
    ax1.set_title(f"Archive screen: {len(ok)} days, "
                  f"{summary['total_transients']:,} transients, "
                  f"{summary['total_trains']} trains, "
                  f"{n_cand} discharge-template candidates")
    handles = [plt.Line2D([], [], marker="o", ls="", mfc="none", mec=c,
                          label=k) for k, c in colors.items()
               if k in {c2.split(" (")[0] for c2 in
                        [t["class"] for r in ok for t in r["trains"]]}
               or k.startswith("DISCHARGE")]
    ax1.legend(handles=handles, fontsize=8, loc="upper left")

    grid = np.zeros((len(classes), len(months)))
    for r in ok:
        mi = months.index(r["date"][:7])
        for t in r["trains"]:
            grid[classes.index(t["class"].split(" (")[0]), mi] += 1
    m = ax2.imshow(grid, aspect="auto", cmap="magma")
    ax2.set_xticks(range(len(months)))
    ax2.set_xticklabels(months, rotation=45, ha="right", fontsize=7)
    ax2.set_yticks(range(len(classes)))
    ax2.set_yticklabels(classes, fontsize=8)
    fig.colorbar(m, ax=ax2, label="trains per month")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_batch_screen.png"), dpi=150)
    print(f"wrote {figdir}/fig_batch_screen.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=2022)
    ap.add_argument("--seconds", type=int, default=60)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--tmp", default=os.environ.get(
        "BATCH_TMP", "data/batch_tmp"))
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--tag", default="", help="output filename suffix, "
                    "e.g. _v2 for the chain-mining rescreen")
    args = ap.parse_args()
    if args.aggregate:
        globals()["TAG"] = args.tag
        aggregate()
    else:
        run_year(args.year, args.seconds, args.workers,
                 f"results/batch_screen_{args.year}{args.tag}.jsonl",
                 args.tmp)


if __name__ == "__main__":
    main()
