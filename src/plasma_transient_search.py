#!/usr/bin/env python3
"""Search full-bandwidth hydrophone data for plasma-discharge transient
signatures.

A craft enveloped in a plasma sheath, or propelled by plasma, must
sustain electrical discharge in water. Underwater discharge acoustics is
well characterized (marine "sparker" seismic sources are exactly this):

* each pulse: a fast broadband shock (µs-scale rise) followed by a
  BUBBLE-OSCILLATION echo — the vapor/plasma cavity re-expands and
  collapses with period T ∝ E^(1/3) / P^(5/6) (sub-ms to a few ms;
  compressed at depth),
* per-pulse spectrum: BROADBAND (fractional bandwidth ~1), unlike the
  narrowband tonal pings of echosounders,
* sustainment: driven by a pulsed power supply -> METRONOMIC repetition
  (inter-pulse interval jitter far below biological trains: snapping
  shrimp are Poisson-random; odontocete inter-click intervals drift
  smoothly and widely as the animal ranges a target),
* a moving source adds slow secular drift of repetition rate (Doppler)
  and received level.

Template for a DISCHARGE-LIKE candidate train, all required:
  n >= 8 pulses, IPI CV < 0.05 (machine-regular), median fractional
  bandwidth > 0.5 (broadband, excludes echosounders), median duration
  < 5 ms (impulsive, excludes tonal bursts).

This script detects transients in the 256 kHz MARS slices (5-120 kHz
band), extracts per-event waveform physics (rise time, duration,
spectral centroid, fractional bandwidth, bubble-echo lag/strength,
amplitude), clusters events into trains, computes train timing
statistics, classifies every train and isolated event, and reports any
template matches.

Usage:
    python3 src/plasma_transient_search.py [--slices data/slices]
"""

import argparse
import glob
import json
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spectral_analysis import load_pressure as _load_upa
from fetch_data import fetch_slice


def load_pressure_pa(path):
    """Calibrated pressure in Pa (part 1's load_pressure returns µPa)."""
    p_upa, fs = _load_upa(path)
    return p_upa * 1e-6, fs

# widen temporal coverage beyond the part-1 survey slices
EXTRA_SLICES = [
    ("2022", "02", "20220214_09"),
    ("2022", "04", "20220404_15"),
    ("2022", "06", "20220616_00"),
    ("2022", "08", "20220808_12"),
    ("2022", "10", "20221010_06"),
    ("2022", "12", "20221220_21"),
]

DET_BAND = (5000.0, 120000.0)
MAD_K = 10.0
MIN_SEP_S = 0.002
EVENT_HALF_S = 0.010          # +/- 10 ms analysis window per event
TRAIN_GAP_S = 2.0
TRAIN_MIN_N = 8


def detect(p, fs):
    sos = signal.butter(6, DET_BAND, btype="bandpass", fs=fs, output="sos")
    x = signal.sosfiltfilt(sos, p)
    env = np.abs(signal.hilbert(x))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    peaks, props = signal.find_peaks(env, height=med + MAD_K * mad,
                                     distance=int(MIN_SEP_S * fs))
    return x, env, peaks, props["peak_heights"]


def event_features(x, env, fs, i):
    """Waveform physics of one transient at sample index i."""
    h = int(EVENT_HALF_S * fs)
    a, b = max(0, i - h), min(len(x), i + h)
    w, e = x[a:b], env[a:b]
    pk = i - a
    peak = e[pk]
    # duration at -20 dB of peak
    above = e >= 0.1 * peak
    lo = pk
    while lo > 0 and above[lo - 1]:
        lo -= 1
    hi = pk
    while hi < len(e) - 1 and above[hi + 1]:
        hi += 1
    dur_ms = (hi - lo) / fs * 1e3
    # rise time 10% -> 90% of peak, searching backwards from the peak
    j10 = lo
    j90 = pk
    while j90 > lo and e[j90 - 1] >= 0.9 * peak:
        j90 -= 1
    rise_us = max(j90 - j10, 1) / fs * 1e6
    # spectrum of the event
    seg = w[max(0, pk - int(0.004 * fs)):pk + int(0.006 * fs)]
    f, pxx = signal.welch(seg, fs=fs, nperseg=min(1024, len(seg)))
    sel = f >= DET_BAND[0]
    fq, s = f[sel], pxx[sel]
    centroid = float(np.sum(fq * s) / np.sum(s))
    cum = np.cumsum(s) / np.sum(s)
    bw = float(fq[np.searchsorted(cum, 0.95)] -
               fq[np.searchsorted(cum, 0.05)])   # 90% energy bandwidth
    # bubble-oscillation echo: envelope autocorrelation 0.3-20 ms lag
    ec = e[pk:pk + int(0.020 * fs)].copy()
    ec -= ec.mean()
    ac = signal.correlate(ec, ec, mode="full")[len(ec) - 1:]
    ac /= ac[0] + 1e-30
    l0 = int(0.0003 * fs)
    if len(ac) > l0 + 2:
        k = l0 + int(np.argmax(ac[l0:]))
        echo_lag_ms, echo_r = k / fs * 1e3, float(ac[k])
    else:
        echo_lag_ms, echo_r = np.nan, 0.0
    return {
        "t_s": round(i / fs, 4),
        "peak_pa": round(float(peak), 2),
        "dur_ms": round(float(dur_ms), 3),
        "rise_us": round(float(rise_us), 1),
        "centroid_hz": round(centroid, 0),
        "frac_bw": round(bw / centroid, 2),
        "echo_lag_ms": round(float(echo_lag_ms), 2),
        "echo_r": round(echo_r, 2),
    }


def classify_train(evs, cv_override=None):
    ipi = np.diff([e["t_s"] for e in evs])
    cv = float(cv_override) if cv_override is not None \
        else float(np.std(ipi) / np.mean(ipi))
    # secular IPI drift (Doppler / range-rate proxy): relative slope
    tt = np.array([e["t_s"] for e in evs[:-1]])
    drift = float(np.polyfit(tt, ipi, 1)[0] / np.mean(ipi)) if len(ipi) > 3 \
        else np.nan
    # tolerate events lacking waveform features (batch mode subsamples
    # feature extraction on very dense click storms; timing uses all)
    med = lambda k: float(np.median([e[k] for e in evs if k in e]))
    stats = {"n": len(evs), "rate_hz": round(1.0 / float(np.mean(ipi)), 3),
             "ipi_cv": round(cv, 3),
             "ipi_drift_per_s": round(drift, 5) if np.isfinite(drift)
             else None,
             "med_dur_ms": round(med("dur_ms"), 3),
             "med_centroid_hz": round(med("centroid_hz"), 0),
             "med_frac_bw": round(med("frac_bw"), 2),
             "med_echo_lag_ms": round(med("echo_lag_ms"), 2)}
    if stats["med_frac_bw"] < 0.2:
        # Narrowband dominates any timing evidence: an engineered tonal
        # ping (e.g. the 38 kHz fisheries-echosounder standard). Timing
        # CV is unreliable here because each ping arrives with surface/
        # bottom echoes that the detector also picks up.
        label = "echosounder / engineered narrowband ping activity"
    elif cv < 0.05 and stats["med_frac_bw"] > 0.5 \
            and stats["med_dur_ms"] < 5:
        label = "DISCHARGE-LIKE CANDIDATE (metronomic, broadband, impulsive)"
    elif cv < 0.05:
        label = "metronomic train, unclassified"
    elif stats["med_centroid_hz"] > 15000 and stats["med_dur_ms"] < 1 \
            and stats["med_frac_bw"] > 0.8:
        # Sub-ms broadband high-frequency clicks; timing jitter is high
        # when several animals' trains overlap in one aggregate.
        label = "odontocete click activity (single or overlapping trains)"
    else:
        label = "irregular impulse cluster (snaps / unclassified)"
    stats["class"] = label
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", default=os.environ.get(
        "PACIFIC_SOUND_DIR", "data/slices"))
    ap.add_argument("--seconds", type=int, default=60)
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()

    for year, month, prefix in EXTRA_SLICES:
        try:
            fetch_slice(year, month, prefix, args.seconds, args.slices)
        except Exception as e:
            print(f"  fetch {prefix} failed ({e}); continuing")

    catalog, gallery = [], []
    for path in sorted(glob.glob(os.path.join(args.slices, "*.wav"))):
        name = os.path.basename(path).split("_first")[0]
        p, fs = load_pressure_pa(path)
        x, env, peaks, heights = detect(p, fs)
        evs = [event_features(x, env, fs, i) for i in peaks]
        # cluster into trains by time gaps
        trains, cur = [], []
        for e in evs:
            if cur and e["t_s"] - cur[-1]["t_s"] > TRAIN_GAP_S:
                trains.append(cur)
                cur = []
            cur.append(e)
        if cur:
            trains.append(cur)
        entry = {"slice": name, "n_events": len(evs), "trains": [],
                 "isolated": 0}
        for tr in trains:
            if len(tr) >= TRAIN_MIN_N:
                st = classify_train(tr)
                st["t_start_s"] = tr[0]["t_s"]
                entry["trains"].append(st)
            else:
                entry["isolated"] += len(tr)
        catalog.append(entry)
        if len(evs):
            k = int(np.argmax([e["peak_pa"] for e in evs]))
            i = peaks[k]
            gallery.append((name, evs[k],
                            x[max(0, i - fs // 250):i + fs // 125].copy(),
                            fs))
        tl = "; ".join(f"{t['class'].split(' (')[0]} n={t['n']} "
                       f"cv={t['ipi_cv']}" for t in entry["trains"]) or "-"
        print(f"{name}: {len(evs)} transients, {len(entry['trains'])} "
              f"trains [{tl}]")

    n_cand = sum(1 for c in catalog for t in c["trains"]
                 if t["class"].startswith("DISCHARGE"))
    print(f"\ndischarge-template candidate trains: {n_cand}")
    with open(os.path.join(args.out, "transient_catalog.json"), "w") as fh:
        json.dump({"template": {"n_min": TRAIN_MIN_N, "ipi_cv_max": 0.05,
                                "frac_bw_min": 0.5, "dur_ms_max": 5},
                   "n_candidates": n_cand, "slices": catalog}, fh, indent=1)

    # ---- figures -------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    allev = [(c["slice"], t) for c in catalog for t in c["trains"]]
    colors = {"DISCHARGE": "tab:red", "odontocete": "tab:blue",
              "echosounder": "tab:green", "irregular": "tab:gray"}
    for name, t in allev:
        key = next(k for k in colors if t["class"].startswith(k))
        ax1.scatter(t["rate_hz"], t["ipi_cv"], s=30 + t["n"] / 3,
                    color=colors[key], alpha=0.75)
        ax2.scatter(t["med_centroid_hz"] / 1e3, t["med_frac_bw"],
                    s=30 + t["n"] / 3, color=colors[key], alpha=0.75)
    ax1.axhspan(1e-4, 0.05, alpha=0.12, color="tab:red")
    ax1.text(0.05, 0.035, "machine-regular zone: a sustained discharge\n"
             "(plasma sheath drive) must sit here", fontsize=8,
             color="tab:red", transform=ax1.get_yaxis_transform())
    ax1.set_yscale("log")
    ax1.set_xlabel("train repetition rate [Hz]")
    ax1.set_ylabel("inter-pulse-interval CV (timing jitter)")
    ax2.axhspan(0.5, 2.0, alpha=0.10, color="tab:red")
    ax2.text(0.55, 0.9, "broadband zone (discharge/cavitation);\n"
             "echosounder pings fall below", fontsize=8, color="tab:red",
             transform=ax2.get_xaxis_transform())
    ax2.set_xlabel("median spectral centroid [kHz]")
    ax2.set_ylabel("median fractional bandwidth")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c,
                          label={"DISCHARGE": "discharge-like candidate",
                                 "odontocete": "odontocete click train",
                                 "echosounder": "echosounder/ping",
                                 "irregular": "irregular/snaps"}[k])
               for k, c in colors.items()]
    ax1.legend(handles=handles, fontsize=8)
    fig.suptitle("Transient trains in 12 x 60 s of MARS 256 kHz data vs "
                 "the plasma-discharge template", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_transient_features.png"),
                dpi=150)

    ng = min(len(gallery), 6)
    order = sorted(range(len(gallery)),
                   key=lambda i: -gallery[i][1]["peak_pa"])[:ng]
    fig, axes = plt.subplots(ng, 1, figsize=(11, 2.1 * ng))
    for ax, gi in zip(np.atleast_1d(axes), order):
        name, ev, w, fs = gallery[gi]
        t = (np.arange(len(w)) - len(w) / 3) / fs * 1e3
        ax.plot(t, w, lw=0.5, color="k")
        ax.set_ylabel("Pa", fontsize=8)
        ax.set_title(f"{name}  t={ev['t_s']} s — dur {ev['dur_ms']} ms, "
                     f"rise {ev['rise_us']} µs, centroid "
                     f"{ev['centroid_hz']/1e3:.0f} kHz, frac-BW "
                     f"{ev['frac_bw']}, echo lag {ev['echo_lag_ms']} ms "
                     f"(r={ev['echo_r']})", fontsize=8)
    np.atleast_1d(axes)[-1].set_xlabel("time [ms]")
    fig.suptitle("Strongest transient per slice — waveform physics",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_transient_gallery.png"),
                dpi=150)
    print(f"wrote {args.figdir}/fig_transient_features.png and "
          "fig_transient_gallery.png")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Periodicity mining (PRI deinterleaving). The gap-based train clustering
# above describes the transient population well, but injection testing
# exposed two structural blind spots for CANDIDATE DETECTION: trains with
# period >= TRAIN_GAP_S fragment into undersized pieces, and trains
# interleaved with biosonar/echosounder activity inherit huge timing CV.
# A periodic source must instead be sought as a periodic CHAIN inside the
# mixed event stream, tolerant of interlopers and missing pulses.
# ---------------------------------------------------------------------------

def find_periodic_chains(times, period_range=(0.05, 5.0), n_min=8,
                         rel_tol=0.012, max_skip=1, max_periods=30):
    """Extract metronomic chains from a mixed stream of event times.

    1. histogram all pairwise time differences within period_range,
    2. take the strongest histogram peaks as candidate periods,
    3. greedily grow chains event->event at the candidate period
       (tolerance max(rel_tol*T, 8 ms), riding over up to max_skip
       missing pulses),
    4. keep chains with >= n_min members whose per-step period estimates
       have CV < 0.05; dedupe chains sharing >50% of members.

    Returns a list of (member_time_array, period_s).
    """
    t = np.sort(np.asarray(times, dtype=float))
    n = len(t)
    if n < n_min:
        return []
    # pairwise diffs within range (bounded per event by period_range[1])
    diffs = []
    j0 = 0
    for i in range(n):
        j = i + 1
        while j < n and t[j] - t[i] <= period_range[1]:
            if t[j] - t[i] >= period_range[0]:
                diffs.append(t[j] - t[i])
            j += 1
    if len(diffs) < n_min - 1:
        return []
    diffs = np.array(diffs)
    binw = 0.002
    hist, edges = np.histogram(diffs, bins=int(
        (period_range[1] - period_range[0]) / binw))
    # windowed SUM (not average): a period concentrated in one 2 ms bin
    # must not be diluted below threshold by its empty neighbors
    score = np.convolve(hist, np.ones(5), mode="same")
    order = np.argsort(score)[::-1]
    periods, used_bins = [], np.zeros(len(score), bool)
    for b in order:
        if len(periods) >= max_periods or score[b] < (n_min - 1) * 0.6:
            break
        if used_bins[max(0, b - 3):b + 4].any():
            continue
        used_bins[b] = True
        periods.append(0.5 * (edges[b] + edges[b + 1]))

    chains = []
    for T in periods:
        tol = max(rel_tol * T, 0.005)
        # seeds: events that begin a T-spaced pair
        seeds = [i for i in range(n - 1)
                 if np.any(np.abs(t[i + 1:] - t[i] - T) < tol)][:60]
        best = None
        for s in seeds:
            chain = [t[s]]
            while True:
                nxt = None
                for skip in range(1, max_skip + 2):
                    pred = chain[-1] + skip * T
                    if pred > t[-1] + tol:
                        break
                    j = np.searchsorted(t, pred)
                    for jj in (j - 1, j):
                        if 0 <= jj < n and abs(t[jj] - pred) < tol * skip \
                                and t[jj] > chain[-1]:
                            nxt = t[jj]
                            break
                    if nxt is not None:
                        break
                if nxt is None:
                    break
                chain.append(nxt)
            if best is None or len(chain) > len(best):
                best = chain
        if best is None or len(best) < n_min:
            continue
        steps = np.diff(best)
        per_step = steps / np.round(steps / T)
        if np.std(per_step) / np.mean(per_step) >= 0.05:
            continue
        chains.append((np.array(best), float(np.mean(per_step))))
    # dedupe by member overlap, longest first
    chains.sort(key=lambda c: -len(c[0]))
    kept = []
    for ch, T in chains:
        if all(len(np.intersect1d(np.round(ch, 4), np.round(k0, 4)))
               <= 0.5 * len(ch) for k0, _ in kept):
            kept.append((ch, T))
    return kept


def _max_chance_chain(times, n_surrogates=3, seed=0):
    """Longest chain the builder produces on Poisson surrogates with the
    same event density: the chance floor for chain significance."""
    rng = np.random.default_rng(seed)
    t = np.asarray(times)
    if len(t) < 2:
        return 0
    span = t[-1] - t[0]
    m = 0
    for _ in range(n_surrogates):
        surr = np.sort(t[0] + rng.uniform(0, span, len(t)))
        for ch, _T in find_periodic_chains(surr):
            m = max(m, len(ch))
    return m


def mine_chains(times, heights, n_min=8):
    """Amplitude-tiered periodic-chain mining with chance-level control.

    Dense biologic activity swamps the pairwise-interval histogram, so
    chains are mined in amplitude tiers (all events, top-600, top-150 by
    peak height) where a strong periodic source stands out. Each tier's
    chains must exceed that tier's Poisson chance floor.
    """
    times = np.asarray(times)
    heights = np.asarray(heights)
    tiers = [np.arange(len(times))]
    for k in (600, 150):
        if len(times) > k:
            tiers.append(np.sort(np.argsort(heights)[-k:]))
    # amplitude BANDS: a machine train repeats at near-constant received
    # level, so it concentrates inside a narrow band even when biologics
    # and strong pings dominate the totals. Slide a 6 dB window in 3 dB
    # steps across the height distribution.
    logh = 20 * np.log10(np.maximum(heights, 1e-12))
    for lo in np.arange(np.floor(logh.min()), logh.max(), 3.0):
        band = np.nonzero((logh >= lo) & (logh < lo + 6.0))[0]
        if n_min <= len(band) <= 2000:
            tiers.append(band)
    kept = []
    for tier in tiers:
        t = times[tier]
        found = find_periodic_chains(t, n_min=n_min)
        if not found:
            continue
        floor = _max_chance_chain(t)
        for ch, T in found:
            if len(ch) < max(n_min, floor + 2):
                continue
            if all(len(np.intersect1d(np.round(ch, 4), np.round(k0, 4)))
                   <= 0.5 * len(ch) for k0, _ in kept):
                kept.append((ch, T))
    return kept


def screen_events(x, env, fs, peaks, heights, feature_cap=1000):
    """Shared screening core: gap-clustered trains (population
    description) plus periodicity-mined chains (candidate detection).
    Returns (gap_train_stats, chain_stats)."""
    feat_idx = set(np.linspace(0, len(peaks) - 1,
                               min(len(peaks), feature_cap)).astype(int)
                   ) if len(peaks) else set()
    evs = [event_features(x, env, fs, i) if j in feat_idx
           else {"t_s": round(i / fs, 4)}
           for j, i in enumerate(peaks)]
    gap_trains, cur = [], []
    for e in evs:
        if cur and e["t_s"] - cur[-1]["t_s"] > TRAIN_GAP_S:
            gap_trains.append(cur)
            cur = []
        cur.append(e)
    if cur:
        gap_trains.append(cur)
    gap_stats = [classify_train(t) | {"t_start_s": t[0]["t_s"]}
                 for t in gap_trains if len(t) >= TRAIN_MIN_N]

    chain_stats = []
    times = peaks / fs
    for ch, T in mine_chains(times, heights):
        idx = np.searchsorted(times, ch - 1e-4)
        members = [event_features(x, env, fs, int(peaks[i]))
                   for i in idx[:200]]
        steps = np.diff(ch)
        per_step = steps / np.round(steps / T)
        chain_cv = float(np.std(per_step) / np.mean(per_step))
        st = classify_train(members, cv_override=chain_cv)
        st["t_start_s"] = round(float(ch[0]), 3)
        st["n"] = len(ch)
        st["rate_hz"] = round(1.0 / T, 3)
        st["path"] = "chain"
        chain_stats.append(st)
    return gap_stats, chain_stats
