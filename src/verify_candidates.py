#!/usr/bin/env python3
"""Adjudicate discharge-template candidate trains found by the batch
screen.

A metronomic broadband impulse train has three plausible explanations:
an electrical-discharge source (the template's target), an engineered
broadband pinger, or an unusually regular biosonar click train. The
per-pulse physics that separates them:

* amplitude stability — a fixed discharge/pinger transmits identical
  pulses, so received levels are stable; an echolocating animal sweeps
  its beam, so received click levels vary strongly pulse to pulse;
* waveform/spectral self-similarity — machine pulses are clones
  (pairwise spectral correlation ~1); animal clicks decorrelate as the
  beam aspect changes;
* echo-lag consistency — a discharge's bubble-oscillation echo has a
  fixed lag every pulse; a surface-reflection ghost varies;
* context — pulses physically identical to surrounding jittery biosonar
  in the same minute argue for the same (biologic) source.

For each candidate: re-fetch the slice, re-detect, isolate the
candidate train by its logged start time and rate, and measure all of
the above. Produces one diagnostic figure per candidate and a verdict
JSON.

Usage:
    python3 src/verify_candidates.py
"""

import json
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fetch_data import fetch_slice
from plasma_transient_search import detect, load_pressure_pa

CANDIDATES = [
    {"id": "A", "year": "2021", "month": "11", "prefix": "20211101_23",
     "t_start_s": 1.148, "rate_hz": 0.506, "n": 14},
    {"id": "B", "year": "2021", "month": "12", "prefix": "20211228_14",
     "t_start_s": 9.195, "rate_hz": 1.555, "n": 16},
]


def isolate_train(times, t0, rate, n):
    """Greedy re-pick of the metronomic train from all detection times."""
    period = 1.0 / rate
    train = [min(times, key=lambda t: abs(t - t0))]
    for _ in range(n - 1):
        pred = train[-1] + period
        cands = [t for t in times if abs(t - pred) < 0.25 * period]
        if not cands:
            break
        train.append(min(cands, key=lambda t: abs(t - pred)))
    return np.array(train)


def analyze(cand, slices_dir, figdir, seconds=60):
    path = fetch_slice(cand["year"], cand["month"], cand["prefix"],
                       seconds, slices_dir)
    p, fs = load_pressure_pa(path)
    x, env, peaks, heights = detect(p, fs)
    times = peaks / fs
    tr = isolate_train(times, cand["t_start_s"], cand["rate_hz"], cand["n"])
    idx = [int(round(t * fs)) for t in tr]

    # per-pulse extracts: 6 ms windows, aligned on envelope peak
    h = int(0.003 * fs)
    wins, amps, specs, echoes = [], [], [], []
    for i in idx:
        k = i - h + int(np.argmax(env[i - h:i + h]))
        w = x[k - h:k + h]
        wins.append(w)
        amps.append(float(env[k]))
        f, s = signal.welch(w, fs=fs, nperseg=512)
        specs.append(s / s.sum())
        e = env[k:k + int(0.01 * fs)] - env[k:k + int(0.01 * fs)].mean()
        ac = signal.correlate(e, e, "full")[len(e) - 1:]
        ac /= ac[0] + 1e-30
        l0 = int(0.0003 * fs)
        echoes.append((l0 + int(np.argmax(ac[l0:]))) / fs * 1e3)
    amps = np.array(amps)
    S = np.array(specs)
    # band dominance: is the pulse energy actually broadband, or confined
    # to the 36-40 kHz echosounder band? (the batch screener's fractional-
    # bandwidth feature measures signal+noise and overestimates bandwidth
    # for low-SNR pulses)
    sos_ping = signal.butter(4, [36000, 40000], btype="bandpass", fs=fs,
                             output="sos")
    sos_high = signal.butter(4, [55000, 95000], btype="bandpass", fs=fs,
                             output="sos")
    ep = np.abs(signal.hilbert(signal.sosfiltfilt(sos_ping, x)))
    eh = np.abs(signal.hilbert(signal.sosfiltfilt(sos_high, x)))
    snr_p, snr_h = [], []
    for i in idx:
        a, b = i - h, i + h
        snr_p.append(ep[a:b].max() / np.median(ep))
        snr_h.append(eh[a:b].max() / np.median(eh))
    ping_snr = float(np.median(snr_p))
    high_snr = float(np.median(snr_h))
    # pairwise spectral correlation
    C = np.corrcoef(S)
    iu = np.triu_indices(len(S), 1)
    metrics = {
        "id": cand["id"], "slice": os.path.basename(path),
        "n_pulses_reacquired": len(tr),
        "ipi_cv": round(float(np.std(np.diff(tr)) / np.mean(np.diff(tr))),
                        4),
        "amp_cv": round(float(np.std(amps) / np.mean(amps)), 3),
        "amp_range_db": round(float(20 * np.log10(amps.max() / amps.min())),
                              1),
        "spec_corr_median": round(float(np.median(C[iu])), 3),
        "echo_lag_cv": round(float(np.std(echoes) / np.mean(echoes)), 3),
        "ping_band_snr": round(ping_snr, 1),
        "high_band_snr": round(high_snr, 1),
    }
    # verdict logic
    if ping_snr > 4 * high_snr and ping_snr > 10:
        # pulse energy confined to the 38 kHz sonar band, nothing above:
        # a narrowband echosounder whose low SNR fooled the batch
        # screener's bandwidth feature
        metrics["verdict"] = ("engineered narrowband ping (38 kHz-class " 
                              "echosounder); batch frac-BW overestimated "
                              "at low SNR")
    elif (metrics["amp_cv"] < 0.15 and metrics["spec_corr_median"] > 0.9
            and metrics["echo_lag_cv"] < 0.15):
        metrics["verdict"] = ("machine-like broadband transmitter "
                              "(discharge or engineered; needs follow-up)")
    else:
        metrics["verdict"] = ("biosonar: click-to-click level/spectrum "
                              "variability of a scanning animal, not a "
                              "fixed transmitter")

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    ax = axes[0, 0]
    ax.plot(range(1, len(tr)), np.diff(tr) * 1e3, "o-")
    ax.set_xlabel("pulse #")
    ax.set_ylabel("inter-pulse interval [ms]")
    ax.set_title(f"IPI series (CV {metrics['ipi_cv']})", fontsize=10)
    ax = axes[0, 1]
    t_ms = (np.arange(2 * h) - h) / fs * 1e3
    for w in wins:
        ax.plot(t_ms, w / np.max(np.abs(w)), lw=0.5, alpha=0.5)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xlabel("time [ms]")
    ax.set_title("all pulses, peak-normalized overlay", fontsize=10)
    ax = axes[1, 0]
    ax.semilogy(range(1, len(amps) + 1), amps, "o-")
    ax.set_xlabel("pulse #")
    ax.set_ylabel("received peak envelope [Pa]")
    ax.set_title(f"levels: CV {metrics['amp_cv']}, range "
                 f"{metrics['amp_range_db']} dB", fontsize=10)
    ax = axes[1, 1]
    f = np.fft.rfftfreq(512, 1 / fs)[:len(S[0])]
    for s in S:
        ax.semilogy(f / 1e3, s, lw=0.5, alpha=0.5)
    ax.set_xlabel("frequency [kHz]")
    ax.set_title(f"per-pulse spectra (median pairwise corr "
                 f"{metrics['spec_corr_median']})", fontsize=10)
    fig.suptitle(f"Candidate {cand['id']} — {os.path.basename(path)}\n"
                 f"verdict: {metrics['verdict']}", fontsize=11)
    fig.tight_layout()
    out = os.path.join(figdir, f"fig_candidate_{cand['id']}.png")
    fig.savefig(out, dpi=150)
    print(json.dumps(metrics, indent=1))
    return metrics


def main():
    slices_dir = "data/batch_tmp"
    os.makedirs(slices_dir, exist_ok=True)
    out = [analyze(c, slices_dir, "figures") for c in CANDIDATES]
    with open("results/candidate_verdicts.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
