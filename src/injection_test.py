#!/usr/bin/env python3
"""Measure the false-negative rate of the discharge-template screen by
signal injection.

"We found no discharge trains" is only meaningful with a measured
detection threshold: this script synthesizes physically modeled
discharge pulse trains (fast broadband shock + bubble-oscillation echoes
at 0.8/1.6 ms, machine-regular timing with 0.3% jitter), injects them
into REAL 60 s MARS slices spanning quiet ambient, biologic click
storms, and echosounder activity, then runs the exact batch-screen
pipeline and asks: was a DISCHARGE-LIKE train recovered at the injected
repetition rate?

Injection amplitudes are expressed as multiples (alpha) of each slice's
own detection threshold (median + 10*MAD of the band envelope), so
results transfer across noise conditions. A zero-amplitude control run
per slice measures the false-positive side.

Outputs results/injection_test.json and figures/fig_injection.png
(recovery matrix + received-level threshold -> detection-range curve).

Usage:
    python3 src/injection_test.py [--slices data/slices] [--workers 3]
"""

import argparse
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plasma_transient_search import (detect, load_pressure_pa,
                                     screen_events, DET_BAND)

ALPHAS = [0.0, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
RATES_HZ = [0.5, 2.0, 8.0]
IPI_JITTER = 0.003
FEATURE_CAP = 1000
RNG_SEED = 20260821


def make_pulse(fs, rng):
    """Discharge pulse model: 0.1 ms shock burst + bubble echoes."""
    n = int(0.004 * fs)
    t = np.arange(n) / fs
    burst = rng.standard_normal(n) * np.exp(-((t - 0.0002) / 8e-5) ** 2)
    for lag, amp in ((0.0008, 0.4), (0.0016, 0.15)):
        k = int(lag * fs)
        burst[k:] += amp * rng.standard_normal(n - k) * \
            np.exp(-((t[:n - k] - 0.0002) / 8e-5) ** 2)
    sos = signal.butter(4, (DET_BAND[0], 110000), btype="bandpass",
                        fs=fs, output="sos")
    p = signal.sosfiltfilt(sos, burst)
    return p / np.abs(signal.hilbert(p)).max()


def screen(p, fs):
    """The v2 batch-screen pipeline: gap trains + periodicity-mined
    chains (see plasma_transient_search.screen_events)."""
    x, env, peaks, heights = detect(p, fs)
    gap, chains = screen_events(x, env, fs, peaks, heights,
                                feature_cap=FEATURE_CAP)
    return gap + chains


def run_slice(path):
    rng = np.random.default_rng(RNG_SEED)
    name = os.path.basename(path).split("_first")[0]
    p, fs = load_pressure_pa(path)
    sos = signal.butter(6, DET_BAND, btype="bandpass", fs=fs, output="sos")
    env0 = np.abs(signal.hilbert(signal.sosfiltfilt(sos, p)))
    thr = float(np.median(env0) + 10 * np.median(
        np.abs(env0 - np.median(env0))))
    pulse = make_pulse(fs, rng)
    rows = []
    for rate in RATES_HZ:
        for alpha in ALPHAS:
            if alpha == 0.0 and rate != RATES_HZ[0]:
                continue                       # one control per slice
            pi = p.copy()
            if alpha > 0:
                t = 0.5
                while t < len(p) / fs - 0.1:
                    i = int(t * fs)
                    amp = alpha * thr * (1 + 0.02 * rng.standard_normal())
                    seg = pulse[:min(len(pulse), len(pi) - i)]
                    pi[i:i + len(seg)] += amp * seg
                    t += (1 + IPI_JITTER * rng.standard_normal()) / rate
            trains = screen(pi, fs)
            hit = any(t["class"].startswith("DISCHARGE")
                      and abs(t["rate_hz"] - rate) / rate < 0.15
                      for t in trains) if alpha > 0 else None
            fp = sum(1 for t in trains if t["class"].startswith("DISCHARGE")
                     ) if alpha == 0.0 else None
            rows.append({"slice": name, "rate_hz": rate, "alpha": alpha,
                         "thr_pa": round(thr, 4), "recovered": hit,
                         "false_pos": fp})
            print(f"  {name} rate={rate} alpha={alpha}: "
                  f"{'FP=' + str(fp) if alpha == 0 else 'hit' if hit else 'miss'}",
                  flush=True)
    return rows


def detection_range_km(rl_pa, sl_db):
    """Range at which received level falls to rl_pa, for source level
    sl_db (dB re 1 uPa @ 1 m), spherical spreading + 10 dB/km absorption
    (~40 kHz)."""
    rl_db = 20 * np.log10(rl_pa * 1e6)
    r = np.logspace(0, 5.7, 4000)                    # 1 m .. 500 km
    tl = 20 * np.log10(r) + 10.0 * (r / 1000.0)
    ok = np.nonzero(sl_db - tl >= rl_db)[0]
    return float(r[ok[-1]] / 1000.0) if len(ok) else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", default=os.environ.get(
        "PACIFIC_SOUND_DIR", "data/slices"))
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    paths = sorted(glob.glob(os.path.join(args.slices, "*.wav")))
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for rr in ex.map(run_slice, paths):
            rows += rr
    with open("results/injection_test.json", "w") as fh:
        json.dump(rows, fh, indent=1)

    inj = [r for r in rows if r["alpha"] > 0]
    ctl = [r for r in rows if r["alpha"] == 0.0]
    n_fp = sum(r["false_pos"] for r in ctl)
    # recovery matrix alpha x rate
    mat = np.zeros((len(ALPHAS) - 1, len(RATES_HZ)))
    for ia, a in enumerate(ALPHAS[1:]):
        for ir, rt in enumerate(RATES_HZ):
            sel = [r for r in inj if r["alpha"] == a and r["rate_hz"] == rt]
            mat[ia, ir] = np.mean([r["recovered"] for r in sel])
    thr_med = float(np.median([r["thr_pa"] for r in ctl]))
    # threshold: lowest alpha with >=90% overall recovery
    overall = mat.mean(axis=1)
    a90 = next((a for a, v in zip(ALPHAS[1:], overall) if v >= 0.9), None)
    print(f"\ncontrols: {n_fp} false positives in {len(ctl)} clean runs")
    print(f"median detection threshold: {thr_med:.3f} Pa; "
          f"90% recovery at alpha={a90}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    im = ax1.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax1.set_xticks(range(len(RATES_HZ)))
    ax1.set_xticklabels([f"{r} Hz" for r in RATES_HZ])
    ax1.set_yticks(range(len(ALPHAS) - 1))
    ax1.set_yticklabels([f"{a}x" for a in ALPHAS[1:]])
    for ia in range(mat.shape[0]):
        for ir in range(mat.shape[1]):
            ax1.text(ir, ia, f"{mat[ia, ir]*100:.0f}%", ha="center",
                     va="center", fontsize=9)
    ax1.set_xlabel("injected repetition rate")
    ax1.set_ylabel("injected amplitude [x detection threshold]")
    ax1.set_title(f"Recovery rate over {len(paths)} real slices "
                  f"(controls: {n_fp} FP)", fontsize=10)
    fig.colorbar(im, ax=ax1, shrink=0.85)
    sls = np.arange(160, 231, 1)
    for a, ls in ((a90 or 2.0, "-"), (5.0, "--")):
        rng_km = [detection_range_km(a * thr_med, sl) for sl in sls]
        ax2.plot(sls, rng_km, ls, label=f"received level {a}x threshold "
                 f"({a*thr_med:.2f} Pa)")
    ax2.set_yscale("log")
    ax2.set_xlabel("source level [dB re 1 µPa @ 1 m]")
    ax2.set_ylabel("detection range [km]")
    ax2.set_title("Implied detection range (spherical spreading + "
                  "10 dB/km @ 40 kHz)", fontsize=10)
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8)
    fig.suptitle("Injection test: measured sensitivity of the "
                 "discharge-template screen", fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/fig_injection.png", dpi=150)
    print("wrote figures/fig_injection.png")


if __name__ == "__main__":
    main()
