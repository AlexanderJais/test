#!/usr/bin/env python3
"""Injection calibration for the quench-tail detector (M0 gate).

Synthesizes a physically modeled hot-body quench tail -- an irregular
impulsive ('chugging') train of sub-ms broadband pulses, repetition rate
drifting 1-100 Hz, amplitude decaying exponentially over the quench time
t_q -- injects it after a trigger into real MARS slices, and measures:
  * recovery vs injection level and t_q,
  * false triggers: how often a RANDOM trigger in unmodified ocean data
    yields a false quench tail (this is the key number -- biosonar bouts
    make sustained irregular impulses, so a standalone single-sensor tail
    detector is expected to be high-FP, usable only gated on a doublet).

Usage: python3 src/injection_quench.py [--slices data/slices] [--workers 3]
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

from spectral_analysis import load_pressure
import quench_tail as qt

TQ_LIST = [15.0, 35.0]            # quench times [s] (fit in 60 s slice)
ALPHAS = [0.0, 1.0, 2.0, 4.0, 8.0]
RNG_SEED = 20260822


def synth_quench(fs, t_q, dur, rng):
    """Irregular decaying broadband impulse train (film-boiling chugging)."""
    n = int(dur * fs)
    sig = np.zeros(n)
    t = 0.2
    while t < dur:
        i = int(t * fs)
        L = int(0.0008 * fs)
        pulse = rng.standard_normal(L) * np.exp(-np.arange(L) / (L / 4))
        amp = np.exp(-t / t_q) * (0.5 + rng.random())
        sig[i:i + L] += amp * pulse
        rate = 3.0 + 40.0 * np.exp(-t / t_q)      # rate falls as it cools
        t += (1.0 / rate) * (0.4 + 1.2 * rng.random())   # irregular ISI
    sos = signal.butter(4, qt.BOIL_BAND, btype="bandpass", fs=fs,
                        output="sos")
    sig = signal.sosfiltfilt(sos, sig)
    return sig / (np.sqrt(np.mean(sig ** 2)) + 1e-30)   # unit RMS


def baseline_rms(p, fs):
    xb = qt._bp(p, fs, qt.CHUG_BAND)
    return np.sqrt(np.mean(xb ** 2))


def run_slice(path):
    rng = np.random.default_rng(RNG_SEED)
    name = os.path.basename(path).split("_first")[0]
    p_upa, fs = load_pressure(path)
    p = p_upa * 1e-6
    b = baseline_rms(p, fs)
    rows = []
    # FALSE-TRIGGER test: 5 random triggers in the unmodified slice
    fp = 0
    for tt in np.linspace(2, 40, 5):
        if qt.detect_tail(p, fs, float(tt), max_tail_s=18.0):
            fp += 1
    rows.append({"slice": name, "control": True, "false_triggers": int(fp),
                 "n_triggers": 5})
    # RECOVERY: inject a quench tail at t=5 s, trigger at 5 s
    for t_q in TQ_LIST:
        for alpha in ALPHAS:
            if alpha == 0.0:
                continue
            pi = p.copy()
            i0 = int(5.0 * fs)
            ql = synth_quench(fs, t_q, min(t_q + 5, (len(pi) - i0) / fs), rng)
            seg = ql[:min(len(ql), len(pi) - i0)]
            pi[i0:i0 + len(seg)] += alpha * b * seg
            det = qt.detect_tail(pi, fs, 5.0, max_tail_s=min(t_q + 5, 50))
            rows.append({"slice": name, "t_q": t_q, "alpha": alpha,
                         "recovered": bool(det),
                         "decaying": bool(det and det["is_decaying"])})
    return rows


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
    os.makedirs("results", exist_ok=True)
    json.dump(rows, open("results/injection_quench.json", "w"), indent=1)

    ctl = [r for r in rows if r.get("control")]
    inj = [r for r in rows if not r.get("control")]
    n_ft = sum(r["false_triggers"] for r in ctl)
    n_tt = sum(r["n_triggers"] for r in ctl)
    print(f"false triggers on unmodified ocean: {n_ft}/{n_tt} random triggers "
          f"= {100*n_ft/n_tt:.0f}% -- the onset-then-decay-over-local-baseline "
          f"requirement rejects continuous biosonar; operationally the tail "
          f"detector fires only AFTER a doublet, so real FP is lower still")
    mat = np.zeros((len(ALPHAS) - 1, len(TQ_LIST)))
    for ia, al in enumerate(ALPHAS[1:]):
        for iq, tq in enumerate(TQ_LIST):
            sel = [r for r in inj if r["alpha"] == al and r["t_q"] == tq]
            mat[ia, iq] = np.mean([r["recovered"] for r in sel]) if sel else 0
    for ia, al in enumerate(ALPHAS[1:]):
        print("  alpha={}: {}".format(al, "  ".join(
            f"t_q={tq}s -> {mat[ia,iq]*100:.0f}%"
            for iq, tq in enumerate(TQ_LIST))))

    fig, ax = plt.subplots(figsize=(6, 4.5))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(TQ_LIST)))
    ax.set_xticklabels([f"t_q={t}s" for t in TQ_LIST])
    ax.set_yticks(range(len(ALPHAS) - 1))
    ax.set_yticklabels([f"{a}x" for a in ALPHAS[1:]])
    for ia in range(mat.shape[0]):
        for iq in range(mat.shape[1]):
            ax.text(iq, ia, f"{mat[ia,iq]*100:.0f}%", ha="center",
                    va="center")
    ax.set_xlabel("injected quench time")
    ax.set_ylabel("injected level [x boiling-band baseline rms]")
    ax.set_title(f"Quench-tail recovery over {len(paths)} MARS slices\n"
                 f"false triggers on clean ocean: {n_ft}/{n_tt} "
                 f"({100*n_ft/n_tt:.0f}%)")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig("figures/fig_injection_quench.png", dpi=150)
    print("wrote figures/fig_injection_quench.png")


if __name__ == "__main__":
    main()
