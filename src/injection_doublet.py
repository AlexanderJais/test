#!/usr/bin/env python3
"""Injection calibration for the entry-doublet detector (M0 gate).

No detector proceeds past M0 without a measured sensitivity curve. This
synthesizes physically modeled water-entry doublets -- slam shock (broadband,
low corner ~240/a Hz) + quiet interlude + pinch-off pulse at
t_p = 2*sqrt(a/g) + Minnaert ring f0 = 3.26*sqrt(1+z/10)/R_b -- injects them
into REAL MARS slices at controlled amplitude (multiples of each slice's own
broadband-impulse detection threshold), across quiet / biosonar / echosounder
conditions, then measures:
  * recovery rate vs body size a (=> t_p) and injection amplitude,
  * false-positive rate on zero-amplitude control runs.

Outputs results/injection_doublet.json and figures/fig_injection_doublet.png.

Usage:
    python3 src/injection_doublet.py [--slices data/slices] [--workers 3]
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
import entry_doublet as ed

G = 9.81
BODY_A = [0.05, 0.15, 0.5]        # m  -> t_p = 0.14, 0.25, 0.45 s
ALPHAS = [0.0, 1.0, 1.5, 2.0, 3.0, 5.0]
DEPTH_M = 891.0                   # MARS depth for the Minnaert ring
RNG_SEED = 20260822


def synth_doublet(fs, a_m, rng):
    """Physically shaped entry doublet at unit peak amplitude."""
    tp = 2.0 * np.sqrt(a_m / G)                 # pinch delay
    total = tp + 1.0
    n = int(total * fs)
    sig = np.zeros(n)

    def shock(centeridx, band_lo, band_hi, amp):
        L = int(0.004 * fs)
        t = np.arange(L) / fs
        burst = rng.standard_normal(L) * np.exp(-((t - 0.0003) / 1.2e-4) ** 2)
        sos = signal.butter(4, (band_lo, band_hi), btype="bandpass",
                            fs=fs, output="sos")
        b = signal.sosfiltfilt(sos, burst)
        b /= np.abs(b).max() + 1e-30
        s, e = centeridx, centeridx + L
        sig[s:min(e, n)] += amp * b[:min(e, n) - s]

    # slam: low-pass with spectral corner f_c ~ 240/a Hz (energy f_c/3..f_c*3)
    fc = max(200.0, 240.0 / a_m)
    shock(int(0.02 * fs), max(120.0, fc / 3), min(fc * 3, fs / 2 * 0.9), 1.0)
    ip = int((0.02 + tp) * fs)
    shock(ip, 150.0, 1400.0, 0.8)                                      # pinch
    # Minnaert ring after pinch
    Rb = max(0.05, a_m)
    f0 = 3.26 * np.sqrt(1 + DEPTH_M / 10.0) / Rb
    tr = np.arange(int(0.6 * fs)) / fs
    ring = 0.2 * np.sin(2 * np.pi * f0 * tr) * np.exp(-tr / 0.15)
    sig[ip:ip + len(ring)] += ring[:max(0, n - ip)][:len(sig[ip:ip + len(ring)])]
    # quiet interlude is just the gap (near zero) between the pulses
    return sig


def thr_of(p, fs):
    x = ed._bandpass(p, fs, ed.DET_BAND)
    env = np.abs(signal.hilbert(x))
    return float(np.median(env) + ed.MAD_K * np.median(np.abs(env - np.median(env))))


def run_slice(path):
    rng = np.random.default_rng(RNG_SEED)
    name = os.path.basename(path).split("_first")[0]
    p_upa, fs = load_pressure(path)
    p = p_upa * 1e-6
    thr = thr_of(p, fs)
    rows = []
    for a_m in BODY_A:
        tp = 2.0 * np.sqrt(a_m / G)
        for alpha in ALPHAS:
            if alpha == 0.0 and a_m != BODY_A[0]:
                continue                          # one control per slice
            pi = p.copy()
            # inject 3 well-separated doublets in the slice
            inj_times = []
            if alpha > 0:
                dl = synth_doublet(fs, a_m, rng)
                for t0 in (8.0, 25.0, 45.0):
                    i0 = int(t0 * fs)
                    seg = dl[:min(len(dl), len(pi) - i0)]
                    pi[i0:i0 + len(seg)] += alpha * thr * seg
                    inj_times.append(t0 + 0.02)
            dets = ed.detect(pi, fs)
            n_imp = int(len(ed.find_impulses(p, fs)[2]))   # ambient density
            if alpha == 0.0:
                rows.append({"slice": name, "a_m": a_m, "alpha": 0.0,
                             "false_pos": len(dets),
                             "false_pos_A": sum(1 for d in dets
                                                if d["grade"] == "A"),
                             "n_imp": n_imp})
            else:
                # recovered if a detection's slam lands within 0.1 s of an
                # injected slam AND tp matches within 20%
                rec = 0
                for it in inj_times:
                    ok = any(abs(d["t_slam_s"] - it) < 0.1
                             and abs(d["tp_s"] - tp) / tp < 0.2 for d in dets)
                    rec += int(ok)
                rows.append({"slice": name, "a_m": a_m, "alpha": alpha,
                             "recovered": rec, "injected": len(inj_times),
                             "n_imp": n_imp})
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
    json.dump(rows, open("results/injection_doublet.json", "w"), indent=1)

    ctl = [r for r in rows if r["alpha"] == 0.0]
    inj = [r for r in rows if r["alpha"] > 0]
    n_fp = sum(r["false_pos"] for r in ctl)
    n_fpA = sum(r.get("false_pos_A", 0) for r in ctl)
    # classify slices quiet vs storm by ambient impulse count (median split)
    imp_by_slice = {r["slice"]: r["n_imp"] for r in ctl}
    med = np.median(list(imp_by_slice.values()))
    quiet = {s for s, v in imp_by_slice.items() if v <= med}
    print(f"controls: {n_fp} false-positive doublets ({n_fpA} grade-A) in "
          f"{len(ctl)} clean 60 s runs across {len(paths)} slices "
          f"-> {n_fp/len(paths):.1f} FP/min single-sensor")
    for label, set_ in [("QUIET slices", quiet),
                         ("STORM slices", set(imp_by_slice) - quiet)]:
        sub = [r for r in inj if r["slice"] in set_]
        for al in [1.5, 3.0, 5.0]:
            ss = [r for r in sub if r["alpha"] == al]
            tot = sum(r["injected"] for r in ss); rec = sum(r["recovered"] for r in ss)
            print(f"  {label} alpha={al}: {100*rec/tot:.0f}% recovered "
                  f"({rec}/{tot})")
    mat = np.zeros((len(ALPHAS) - 1, len(BODY_A)))
    for ia, al in enumerate(ALPHAS[1:]):
        for ib, a_m in enumerate(BODY_A):
            sel = [r for r in inj if r["alpha"] == al and r["a_m"] == a_m]
            tot = sum(r["injected"] for r in sel)
            rec = sum(r["recovered"] for r in sel)
            mat[ia, ib] = rec / tot if tot else np.nan
    for ia, al in enumerate(ALPHAS[1:]):
        print(f"  alpha={al}: " + "  ".join(
            f"a={a}m t_p={2*np.sqrt(a/G):.2f}s -> {mat[ia,ib]*100:.0f}%"
            for ib, a in enumerate(BODY_A)))

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(BODY_A)))
    ax.set_xticklabels([f"a={a} m\nt_p={2*np.sqrt(a/G):.2f} s" for a in BODY_A])
    ax.set_yticks(range(len(ALPHAS) - 1))
    ax.set_yticklabels([f"{a}x" for a in ALPHAS[1:]])
    for ia in range(mat.shape[0]):
        for ib in range(mat.shape[1]):
            ax.text(ib, ia, f"{mat[ia,ib]*100:.0f}%", ha="center",
                    va="center", fontsize=10)
    ax.set_xlabel("injected body size (=> pinch delay)")
    ax.set_ylabel("injected amplitude [x impulse threshold]")
    ax.set_title(f"Entry-doublet detector recovery over {len(paths)} real "
                 f"MARS slices\ncontrols: {n_fp} false positives in "
                 f"{len(ctl)} clean runs")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig("figures/fig_injection_doublet.png", dpi=150)
    print("wrote figures/fig_injection_doublet.png")


if __name__ == "__main__":
    main()
