#!/usr/bin/env python3
"""Injection calibration for the Doppler-drift tracker (M0 gate).

Synthesizes a moving narrowband tone -- frequency following the CPA
S-curve, amplitude following the range hyperbola RL(t) = -20 log10(R(t)) --
injects it into real MARS slices in a quiet sub-band (28-44 kHz, above most
biosonar, avoiding the 50 kHz instrument line), and measures:
  * recovery: fraction of injections whose fitted speed is within 25% of
    truth, vs injected SNR and true speed;
  * false transits: how often UNMODIFIED ocean yields a spurious S-curve
    fit flagged is_transit (drifting ship lines / whistles).

Usage: python3 src/injection_doppler.py [--slices data/slices] [--workers 3]
"""

import argparse
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spectral_analysis import load_pressure
import doppler_track as dt

BAND = (28000.0, 44000.0)
F0 = 36000.0
SPEEDS_KN = [50.0, 150.0, 400.0]
SNRS_DB = [0.0, 6.0, 12.0, 20.0]   # in-band tone SNR vs local noise
R_CPA = 700.0
RNG_SEED = 20260822


def band_noise_rms(p, fs):
    xb = dt.signal.sosfiltfilt(
        dt.signal.butter(4, BAND, btype="bandpass", fs=fs, output="sos"), p)
    return np.sqrt(np.mean(xb ** 2))


def synth_transit(fs, dur, v_ms, snr_db, noise_rms, rng):
    t = np.arange(int(dur * fs)) / fs
    t_c = dur / 2
    f_inst = dt.scurve(t, F0, v_ms, R_CPA, t_c)
    phase = 2 * np.pi * np.cumsum(f_inst) / fs
    R = np.sqrt(R_CPA ** 2 + (v_ms * (t - t_c)) ** 2)
    amp = (R_CPA / R)              # spherical-spreading envelope (norm at CPA)
    tone = amp * np.sin(phase)
    tone /= np.sqrt(np.mean(tone ** 2)) + 1e-30
    return tone * noise_rms * 10 ** (snr_db / 20)


def run_slice(path):
    rng = np.random.default_rng(RNG_SEED)
    name = os.path.basename(path).split("_first")[0]
    p_upa, fs = load_pressure(path)
    p = p_upa * 1e-6
    nrms = band_noise_rms(p, fs)
    rows = []
    # false-transit control on unmodified ocean
    fit0 = dt.detect(p, fs, BAND)
    rows.append({"slice": name, "control": True,
                 "false_transit": bool(fit0 and fit0["is_transit"]),
                 "false_v_kn": (fit0 or {}).get("v_kn")})
    dur = len(p) / fs
    for v_kn in SPEEDS_KN:
        v = v_kn * dt.KN
        for snr in SNRS_DB:
            pi = p + synth_transit(fs, dur, v, snr, nrms, rng)
            fit = dt.detect(pi, fs, BAND)
            ok = bool(fit and fit["is_transit"]
                      and abs(fit["v_kn"] - v_kn) / v_kn < 0.25)
            rows.append({"slice": name, "v_kn": v_kn, "snr_db": snr,
                         "recovered": ok,
                         "fit_v_kn": (fit or {}).get("v_kn")})
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
    json.dump(rows, open("results/injection_doppler.json", "w"), indent=1)

    ctl = [r for r in rows if r.get("control")]
    inj = [r for r in rows if not r.get("control")]
    n_ft = sum(1 for r in ctl if r["false_transit"])
    print(f"false transits on unmodified ocean: {n_ft}/{len(ctl)} slices")
    mat = np.zeros((len(SNRS_DB), len(SPEEDS_KN)))
    for i, snr in enumerate(SNRS_DB):
        for j, v in enumerate(SPEEDS_KN):
            sel = [r for r in inj if r["snr_db"] == snr and r["v_kn"] == v]
            mat[i, j] = np.mean([r["recovered"] for r in sel]) if sel else 0
    for i, snr in enumerate(SNRS_DB):
        print(f"  SNR={snr:+.0f}dB: " + "  ".join(
            f"{v}kn->{mat[i,j]*100:.0f}%" for j, v in enumerate(SPEEDS_KN)))

    fig, ax = plt.subplots(figsize=(6.5, 5))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(SPEEDS_KN)))
    ax.set_xticklabels([f"{v} kn" for v in SPEEDS_KN])
    ax.set_yticks(range(len(SNRS_DB)))
    ax.set_yticklabels([f"{s:+.0f} dB" for s in SNRS_DB])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i,j]*100:.0f}%", ha="center", va="center")
    ax.set_xlabel("true transit speed")
    ax.set_ylabel("injected tone SNR (in 28-44 kHz band)")
    ax.set_title(f"Doppler-tracker speed recovery (within 25%) over "
                 f"{len(paths)} MARS slices\nfalse transits on clean ocean: "
                 f"{n_ft}/{len(ctl)}")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.tight_layout()
    fig.savefig("figures/fig_injection_doppler.png", dpi=150)
    print("wrote figures/fig_injection_doppler.png")


if __name__ == "__main__":
    main()
