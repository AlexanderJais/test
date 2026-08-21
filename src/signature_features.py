#!/usr/bin/env python3
"""Propulsion-signature screening for the array windows.

The acoustic discriminants that would separate a plasma-propelled craft
from everything else the ocean contains:

* conventional propeller craft — cavitation broadband whose envelope is
  modulated at shaft/blade rate: DEMON analysis shows discrete modulation
  lines; machinery adds stable narrowband tonals.
* MHD "caterpillar" drive (Yamato-1 class) — NO blade-rate modulation,
  but electrolysis bubble hiss (sustained broadband) and power-converter
  harmonic combs.
* supercavitating vehicle — extreme sustained broadband, no blade lines.
* plasma sheath / discharge propulsion — sustained broadband PLUS high
  impulsiveness (kHz-rate micro-discharge shocks -> heavy-tailed,
  high-kurtosis amplitude statistics), no blade lines.
* T-phase (geophysical) — a single emergent packet, not sustained.
* fin/blue whale calls (biologic) — regular ~1 s pulses at 15-25 Hz with
  stereotyped inter-pulse intervals.

This script computes, for every 30 s window of every station in the
fetched OOI array data: broadband level, kurtosis (impulsiveness), DEMON
modulation-line SNR, narrowband tonal count, and pulse-train regularity —
then applies rule-based classification, including the anomalous-
propulsion criterion (sustained broadband + high kurtosis + no blade
lines + no pulse regularity). Outputs a feature table, class counts, and
a feature-space figure with the hypothetical plasma-craft regions marked.

Usage:
    python3 src/signature_features.py [--data data/array] [--out results]
"""

import argparse
import glob
import json
import os

import numpy as np
from scipy import signal, stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WIN_S, HOP_S = 30.0, 15.0
BB_BAND = (10.0, 90.0)          # broadband metric band
DEMON_BAND = (30.0, 90.0)       # carrier band for envelope demodulation
DEMON_MOD = (0.5, 15.0)         # shaft/blade modulation search range
CALL_BAND = (14.0, 28.0)        # fin/blue whale pulse band


def load_windows(path):
    z = np.load(path)
    stations = sorted({k.split("_")[0] for k in z.files if k.endswith("_x")})
    return {s: (z[f"{s}_x"].astype(float), float(z[f"{s}_fs"]))
            for s in stations}


def bandpass(x, fs, band):
    sos = signal.butter(4, band, btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x)


def demon_snr(xw, fs):
    """Peak SNR (dB) of the envelope-modulation (DEMON) spectrum."""
    env = np.abs(signal.hilbert(bandpass(xw, fs, DEMON_BAND)))
    env = env - env.mean()
    f, p = signal.welch(env, fs=fs, nperseg=1024)
    sel = (f >= DEMON_MOD[0]) & (f <= DEMON_MOD[1])
    fm, pm = f[sel], p[sel]
    base = signal.medfilt(pm, 21)
    snr = 10 * np.log10(pm / np.maximum(base, 1e-30))
    k = int(np.argmax(snr))
    return float(snr[k]), float(fm[k])


def tonal_count(xw, fs):
    f, p = signal.welch(xw, fs=fs, nperseg=2048)
    sel = f >= 5.0
    base = signal.medfilt(p[sel], 101)
    snr = 10 * np.log10(p[sel] / np.maximum(base, 1e-30))
    peaks, _ = signal.find_peaks(snr, height=8.0, distance=5)
    return int(len(peaks))


def spike_count(xw, fs):
    """Number of distinct broadband impulse excursions (>8 MAD) in the
    window. A discharge-train propulsion signature needs MANY per second;
    a data glitch or single close transient gives only a few."""
    xb = bandpass(xw, fs, BB_BAND)
    env = np.abs(signal.hilbert(xb))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    peaks, _ = signal.find_peaks(env, height=med + 8 * mad,
                                 distance=int(0.05 * fs))
    return int(len(peaks))


def pulse_regularity(xw, fs):
    """Detect regular biologic pulse trains (e.g. fin whale 20 Hz calls):
    returns (n_pulses, coefficient of variation of inter-pulse intervals).
    """
    env = np.abs(signal.hilbert(bandpass(xw, fs, CALL_BAND)))
    sos = signal.butter(4, 2.0, btype="low", fs=fs, output="sos")
    env = signal.sosfiltfilt(sos, env)
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    peaks, _ = signal.find_peaks(env, height=med + 6 * mad,
                                 distance=int(3 * fs))
    if len(peaks) < 3:
        return int(len(peaks)), np.nan
    ipi = np.diff(peaks) / fs
    return int(len(peaks)), float(np.std(ipi) / np.mean(ipi))


def classify(feat, bb_baseline_db):
    bb_up = feat["bb_db"] - bb_baseline_db
    if feat["n_pulses"] >= 3 and feat["ipi_cv"] < 0.25:
        return "biologic pulse train (fin/blue whale class)"
    if feat["kurtosis"] > 50 and feat["n_spikes"] < 10:
        return "isolated spike (data glitch / single close transient)"
    if bb_up > 15:
        return "geophysical energy burst (T-phase class)"
    if bb_up > 6 and feat["demon_snr"] >= 10:
        return "conventional propeller craft (DEMON blade lines)"
    if bb_up > 6 and feat["n_tonals"] >= 3:
        return "vessel machinery (tonal comb)"
    if bb_up > 6 and feat["kurtosis"] > 5 and feat["n_spikes"] >= 30 \
            and feat["demon_snr"] < 10:
        return "ANOMALOUS-PROPULSION CANDIDATE (sustained impulse train, " \
               "no blade lines)"
    if bb_up > 6:
        return "broadband elevation, unclassified (weather/distant ship)"
    return "ambient"


def context_pass(rows, base):
    """T-phase packets are minutes long: any elevated window within 150 s
    of a station's dominant (>15 dB) energy burst in the same file
    belongs to that burst's packet/coda, not to a sustained source."""
    for wname in {r["window"] for r in rows}:
        for sta in {r["station"] for r in rows}:
            grp = [r for r in rows
                   if r["window"] == wname and r["station"] == sta]
            if not grp:
                continue
            peaks = [r for r in grp
                     if r["bb_db"] - base[sta] > 15]
            if not peaks:
                continue
            tpk = max(peaks, key=lambda r: r["bb_db"])["t_s"]
            for r in grp:
                if abs(r["t_s"] - tpk) <= 150 and \
                        r["bb_db"] - base[sta] > 6 and \
                        not r["class"].startswith("biologic"):
                    r["class"] = "geophysical energy burst (T-phase class)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/array")
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()

    rows = []
    for path in sorted(glob.glob(os.path.join(args.data, "*.npz"))):
        wname = os.path.splitext(os.path.basename(path))[0]
        for sta, (x, fs) in load_windows(path).items():
            n, h = int(WIN_S * fs), int(HOP_S * fs)
            for i0 in range(0, len(x) - n, h):
                xw = x[i0:i0 + n]
                xb = bandpass(xw, fs, BB_BAND)
                dsnr, dfreq = demon_snr(xw, fs)
                npul, cv = pulse_regularity(xw, fs)
                rows.append({
                    "window": wname, "station": sta,
                    "t_s": round(i0 / fs, 1),
                    "bb_db": round(float(
                        20 * np.log10(np.std(xb) / 1e-6)), 1),
                    "kurtosis": round(float(stats.kurtosis(
                        xb, fisher=False)), 2),
                    "demon_snr": round(dsnr, 1),
                    "demon_freq": round(dfreq, 2),
                    "n_tonals": tonal_count(xw, fs),
                    "n_spikes": spike_count(xw, fs),
                    "n_pulses": npul,
                    "ipi_cv": round(cv, 3) if np.isfinite(cv) else None,
                })

    # per-station ambient baseline = 20th percentile of broadband level
    base = {}
    for sta in {r["station"] for r in rows}:
        base[sta] = float(np.percentile(
            [r["bb_db"] for r in rows if r["station"] == sta], 20))
    for r in rows:
        r["ipi_cv"] = r["ipi_cv"] if r["ipi_cv"] is not None else np.nan
        r["class"] = classify(r, base[r["station"]])
        r["ipi_cv"] = None if not np.isfinite(r["ipi_cv"]) else r["ipi_cv"]
    context_pass(rows, base)

    counts = {}
    for r in rows:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    print(f"{len(rows)} windows classified:")
    for c, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {c}")

    with open(os.path.join(args.out, "signature_windows.json"), "w") as fh:
        json.dump({"baseline_db": base, "counts": counts, "windows": rows},
                  fh, indent=1)

    # ---- feature-space figure -----------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    classes = sorted(counts)
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(classes):
        sel = [r for r in rows if r["class"] == c]
        up = [r["bb_db"] - base[r["station"]] for r in sel]
        ax1.scatter(up, [r["demon_snr"] for r in sel], s=14, alpha=0.7,
                    color=cmap(i % 10), label=f"{c} (n={len(sel)})")
        ax2.scatter(up, [r["kurtosis"] for r in sel], s=14, alpha=0.7,
                    color=cmap(i % 10))
    # hypothetical plasma-craft regions
    ax1.axhspan(10, 30, xmin=0.55, alpha=0.12, color="tab:green")
    ax1.text(23, 26, "conventional propeller craft\n(blade-rate DEMON lines)",
             fontsize=8, color="tab:green")
    ax1.axhspan(-2, 10, xmin=0.55, alpha=0.12, color="tab:red")
    ax1.text(23, 7.5, "MHD drive / supercavitation /\nplasma sheath: loud "
             "but NO\nblade lines -> lands here", fontsize=8, color="tab:red")
    ax2.axhspan(3.5, 60, xmin=0.55, alpha=0.12, color="tab:red")
    ax2.text(22, 30, "plasma-sheath discharge:\nsustained broadband +\n"
             "impulse train (high kurtosis)", fontsize=8, color="tab:red")
    ax1.set_xlabel("broadband elevation above ambient baseline [dB]")
    ax1.set_ylabel("DEMON modulation-line SNR [dB]")
    ax2.set_xlabel("broadband elevation above ambient baseline [dB]")
    ax2.set_ylabel("amplitude kurtosis (impulsiveness)")
    ax2.set_yscale("log")
    ax1.legend(fontsize=7, loc="upper left")
    fig.suptitle("Propulsion-signature feature space — OOI array windows "
                 "vs where plasma-propulsion classes would land",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_signatures.png"), dpi=150)
    print(f"wrote {args.figdir}/fig_signatures.png")


if __name__ == "__main__":
    main()
