#!/usr/bin/env python3
"""Cross-slice catalog of narrowband lines + physics band-mismatch chart.

Clusters the per-slice line detections from ``spectral_analysis.py`` by
frequency, measures how many of the six survey slices each spectral
feature appears in (persistence), classifies the persistent ones, and
renders:

* ``figures/fig_psd_lines.png``  — all six calibrated PSDs overlaid with
  the persistent features annotated,
* ``figures/fig_band_mismatch.png`` — the hydrophone band next to real
  electron-plasma-frequency regimes on one log axis,
* ``results/line_catalog.md``    — the feature table used in REPORT.md.
"""

import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FINE_TOL = lambda f: max(3.0, 5e-4 * f)     # cluster sidebands of one tone
GROUP_REL = 5e-3                            # merge clusters into features


def load_lines(results_dir):
    per_slice = {}
    for path in sorted(glob.glob(os.path.join(results_dir, "*_lines.json"))):
        with open(path) as fh:
            d = json.load(fh)
        per_slice[d["slice"]] = d["lines"]
    return per_slice


def cluster(per_slice):
    """Two-stage clustering -> list of features."""
    rows = [(L["freq_hz"], L["snr_db"], s)
            for s, lines in per_slice.items() for L in lines]
    rows.sort()
    clusters = []
    for f, snr, s in rows:
        if clusters and f - clusters[-1]["fmax"] <= FINE_TOL(f):
            c = clusters[-1]
            c["fmax"] = f
            c["members"].append((f, snr, s))
        else:
            clusters.append({"fmax": f, "members": [(f, snr, s)]})
    feats = []
    for c in clusters:
        fc = float(np.average([m[0] for m in c["members"]],
                              weights=[m[1] for m in c["members"]]))
        if feats and abs(fc - feats[-1]["freq"]) / fc <= GROUP_REL:
            feats[-1]["members"] += c["members"]
            m = feats[-1]["members"]
            feats[-1]["freq"] = float(np.average([x[0] for x in m],
                                                 weights=[x[1] for x in m]))
        else:
            feats.append({"freq": fc, "members": c["members"]})
    n_slices = len(per_slice)
    for ft in feats:
        ft["slices"] = sorted({m[2] for m in ft["members"]})
        ft["persistence"] = f"{len(ft['slices'])}/{n_slices}"
        ft["max_snr"] = max(m[1] for m in ft["members"])
    return feats


def classify(ft, feats):
    f, n = ft["freq"], len(ft["slices"])
    if f < 100:
        # Lines here sit at different frequencies in every slice: broadband
        # shipping/flow noise resolved into wandering peaks, not one tone.
        return "low-frequency ambient / vessel band (frequency wanders)"
    persistent = [x["freq"] for x in feats
                  if len(x["slices"]) >= 5 and x["freq"] >= 1000]
    for base in persistent:
        for k in (2, 3):
            if base < f and abs(f - k * base) / f < 6e-3:
                return f"harmonic {k}x of {base/1e3:.2f} kHz system tone"
    if n >= 5:
        return "persistent system/instrument tone (fixed EMI or self-noise)"
    if f < 1000:
        return "transient ship tonal"
    # Within-slice harmonic stack (f, 2f, 3f...) -> tonal call or whistle.
    others = [x["freq"] for x in feats if x is not ft
              and set(x["slices"]) & set(ft["slices"])]
    for base in [x for x in others if 100 < x < f]:
        for k in (2, 3):
            if abs(f - k * base) / f < 6e-3:
                return (f"harmonic stack on {base/1e3:.2f} kHz fundamental "
                        "(biological call or vessel whistle)")
    return "intermittent engineering tone / echosounder-class"


def write_catalog(feats, out_md):
    lines = [
        "| feature (Hz) | persistence | max SNR (dB) | classification |",
        "|---:|:---:|---:|---|",
    ]
    for ft in feats:
        lines.append(f"| {ft['freq']:,.1f} | {ft['persistence']} | "
                     f"{ft['max_snr']:.1f} | {ft['class']} |")
    with open(out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def log_max_downsample(f, pxx, nbins=4000, fmin=10.0):
    """Max-pool PSD into log-spaced bins (preserves narrow peaks)."""
    edges = np.logspace(np.log10(fmin), np.log10(f[-1]), nbins + 1)
    idx = np.searchsorted(f, edges)
    fo, po = [], []
    for a, b in zip(idx[:-1], idx[1:]):
        if b > a:
            j = a + int(np.argmax(pxx[a:b]))
            fo.append(f[j]); po.append(pxx[j])
    return np.array(fo), np.array(po)


def fig_psd(results_dir, feats, figpath):
    fig, ax = plt.subplots(figsize=(13, 6))
    for path in sorted(glob.glob(os.path.join(results_dir, "*_psd.npz"))):
        z = np.load(path)
        fo, po = log_max_downsample(z["f"], z["pxx"])
        name = os.path.basename(path).replace("_psd.npz", "")
        ax.plot(fo, 10 * np.log10(po), lw=0.6, alpha=0.75,
                label=name.replace("MARS_", ""))
    for ft in feats:
        if len(ft["slices"]) >= 5 and ft["freq"] > 1000:
            ax.axvline(ft["freq"], color="k", ls=":", lw=0.8, alpha=0.6)
            ax.annotate(f"{ft['freq']/1e3:.1f} kHz", (ft["freq"], 92),
                        rotation=90, fontsize=7, ha="right", va="top")
    ax.set_xscale("log")
    ax.set_xlim(10, 1.28e5)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("PSD [dB re 1 µPa²/Hz]")
    ax.set_title("MARS hydrophone Welch PSDs (0.98 Hz resolution) — "
                 "persistent narrowband features marked")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(figpath, dpi=150)


def fig_band_mismatch(feats, figpath):
    """Hydrophone band vs electron-plasma-frequency regimes, one log axis."""
    fig, ax = plt.subplots(figsize=(13, 4.5))
    # Acoustic side
    ax.axvspan(10, 1.28e5, color="tab:blue", alpha=0.25,
               label="MARS hydrophone band (10 Hz – 128 kHz)")
    for ft in feats:
        if len(ft["slices"]) >= 5:
            ax.axvline(ft["freq"], ymin=0, ymax=0.35, color="tab:blue",
                       lw=1.0)
    ax.text(1.1e3, 0.06, "persistent lines found in the data\n"
            "(all instrument/system tones)", fontsize=8, color="tab:blue",
            transform=ax.get_xaxis_transform())
    # Plasma side: f_pe = 8980 * sqrt(n_e[cm^-3]) Hz
    regimes = [
        (2.0e4, "solar wind at 1 AU\n(nₑ≈5 cm⁻³) — E-field\noscillation in space,\nnot sound in water"),
        (2.8e6, "ionosphere F-layer\n(nₑ=1e5 cm⁻³)"),
        (9.0e6, "ionosphere peak\n(nₑ=1e6 cm⁻³)"),
        (9.0e10, "tokamak core\n(nₑ=1e14 cm⁻³)"),
        (2.8e15, "solid-density laser\nplasma (nₑ=1e23 cm⁻³)"),
        (5.1e15, "water's bound-electron\nplasmon (~21 eV, UV)"),
    ]
    for i, (f, label) in enumerate(regimes):
        ax.axvline(f, ymin=0.45, ymax=0.72, color="tab:red", lw=2)
        ax.text(f, 0.75 + 0.11 * (i % 2), label, fontsize=7.5,
                color="tab:red", ha="center",
                transform=ax.get_xaxis_transform())
    ax.annotate("", xy=(2.8e6, 0.42), xytext=(1.28e5, 0.42),
                xycoords=ax.get_xaxis_transform(),
                arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text(6.5e5, 0.33, "22–70× gap to the nearest\nreal plasma environment",
            fontsize=8, ha="center", transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_xlim(1, 1e17)
    ax.set_yticks([])
    ax.set_xlabel("frequency [Hz]  —  fₚₑ = 8.98 kHz × √(nₑ [cm⁻³])")
    ax.set_title("Why no electron plasma oscillation can appear in hydrophone "
                 "data: frequency scales (and the acoustic-vs-electrostatic "
                 "field mismatch)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(figpath, dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()

    per_slice = load_lines(args.results)
    feats = cluster(per_slice)
    for ft in feats:
        ft["class"] = classify(ft, feats)
    write_catalog(feats, os.path.join(args.results, "line_catalog.md"))
    print(open(os.path.join(args.results, "line_catalog.md")).read())
    fig_psd(args.results, feats, os.path.join(args.figdir,
                                              "fig_psd_lines.png"))
    fig_band_mismatch(feats, os.path.join(args.figdir,
                                          "fig_band_mismatch.png"))
    print(f"wrote {args.figdir}/fig_psd_lines.png and fig_band_mismatch.png")


if __name__ == "__main__":
    main()
