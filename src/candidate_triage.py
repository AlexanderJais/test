#!/usr/bin/env python3
"""Triage the discharge-template candidate chains from the archive
re-screen into a persistent-source catalog and a residual watchlist.

A plasma-propelled CRAFT is, by definition, transient: it passes through.
A candidate chain that recurs across many days at the SAME acoustic
fingerprint (repetition rate, spectral centroid, bubble-echo lag) is
therefore a fixed installation — here, the Monterey Bay 38 kHz
echosounder — not a craft. This is standard passive-sonar practice:
maintain a catalog of known persistent contacts and reject them, so the
residual is the small set of genuinely novel, one-off events that merit
array follow-up.

Clusters DISCHARGE-LIKE chains from batch_screen_*_v{2,3}.jsonl by
(log-rate, centroid, echo-lag); any cluster present on >= PERSIST_DAYS
distinct days is catalogued as a fixed source. Everything else is the
residual watchlist, ranked by how machine-like and how isolated it is.

Usage:
    python3 src/candidate_triage.py [--tag _v3]
"""

import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PERSIST_DAYS = 8          # a fingerprint on >= 8 distinct days is a fixture
# cluster tolerances
TOL_LOGRATE = 0.05        # PRF is very stable
TOL_CEN_KHZ = 9.0         # centroid drifts with SNR/propagation
TOL_ECHO_MS = 0.10        # pulse/bubble structure is the fixed tell


def load_candidates(tag):
    rows = []
    for path in sorted(glob.glob(f"results/batch_screen_*{tag}.jsonl")):
        for l in open(path):
            r = json.loads(l)
            if r.get("status") != "ok":
                continue
            for t in r.get("chains", []) + r.get("trains", []):
                if t["class"].startswith("DISCHARGE"):
                    rows.append({
                        "date": r["date"], "n": t["n"],
                        "rate": t["rate_hz"], "cv": t["ipi_cv"],
                        "cen": t["med_centroid_hz"],
                        "echo": t["med_echo_lag_ms"],
                        "fbw": t["med_frac_bw"],
                        "hf": t.get("med_hf_frac", 0.0)})
    return rows


def fingerprint(r):
    return (np.log10(max(r["rate"], 1e-3)), r["cen"] / 1000.0, r["echo"])


def cluster(rows):
    """Single-link clustering in fingerprint space."""
    clusters = []
    for r in sorted(rows, key=lambda x: x["rate"]):
        fp = fingerprint(r)
        placed = False
        for c in clusters:
            cf = c["centroid_fp"]
            if (abs(fp[0] - cf[0]) < TOL_LOGRATE
                    and abs(fp[1] - cf[1]) < TOL_CEN_KHZ
                    and abs(fp[2] - cf[2]) < TOL_ECHO_MS):
                c["members"].append(r)
                m = c["members"]
                c["centroid_fp"] = tuple(np.mean(
                    [fingerprint(x) for x in m], axis=0))
                placed = True
                break
        if not placed:
            clusters.append({"members": [r], "centroid_fp": fp})
    for c in clusters:
        c["days"] = sorted({m["date"] for m in c["members"]})
        c["n_days"] = len(c["days"])
        c["rate"] = float(np.median([m["rate"] for m in c["members"]]))
        c["cen_khz"] = float(np.median([m["cen"] for m in c["members"]])
                             / 1000.0)
        c["echo_ms"] = float(np.median([m["echo"] for m in c["members"]]))
        c["min_cv"] = float(np.min([m["cv"] for m in c["members"]]))
    return clusters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="_v3")
    ap.add_argument("--out", default="results")
    ap.add_argument("--figdir", default="figures")
    args = ap.parse_args()

    rows = load_candidates(args.tag)
    clusters = cluster(rows)
    catalog = sorted([c for c in clusters if c["n_days"] >= PERSIST_DAYS],
                     key=lambda c: -c["n_days"])
    residual = sorted([c for c in clusters if c["n_days"] < PERSIST_DAYS],
                      key=lambda c: (c["min_cv"], -len(c["members"])))

    n_cat_chains = sum(len(c["members"]) for c in catalog)
    n_res_chains = sum(len(c["members"]) for c in residual)
    print(f"{len(rows)} DISCHARGE chains -> {len(clusters)} fingerprints")
    print(f"CATALOGUED fixed sources (>= {PERSIST_DAYS} days): "
          f"{len(catalog)} fingerprints, {n_cat_chains} chains")
    for c in catalog:
        print(f"  {c['n_days']:3d} days  rate={c['rate']:.2f}Hz  "
              f"cen={c['cen_khz']:.1f}kHz  echo={c['echo_ms']:.2f}ms  "
              f"(min cv {c['min_cv']:.3f})")
    print(f"RESIDUAL watchlist (< {PERSIST_DAYS} days): {len(residual)} "
          f"fingerprints, {n_res_chains} chains on "
          f"{len(set(m['date'] for c in residual for m in c['members']))} "
          f"days")
    print("  most machine-like residuals (lowest CV, isolated):")
    for c in residual[:15]:
        d = c["members"][0]
        print(f"    {c['days'][0]}  rate={c['rate']:.2f}Hz  "
              f"cen={c['cen_khz']:.1f}kHz  echo={c['echo_ms']:.2f}ms  "
              f"cv={c['min_cv']:.3f}  hf={d['hf']}  n_days={c['n_days']}")

    out = {"tag": args.tag, "persist_days": PERSIST_DAYS,
           "n_chains": len(rows), "n_fingerprints": len(clusters),
           "catalog": [{k: c[k] for k in ("n_days", "rate", "cen_khz",
                                          "echo_ms", "min_cv", "days")}
                       for c in catalog],
           "residual": [{"date": c["days"][0], "n_days": c["n_days"],
                         "rate": c["rate"], "cen_khz": c["cen_khz"],
                         "echo_ms": c["echo_ms"], "min_cv": c["min_cv"],
                         "n_chains": len(c["members"])}
                        for c in residual]}
    with open(os.path.join(args.out, "candidate_triage.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    # figure: fingerprint scatter, catalogued vs residual
    fig, ax = plt.subplots(figsize=(11, 6))
    for c in clusters:
        cat = c["n_days"] >= PERSIST_DAYS
        ax.scatter(c["rate"], c["cen_khz"],
                   s=20 + 6 * c["n_days"],
                   c="tab:red" if cat else "tab:green",
                   alpha=0.6, edgecolors="k", linewidths=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("repetition rate [Hz]")
    ax.set_ylabel("spectral centroid [kHz]")
    ax.set_title(f"Discharge candidates by fingerprint: "
                 f"{len(catalog)} catalogued fixed sources (red, "
                 f">={PERSIST_DAYS} days) vs {len(residual)} residual "
                 f"watchlist (green)", fontsize=11)
    ax.scatter([], [], c="tab:red", label="catalogued fixed installation")
    ax.scatter([], [], c="tab:green", label="residual (novel/transient)")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(args.figdir, "fig_triage.png"), dpi=150)
    print(f"wrote {args.figdir}/fig_triage.png")


if __name__ == "__main__":
    main()
