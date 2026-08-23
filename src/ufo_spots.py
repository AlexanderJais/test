#!/usr/bin/env python3
"""Use a geocoded UFO event dataset to rank 'good' hydrophone spots (and times).

Data: Larry Hatch's *U database (~18k geocoded records), from
richgel999/ufo_data (bin/hatch_udb.json). Each record has LatLong + date +
an attribute set that includes 'SND: UFO sounds heard or recorded'.

Method: for each public, accessible hydrophone (location + active window),
count Hatch events within 25/50/100 km -> spatial 'good spots'. Then test
temporal overlap with the sensor's recording window -> 'good times'.

HARD FINDING: Hatch's data ends in 2003; every accessible public hydrophone
starts 2007+ (most 2014+). So this dataset yields good SPOTS (where to watch)
but NOT good TIMES (no event falls inside any sensor's window). For times, a
current dataset (e.g. NUFORC, which runs to the present) must be intersected
with a continuous sensor's coverage (MARS: 256 kHz, 2015-present).
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/user/richgel999/ufo_data/bin/hatch_udb.json"

# public hydrophones: name, lat, lon, active window
HYDRO = [
    ("MARS Monterey Bay", 36.713, -122.186, "2015-present, 256kHz continuous"),
    ("SanctSound Monterey Bay", 36.80, -121.98, "2018-2021"),
    ("SanctSound Channel Is.", 34.04, -120.08, "2018-2021"),
    ("SanctSound Stellwagen", 42.4, -70.4, "2018-2021"),
    ("SanctSound Florida Keys", 24.6, -81.9, "2018-2021"),
    ("SanctSound Grays Reef", 31.4, -80.9, "2018-2021"),
    ("SanctSound Olympic Coast", 47.9, -124.9, "2018-2021"),
    ("SanctSound Hawaii", 20.7, -157.0, "2018-2021"),
    ("NRS05 SoCal", 33.9, -119.58, "2014-2015"),
    ("NRS08 Mid-Atlantic", 39.01, -67.27, "2016-2018"),
    ("OOI Oregon (cabled)", 44.55, -125.2, "2015-present"),
    ("HARP SOCAL H", 32.942, -119.170, "2006-2012"),
    ("HARP SOCAL N", 32.370, -118.563, "2006-2012"),
]


def load(path):
    d = json.load(open(path, encoding="utf-8-sig"))
    recs = d if isinstance(d, list) else next(v for v in d.values() if isinstance(v, list))
    lat, lon, snd = [], [], []
    for r in recs:
        ll = r.get("key_vals", {}).get("LatLong")
        if not ll:
            continue
        try:
            a, b = map(float, ll.split()[:2])
        except ValueError:
            continue
        lat.append(a); lon.append(b)
        snd.append(any("SND" in x for x in r.get("attributes", [])))
    return np.array(lat), np.array(lon), np.array(snd)


def hav(la, lo, la2, lo2):
    la, lo, la2, lo2 = map(np.radians, [la, lo, la2, lo2])
    return 2 * 6371 * np.arcsin(np.sqrt(
        np.sin((la2 - la) / 2) ** 2 + np.cos(la) * np.cos(la2) * np.sin((lo2 - lo) / 2) ** 2))


def main():
    lat, lon, snd = load(DATA)
    print(f"{len(lat)} geocoded events\n")
    print(f"{'hydrophone':26s} {'<25':>4}{'<50':>5}{'<100':>6}{'SND<50':>7}  window")
    rows = []
    for name, hla, hlo, note in HYDRO:
        dk = hav(lat, lon, hla, hlo)
        rows.append((int((dk < 50).sum()), name, int((dk < 25).sum()),
                     int((dk < 100).sum()), int(((dk < 50) & snd).sum()), note))
    for n50, name, n25, n100, nsnd, note in sorted(rows, reverse=True):
        print(f"{name:26s} {n25:4d}{n50:5d}{n100:6d}{nsnd:7d}  {note}")

    fig, ax = plt.subplots(figsize=(13, 8))
    m = (lon > -130) & (lon < -65) & (lat > 23) & (lat < 50)
    ax.scatter(lon[m & ~snd], lat[m & ~snd], s=3, c="0.6", alpha=0.35,
               label="Hatch UFO event (pre-2003)")
    ax.scatter(lon[m & snd], lat[m & snd], s=12, c="tab:orange", alpha=0.8,
               label="...with SOUND reported")
    for name, hla, hlo, note in HYDRO:
        n = int((hav(lat, lon, hla, hlo) < 50).sum())
        ax.scatter(hlo, hla, marker="*", s=120 + n * 12, c="tab:blue",
                   edgecolor="k", zorder=5)
        ax.annotate(f"{name} ({n})", (hlo, hla), fontsize=7.5, xytext=(4, 4),
                    textcoords="offset points", fontweight="bold")
    ax.set_xlim(-130, -65); ax.set_ylim(23, 50)
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ax.set_title("'Good spots': public hydrophones (star, size = UFO events within 50 km)\n"
                 "over Hatch UFO events. Events end 2003; sensors start 2007+ (no time overlap)")
    ax.legend(loc="lower left"); ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig("figures/fig_ufo_hydrophone_spots.png", dpi=115)
    print("\nwrote figures/fig_ufo_hydrophone_spots.png")


if __name__ == "__main__":
    main()
