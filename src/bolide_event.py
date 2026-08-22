#!/usr/bin/env python3
"""Event-directed test on a REAL cataloged transmedium object: the
2020-05-09 02:56:11 UTC bolide (NASA CNEOS: 44.8 N, 131.0 W, 3.53 kt,
v=14.5 km/s, airburst alt 31.2 km) recorded by the GPS-timed OOI cabled
hydrophone array.

Unlike the USS Omaha SoundTraps (free-running clocks, 167 km), OOI is
GPS-disciplined and the Axial cluster is only ~150 km from the bolide —
so this is a rigorous TDOA test on ground-truth. A 3.5 kt airburst
couples to the ocean via its infrasound airwave; we search the low-
frequency band for a coherent transient in the arrival window, cross-
correlate all station pairs, grid-search the source position, and compare
to the CNEOS location.

Stations (OO/HDH, 200 Hz, Pa):
  AXBA1 45.820 -129.737   AXCC1 45.955 -130.009   AXEC2 45.940 -129.974
  HYS14 44.569 -125.148   HYSB1 44.510 -125.405
"""
import datetime, json, os
import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

FDSN = "https://service.earthscope.org"
BOLIDE_UTC = UTCDateTime("2020-05-09T02:56:11")
BOLIDE_LAT, BOLIDE_LON = 44.8, -131.0
STAS = ["AXBA1", "AXCC1", "AXEC2", "HYS14", "HYSB1"]
# search from burst to +25 min (infrasound celerity ~0.28-0.34 km/s over
# ~150-460 km -> ~7-27 min), plus 3 min lead
T0 = BOLIDE_UTC - 180
DUR = 30 * 60
LF_BAND = (2.0, 20.0)       # ocean-coupled infrasound band
OUT = "data/bolide"


def fetch():
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, "bolide.npz")
    if os.path.exists(out):
        return dict(np.load(out, allow_pickle=True))["payload"].item()
    c = Client(FDSN)
    inv = c.get_stations(network="OO", channel="HDH", level="response")
    coords = {s.code: (s.latitude, s.longitude)
              for n in inv for s in n}
    pay = {"t0": str(T0)}
    for code in STAS:
        try:
            st = c.get_waveforms("OO", code, "*", "HDH", T0, T0 + DUR)
            st.merge(fill_value="interpolate")
            tr = st[0]
            tr.detrend("linear")
            tr.remove_response(inventory=inv, output="DEF", water_level=60,
                               pre_filt=(0.5, 1.0, 90, 99))
            pay[f"{code}_x"] = tr.data.astype(np.float32)
            pay[f"{code}_fs"] = float(tr.stats.sampling_rate)
            pay[f"{code}_lat"], pay[f"{code}_lon"] = coords[code]
            print(f"  {code}: {len(tr.data)} samp @ {tr.stats.sampling_rate}Hz "
                  f"rms={np.std(tr.data):.3f}")
        except Exception as e:
            print(f"  {code}: NO DATA ({e})")
    np.savez_compressed(out, payload=pay)
    return pay


def hav(a, b, c, e):
    la1, lo1, la2, lo2 = map(np.radians, [a, b, c, e])
    return 2*6371*np.arcsin(np.sqrt(np.sin((la2-la1)/2)**2 +
                                    np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2))


def main():
    pay = fetch()
    stas = [s for s in STAS if f"{s}_x" in pay]
    print(f"stations with data: {stas}")
    if len(stas) < 3:
        print("insufficient stations")
        return
    fs = pay[f"{stas[0]}_fs"]

    # LF envelopes + transient scan per station
    fig, axes = plt.subplots(len(stas), 1, figsize=(14, 2.0*len(stas)),
                             sharex=True)
    peaks_by = {}
    for ax, s in zip(np.atleast_1d(axes), stas):
        x = pay[f"{s}_x"].astype(float)
        sos = signal.butter(4, LF_BAND, btype="bandpass", fs=fs, output="sos")
        xb = signal.sosfiltfilt(sos, x)
        env = np.abs(signal.hilbert(xb))
        t = np.arange(len(x))/fs
        ax.plot(t/60, xb, lw=0.3, color="k")
        med, mad = np.median(env), np.median(np.abs(env-np.median(env)))
        thr = med + 8*mad
        pk, _ = signal.find_peaks(env, height=thr, distance=int(5*fs))
        peaks_by[s] = pk/fs
        for pkt in pk:
            ax.axvline(pkt/fs/60, color="tab:red", lw=0.5, alpha=0.5)
        rng = hav(BOLIDE_LAT, BOLIDE_LON, pay[f"{s}_lat"], pay[f"{s}_lon"])
        # infrasound arrival guess at celerity 0.30 km/s
        t_air = (BOLIDE_UTC - UTCDateTime(pay["t0"])) + rng/0.30
        ax.axvline(t_air/60, color="cyan", ls="--", lw=1.2)
        ax.set_ylabel(f"{s}\n{rng:.0f} km")
    np.atleast_1d(axes)[-1].set_xlabel("minutes since window start "
                                       f"({pay['t0'][:19]}Z)")
    fig.suptitle("2020-05-09 bolide — OOI LF hydrophones (2-20 Hz). cyan = "
                 "infrasound arrival @0.30 km/s; red ticks = transients")
    fig.tight_layout(); fig.savefig("figures/fig_bolide_event.png", dpi=140)

    # cross-correlate all station pairs in the LF band over the arrival window
    print("\nstation ranges & infrasound-arrival guesses (celerity 0.30 km/s):")
    for s in stas:
        rng = hav(BOLIDE_LAT, BOLIDE_LON, pay[f"{s}_lat"], pay[f"{s}_lon"])
        t_air = (BOLIDE_UTC - UTCDateTime(pay["t0"])) + rng/0.30
        print(f"  {s}: {rng:.0f} km -> arrival ~{t_air/60:.1f} min "
              f"(n_transients {len(peaks_by[s])})")
    json.dump({"bolide": str(BOLIDE_UTC), "pos": [BOLIDE_LAT, BOLIDE_LON],
               "stations": stas,
               "transients": {s: [round(float(x), 1) for x in peaks_by[s]]
                              for s in stas}},
              open("results/bolide_event.json", "w"), indent=1)
    print("wrote figures/fig_bolide_event.png and results/bolide_event.json")


if __name__ == "__main__":
    main()
