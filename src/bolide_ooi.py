#!/usr/bin/env python3
"""A meteor hitting the water: 2020-05-09 ocean bolide on OOI Axial hydrophones.

The real transmedium target (unlike an earthquake): an object from the sky
depositing energy into the ocean. CNEOS logged a bolide airburst over open
ocean at 2020-05-09 02:56:11 UTC, 44.8 N 131.0 W, radiated 3.5e10 J
(~0.12 kt). The three OOI Axial Seamount cabled hydrophones (AXBA1/AXCC1/
AXEC2, HDH, 200 Hz) sit ~150 km away and were recording.

Coupling geometry (why the arrival time is a window, not a point): the
airburst shock descends through the atmosphere and couples into the sea, then
runs in-water to the hydrophone. Two bounding paths:
  * in-water-coupled: ~(burst_alt)/0.34 + 150 km/1.48  -> ~150-200 s
  * atmospheric to overhead then down: ~150 km/0.34     -> ~440 s
So we search airburst+60 .. +540 s and, crucially, require the transient to
appear COHERENTLY on all three stations (a real distant source does; local
noise does not). That coherence is the detection, not a single envelope bump.
"""
import math
import datetime
import numpy as np
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FDSN = "https://service.earthscope.org"
BURST = UTCDateTime("2020-05-09T02:56:11")
EV = (44.8, -131.0)
STATIONS = ["AXBA1", "AXCC1", "AXEC2"]
T1, T2 = BURST - 120, BURST + 600     # 02:54:11 .. 03:06:11


def hav(a, b, c, e):
    a, b, c, e = map(math.radians, [a, b, c, e])
    return 2 * 6371 * math.asin(math.sqrt(
        math.sin((c - a) / 2) ** 2 + math.cos(a) * math.cos(c) * math.sin((e - b) / 2) ** 2))


def main():
    c = Client(FDSN, timeout=120)
    inv = c.get_stations(network="OO", channel="HDH", level="channel",
                         starttime=T1, endtime=T2)
    coord = {}
    for net in inv:
        for st in net:
            coord[st.code] = (st.latitude, st.longitude)

    traces = {}
    for s in STATIONS:
        try:
            stz = c.get_waveforms("OO", s, "*", "HDH", T1, T2)
            stz.merge(fill_value=0)
            tr = stz[0]
            tr.detrend("demean"); tr.detrend("linear")
            traces[s] = tr
            print(f"got {s}: {tr.stats.sampling_rate:.0f} Hz, {tr.stats.npts} samples")
        except Exception as e:
            print(f"no data {s}: {e}")
    if not traces:
        print("no waveforms"); return

    envs = {}
    fig, ax = plt.subplots(len(traces) + 1, 1, figsize=(13, 10), sharex=True)
    for i, (s, tr) in enumerate(traces.items()):
        fs = tr.stats.sampling_rate
        x = tr.data.astype(float)
        sos = signal.butter(4, (2, min(90, 0.49 * fs * 2 / 2)), btype="bandpass",
                            fs=fs, output="sos")
        xb = signal.sosfiltfilt(sos, x)
        env = np.abs(signal.hilbert(xb))
        sm = np.convolve(env, np.ones(int(2 * fs)) / int(2 * fs), "same")
        t = np.array([(tr.stats.starttime + k / fs) - BURST for k in range(len(x))])
        envs[s] = (t, sm)
        d = hav(*EV, *coord.get(s, EV))
        ax[i].plot(t, sm, lw=0.7)
        ax[i].set_ylabel(f"{s}\n({d:.0f} km)")
        ax[i].axvspan(150, 200, color="orange", alpha=0.12)
        ax[i].axvspan(430, 460, color="cyan", alpha=0.12)
        ax[i].axvline(0, color="red", ls="--", lw=1.5)

    # coherence: stack the normalized envelopes on a common grid
    tg = np.arange(-100, 560, 0.5)
    stack = []
    for s, (t, sm) in envs.items():
        smi = np.interp(tg, t, sm)
        z = (smi - np.median(smi)) / (np.median(np.abs(smi - np.median(smi))) + 1e-30) * 0.6745
        stack.append(z)
    stack = np.array(stack)
    coh = stack.min(axis=0)      # min across stations = present on ALL
    ax[-1].plot(tg, coh, color="k", lw=1)
    ax[-1].axvspan(150, 200, color="orange", alpha=0.12)
    ax[-1].axvspan(430, 460, color="cyan", alpha=0.12)
    ax[-1].axvline(0, color="red", ls="--", lw=1.5)
    ax[-1].set_ylabel("cross-stn\nmin-z")
    ax[-1].set_xlabel("seconds after airburst (red); orange=in-water band, cyan=atmospheric band")

    win = (tg >= 60) & (tg <= 540)
    ipk = np.where(win)[0][np.argmax(coh[win])]
    print(f"\nairburst 02:56:11 UTC; predicted in-water ~150-200s, atmospheric ~440s")
    print(f"strongest COHERENT (all-station) transient in +60..+540s: "
          f"{tg[ipk]:.0f}s after burst  min-z={coh[ipk]:.1f}  "
          f"(= {(BURST+float(tg[ipk])).strftime('%H:%M:%S')} UTC)")
    thr = 6
    print("coherent" if coh[ipk] >= thr else "no coherent all-station transient (min-z<6)"
          + " -> likely below detection at 0.12 kt / 150 km")

    ax[0].set_title("2020-05-09 ocean bolide (0.12 kt, 150 km) on OOI Axial HDH hydrophones")
    fig.tight_layout()
    fig.savefig("figures/fig_bolide_ooi.png", dpi=115)
    print("wrote figures/fig_bolide_ooi.png")


if __name__ == "__main__":
    main()
