#!/usr/bin/env python3
"""What could the CI01 snapping shrimp be reacting to at the UAP event?

The snap-rate rise at CI01 coincident with the 2019-07-16 (USS Omaha) event
is real and temperature-independent (see docs/CI01_TEMPERATURE_CORRECTION.md).
This asks the mechanistic follow-up the user posed: are the shrimp reacting
to an ENERGY input, and if so, in what channel?

Honest scope note, stated up front: a hydrophone measures acoustic PRESSURE.
It records an electromagnetic field ONLY if that field couples electrically
into the cable/preamp (e.g. mains hum). It cannot measure light, a static or
low-frequency magnetic field, water chemistry, or a pressure/current change
below its high-pass. So this can only test energy inputs the sensor can see:
  * acoustic bands (infrasound -> snap band),
  * mains-frequency electrical pickup (50/60 Hz) as an EM-coupling proxy.
If none of these leads the snap rise, the trigger -- if there is one -- was
in a channel this instrument does not capture, OR the rise is not a direct
reaction to a recorded stimulus.

For each candidate band we compute, over the event window (t0 +/- 90 min):
  1. band energy per 30 s (dB re the window median),
  2. the snap rate per 30 s,
  3. the lead-lag cross-correlation of each band against the snap rate
     (positive lag = band LEADS snaps = stimulus -> response),
  4. the band energy step across the event (mean over t0..t0+20 min minus
     t0-20..t0 min).
A genuine "the shrimp reacted to energy X" signature would be: band X shows a
coincident energy INCREASE that LEADS the snap-rate rise. The output says
whatever the data say.
"""
import datetime
import numpy as np
import soundfile as sf
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EVENT_FILE = "data/omaha/SanctSound_CI01_02_671883305_190716033643.flac"
FSTART = datetime.datetime(2019, 7, 16, 3, 36, 43)          # from filename
EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)         # Omaha encounter
CI01 = (34.0438, -120.0811)
EVENT_LATLON = (32.4894, -119.3647)
C_KM_S = 1.487
WIN_S = 30.0                 # analysis bin
SPAN_MIN = 90               # +/- around the event arrival
SNAP_BAND = (2000.0, 20000.0)
MAX_LAG_MIN = 10.0

# (label, kind, params). kind 'bp' = Butterworth band-pass; 'mains' = narrow
# band around a mains frequency used as an electrical-coupling proxy.
BANDS = [
    ("infrasound 1-20Hz", "bp", (1.0, 20.0)),
    ("low 20-100Hz",      "bp", (20.0, 100.0)),
    ("ship 100-500Hz",    "bp", (100.0, 500.0)),
    ("mid 0.5-2kHz",      "bp", (500.0, 2000.0)),
    ("SNAP 2-20kHz",      "bp", SNAP_BAND),
    ("EM 60Hz",           "mains", 60.0),
    ("EM 50Hz",           "mains", 50.0),
]


def haversine_km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, [a[0], a[1], b[0], b[1]])
    return 2 * 6371 * np.arcsin(np.sqrt(
        np.sin((la2 - la1) / 2) ** 2 +
        np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2))


def band_rms(x, fs, kind, p):
    """RMS of x restricted to a band. 'mains' uses a +/-1.5 Hz band-pass."""
    if kind == "bp":
        lo, hi = p
        hi = min(hi, 0.999 * fs / 2)
        sos = signal.butter(4, (lo, hi), btype="bandpass", fs=fs, output="sos")
    else:  # narrow band around a mains line
        f0 = p
        sos = signal.butter(4, (f0 - 1.5, f0 + 1.5), btype="bandpass",
                            fs=fs, output="sos")
    xb = signal.sosfiltfilt(sos, x)
    return np.sqrt(np.mean(xb ** 2)) + 1e-30


def snap_count(x, fs):
    sos = signal.butter(4, SNAP_BAND, btype="bandpass", fs=fs, output="sos")
    xb = signal.sosfiltfilt(sos, x)
    env = np.abs(signal.hilbert(xb))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    pk, _ = signal.find_peaks(env, height=med + 8 * mad, distance=int(0.01 * fs))
    return len(pk)


def main():
    info = sf.info(EVENT_FILE)
    fs = info.samplerate
    n = int(WIN_S * fs)

    range_km = haversine_km(CI01, EVENT_LATLON)
    t0_s = (EVENT_UTC - FSTART).total_seconds() + range_km / C_KM_S
    lo_s = max(0.0, t0_s - SPAN_MIN * 60)
    hi_s = min(info.frames / fs, t0_s + SPAN_MIN * 60)

    offs = np.arange(int(lo_s * fs), int(hi_s * fs) - n, n)
    tmin = (offs / fs) / 60.0                       # min since FSTART
    t0_min = t0_s / 60.0

    energy = {lbl: np.full(len(offs), np.nan) for lbl, _, _ in BANDS}
    snaprate = np.full(len(offs), np.nan)

    with sf.SoundFile(EVENT_FILE) as fh:
        for i, off in enumerate(offs):
            fh.seek(int(off))
            x = fh.read(n, dtype="float64")
            if x.ndim > 1:
                x = x[:, 0]
            for lbl, kind, p in BANDS:
                energy[lbl][i] = band_rms(x, fs, kind, p)
            snaprate[i] = snap_count(x, fs) * (60.0 / WIN_S)   # per minute

    # dB re window median
    edb = {}
    for lbl in energy:
        e = energy[lbl]
        med = np.nanmedian(e)
        edb[lbl] = 20 * np.log10(e / med) if med > 0 else np.full_like(e, np.nan)

    # lead-lag cross-correlation of each band vs snap rate
    lags = np.arange(-int(MAX_LAG_MIN * 60 / WIN_S),
                     int(MAX_LAG_MIN * 60 / WIN_S) + 1)
    sr = snaprate - np.nanmean(snaprate)
    print(f"event arrival t0 = {t0_min:.1f} min since {FSTART:%H:%M} UTC "
          f"(range {range_km:.0f} km); {len(offs)} x {WIN_S:.0f}s bins\n")
    print("lead-lag of each band vs SNAP rate "
          "(positive lag = band LEADS snaps = stimulus->response):")
    for lbl, _, _ in BANDS:
        if lbl.startswith("SNAP"):
            continue
        b = edb[lbl] - np.nanmean(edb[lbl])
        best_c, best_l = 0.0, 0
        for L in lags:
            if L >= 0:
                a, c = b[:len(b) - L], sr[L:]
            else:
                a, c = b[-L:], sr[:len(sr) + L]
            m = np.isfinite(a) & np.isfinite(c)
            if m.sum() < 10 or np.std(a[m]) == 0 or np.std(c[m]) == 0:
                continue
            cc = np.corrcoef(a[m], c[m])[0, 1]
            if abs(cc) > abs(best_c):
                best_c, best_l = cc, L
        lead = "band leads snaps" if best_l > 0 else "snaps lead band"
        print(f"  {lbl:18s}: peak |corr|={best_c:+.2f} at lag "
              f"{best_l * WIN_S / 60:+.1f} min ({lead})")

    # energy step across the event (+/-20 min)
    before = (tmin >= t0_min - 20) & (tmin < t0_min)
    after = (tmin >= t0_min) & (tmin < t0_min + 20)
    print("\nper-band energy change at event (after-before, +/-20 min):")
    for lbl, _, _ in BANDS:
        d = np.nanmean(edb[lbl][after]) - np.nanmean(edb[lbl][before])
        print(f"  {lbl:18s}: {d:+.2f} dB")

    # figure
    fig, axL = plt.subplots(figsize=(15, 7))
    colors = {"infrasound 1-20Hz": "tab:blue", "low 20-100Hz": "tab:orange",
              "ship 100-500Hz": "tab:green", "EM 60Hz": "tab:red"}
    for lbl in colors:
        axL.plot(tmin, edb[lbl], color=colors[lbl], lw=0.8, label=lbl)
    axL.axvline(t0_min, color="red", ls="--", lw=2)
    axL.set_xlabel(f"min since {FSTART:%H:%M} UTC")
    axL.set_ylabel("band energy (dB re median)")
    axL.legend(loc="upper left")
    axR = axL.twinx()
    axR.plot(tmin, snaprate, color="black", lw=1.6, label="SNAP rate")
    axR.set_ylabel("snap rate/min")
    axR.legend(loc="upper right")
    axL.set_title("CI01 candidate-stimulus bands vs snap rate around the UAP event")
    fig.tight_layout()
    fig.savefig("figures/fig_mechanism.png", dpi=110)
    print("\nwrote figures/fig_mechanism.png")


if __name__ == "__main__":
    main()
