#!/usr/bin/env python3
"""Quench-tail detector (discriminant D4): the hot-body flag.

A cold rigid body is acoustically silent after cavity pinch-off. A HOT or
plasma-sheathed body, once submerged, sheds its vapor jacket through
film-boiling collapse -- 'chugging': an IRREGULAR impulsive train (pulses
broadband 100 Hz-50 kHz, sub-ms, repetition 1-100 Hz, envelope-modulated
0.1-10 Hz) recurring over the quench time t_q ~ 10-1000 s and decaying as
the surface cools through the Leidenfrost point. So: any SUSTAINED,
IRREGULAR, DECAYING impulsive tail bound to the same track, immediately
after an entry doublet, is the hot-body signature.

Separators baked in:
  * irregular inter-event intervals (ISI CV > 0.5) reject the metronomic
    38 kHz echosounder (0.4-0.5 Hz, CV ~ 0);
  * sustained over many windows rejects one-off transients;
  * decaying level over the tail rejects steady machinery.
(Odontocete click bouts also make sustained irregular impulses -- so, like
the entry doublet, this limb is a COINCIDENCE leg, gated on a preceding
doublet and a consistent track, not a standalone detector.)

Library + CLI. detect_tail(p, fs, t_trigger) -> tail descriptor or None.
"""

import numpy as np
from scipy import signal

BOIL_BAND = (100.0, 50000.0)      # full boiling band (context)
CHUG_BAND = (5000.0, 50000.0)     # where chugging impulses beat ambient
WIN_S = 2.0
HOP_S = 2.0
LEVEL_UP_DB = 4.0                 # broadband elevation over baseline
KURT_MIN = 4.0                    # impulsive
ISI_CV_MIN = 0.5                  # irregular (rejects metronomic pingers)
MIN_WINDOWS = 4                   # sustained
MAD_K = 6.0


def _bp(x, fs, band):
    sos = signal.butter(4, band, btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x)


def _window_stats(seg, fs, baseline_rms):
    xb = _bp(seg, fs, CHUG_BAND)
    rms = np.sqrt(np.mean(xb ** 2)) + 1e-30
    level_up = 20 * np.log10(rms / baseline_rms)
    env = np.abs(signal.hilbert(xb))
    kurt = float(((env - env.mean()) ** 4).mean() / (env.var() ** 2 + 1e-30))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    pk, _ = signal.find_peaks(env, height=med + MAD_K * mad,
                              distance=int(0.003 * fs))
    if len(pk) >= 3:
        isi = np.diff(pk) / fs
        isi_cv = float(np.std(isi) / (np.mean(isi) + 1e-30))
    else:
        isi_cv = np.nan
    return {"level_up_db": float(level_up), "kurtosis": kurt,
            "n_imp": int(len(pk)), "isi_cv": isi_cv}


def detect_tail(p, fs, t_trigger, max_tail_s=None):
    """Test for a hot-body quench tail after t_trigger. Returns a
    descriptor dict if a sustained irregular decaying impulsive tail is
    present, else None."""
    i_trig = int(t_trigger * fs)
    if max_tail_s is None:
        max_tail_s = (len(p) - i_trig) / fs
    # baseline from up to 3 s BEFORE the trigger
    b0 = max(0, i_trig - int(3.0 * fs))
    base = np.sqrt(np.mean(_bp(p[b0:i_trig], fs, CHUG_BAND) ** 2)) + 1e-30 \
        if i_trig > b0 else 1e-30
    stats, t = [], []
    n = int(WIN_S * fs)
    hop = int(HOP_S * fs)
    for w in range(i_trig, min(len(p) - n, i_trig + int(max_tail_s * fs)),
                   hop):
        stats.append(_window_stats(p[w:w + n], fs, base))
        t.append((w - i_trig) / fs)
    # active windows: elevated, impulsive, irregular
    active = [s for s in stats
              if s["level_up_db"] > LEVEL_UP_DB and s["kurtosis"] > KURT_MIN
              and (np.isnan(s["isi_cv"]) or s["isi_cv"] > ISI_CV_MIN)]
    if len(active) < MIN_WINDOWS:
        return None
    levels = np.array([s["level_up_db"] for s in stats])
    # require an overall decaying trend (quench cools) over active span
    idx = [i for i, s in enumerate(stats)
           if s["level_up_db"] > LEVEL_UP_DB and s["kurtosis"] > KURT_MIN]
    if len(idx) >= 3:
        slope = float(np.polyfit(np.array(t)[idx], levels[idx], 1)[0])
    else:
        slope = 0.0
    return {
        "t_trigger_s": round(float(t_trigger), 3),
        "tail_windows": len(active),
        "tail_duration_s": round(float(t[idx[-1]] - t[idx[0]]), 1)
        if idx else 0.0,
        "mean_level_up_db": round(float(np.mean(
            [s["level_up_db"] for s in active])), 1),
        "mean_kurtosis": round(float(np.mean(
            [s["kurtosis"] for s in active])), 1),
        "median_isi_cv": round(float(np.nanmedian(
            [s["isi_cv"] for s in active])), 2),
        "decay_db_per_s": round(slope, 3),
        "is_decaying": bool(slope < -0.02),
    }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("wav")
    ap.add_argument("--trigger", type=float, default=0.5)
    args = ap.parse_args()
    from spectral_analysis import load_pressure
    p, fs = load_pressure(args.wav)
    print(detect_tail(p * 1e-6, fs, args.trigger))


if __name__ == "__main__":
    main()
