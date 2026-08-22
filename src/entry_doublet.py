#!/usr/bin/env python3
"""Entry-doublet detector (discriminant D2): the water-entry signature of
a transmedium object crossing the air-water interface.

A rigid (or vapor-jacketed) body entering water produces a stereotyped
acoustic DOUBLET (Truscott et al. 2014; JASA 159:2061, 2026):

  1. SLAM  -- a broadband impact shock at t=0 (100 Hz-30 kHz, low spectral
             corner f_c ~= 240/a Hz), rise 10 us-1 ms;
  2. quiet cavity-open interlude (40-60 dB below the pulses);
  3. PINCH-OFF pulse at t_p = (1.7-2.3)*sqrt(a/g) after the slam, followed
     by Minnaert bubble ringing f0 = 3.26*sqrt(1+z/10)/R_b Hz.

The slam-to-pinch delay encodes body scale:  a ~= g*(t_p/2)^2.
t_p = 0.12-0.16 s (a=5 cm) ... 0.38-0.52 s (a=0.5 m).

The detector finds ISOLATED pairs of broadband, low-corner, impulsive
transients separated by t_p in [0.10, 0.55] s, optionally with trailing
low-frequency ringing. Isolation + low-frequency energy are what separate
an entry doublet from odontocete click pairs (many clicks, HF-dominated,
no low corner) and from echosounder/airgun PRIs (fixed long repetition).

This module is a library (used by inject/calibrate and the pilot screen)
and a CLI for a single WAV/array of pressure.
"""

import numpy as np
from scipy import signal

DET_BAND = (300.0, 120000.0)      # broadband impulse band
LF_BAND = (100.0, 5000.0)         # slam low-corner band (clicks lack this)
RING_BAND = (3.0, 150.0)          # Minnaert pinch-off ring
TP_MIN, TP_MAX = 0.10, 0.55       # slam->pinch delay window [s]
MAD_K = 8.0                       # impulse detection threshold
MIN_SEP_S = 0.02                  # min separation between distinct impulses
GUARD_S = 1.5                     # isolation guard each side of a pair
PRE_QUIET_S = 0.6                 # required quiet interval before the slam
MAX_NEIGHBORS = 2                 # comparable-amplitude impulses near pair
NEIGH_FRAC = 0.4                  # "comparable" = >= this frac of slam height
DUR_MAX_MS = 6.0                  # impulses must be short
LOWHIGH_MIN = 2.0                # slam low-pass; clicks HF (0% of clicks pass)
G = 9.81


def _bandpass(x, fs, band):
    sos = signal.butter(4, band, btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x)


def find_impulses(p, fs):
    x = _bandpass(p, fs, DET_BAND)
    env = np.abs(signal.hilbert(x))
    med = np.median(env)
    mad = np.median(np.abs(env - med)) + 1e-30
    thr = med + MAD_K * mad
    peaks, props = signal.find_peaks(env, height=thr,
                                     distance=int(MIN_SEP_S * fs))
    return x, env, peaks, props["peak_heights"], thr


def impulse_features(p, x, fs, i):
    h = int(0.010 * fs)
    a, b = max(0, i - h), min(len(x), i + h)
    w, e = x[a:b], np.abs(signal.hilbert(x[a:b]))
    pk = i - a
    peak = e[pk] if pk < len(e) else e.max()
    above = e >= 0.1 * peak
    lo = hi = pk
    while lo > 0 and above[lo - 1]:
        lo -= 1
    while hi < len(e) - 1 and above[hi + 1]:
        hi += 1
    dur_ms = (hi - lo) / fs * 1e3
    seg = w[max(0, pk - int(0.003 * fs)):pk + int(0.005 * fs)]
    f, pxx = signal.welch(seg, fs=fs, nperseg=min(1024, len(seg)))
    tot = pxx.sum() + 1e-30
    lf_frac = float(pxx[(f >= LF_BAND[0]) & (f <= LF_BAND[1])].sum() / tot)
    centroid = float((f * pxx).sum() / tot)
    # low/high band ratio: slam is low-pass (corner ~240/a Hz) -> ratio >> 1;
    # odontocete clicks are HF (15-50 kHz) -> ratio << 1. This is the
    # empirically validated separator (0% of real clicks exceed 2.0).
    lo = pxx[(f >= 300) & (f <= 8000)].sum()
    hi = pxx[(f >= 15000) & (f <= 60000)].sum() + 1e-30
    return {"dur_ms": round(dur_ms, 3), "lf_frac": round(lf_frac, 3),
            "centroid_hz": round(centroid, 1), "peak": float(peak),
            "lowhigh": float(lo / hi)}


def has_ringing(p, fs, i_pinch):
    """Detect a decaying low-frequency (3-150 Hz) tone after the pinch pulse."""
    a = i_pinch
    b = min(len(p), a + int(1.0 * fs))
    seg = _bandpass(p[a:b], fs, RING_BAND)
    env = np.abs(signal.hilbert(seg))
    if len(env) < int(0.1 * fs):
        return 0.0, 0.0
    # ring = strong, decaying LF energy vs pre-pinch baseline
    pre = _bandpass(p[max(0, a - int(0.5 * fs)):a], fs, RING_BAND)
    base = np.median(np.abs(signal.hilbert(pre))) + 1e-30 if len(pre) else 1e-30
    ring_snr = float(20 * np.log10(env[:int(0.2 * fs)].max() / base))
    # decay: first quarter louder than last quarter
    q = len(env) // 4
    decay = float(20 * np.log10((env[:q].mean() + 1e-30) /
                                (env[-q:].mean() + 1e-30)))
    f, pxx = signal.welch(seg, fs=fs, nperseg=min(4096, len(seg)))
    sel = (f >= RING_BAND[0]) & (f <= RING_BAND[1])
    f0 = float(f[sel][np.argmax(pxx[sel])]) if sel.any() else 0.0
    return ring_snr, f0 if decay > 2 else 0.0


def detect(p, fs):
    """Return a list of entry-doublet candidates with grades."""
    x, env, peaks, heights, thr = find_impulses(p, fs)
    times = peaks / fs
    feats = [impulse_features(p, x, fs, int(i)) for i in peaks]
    out = []
    for a in range(len(peaks)):
        fa = feats[a]
        # SLAM gate: impulsive, low-corner, carries low-frequency energy
        if fa["dur_ms"] > DUR_MAX_MS or fa["lowhigh"] < LOWHIGH_MIN:
            continue
        ha = heights[a]
        # pre-quiet: no comparable impulse in the window before the slam
        pre = (times >= times[a] - PRE_QUIET_S) & (times < times[a] - 0.005)
        if np.any(heights[pre] >= NEIGH_FRAC * ha):
            continue
        for b in range(a + 1, len(peaks)):
            tp = times[b] - times[a]
            if tp < TP_MIN:
                continue
            if tp > TP_MAX:
                break
            fb = feats[b]
            if fb["dur_ms"] > DUR_MAX_MS:
                continue
            # slam must be the local maximum; pinch <= slam (impact loudest)
            if heights[b] > ha:
                continue
            # isolation: count only COMPARABLE-amplitude impulses in the
            # window (ring sub-peaks and small snaps are ignored), excluding
            # the pair itself
            lo, hi = times[a] - GUARD_S, times[b] + GUARD_S
            win = (times >= lo) & (times <= hi)
            comparable = int(np.sum((heights >= NEIGH_FRAC * ha) & win) - 2)
            if comparable > MAX_NEIGHBORS:
                continue
            # interlude quieter than both pulses
            imid0 = int((times[a] + tp * 0.35) * fs)
            imid1 = int((times[a] + tp * 0.65) * fs)
            interlude = np.abs(env[imid0:imid1]).mean() if imid1 > imid0 else 0
            if interlude >= 0.5 * min(ha, heights[b]):
                continue
            ring_snr, f0 = has_ringing(p, fs, int(peaks[b]))
            implied_a_m = G * (tp / 2) ** 2
            grade = ("A" if (ring_snr > 8 and f0 > 0) else "B")
            out.append({
                "t_slam_s": round(float(times[a]), 4),
                "t_pinch_s": round(float(times[b]), 4),
                "tp_s": round(float(tp), 4),
                "implied_body_a_m": round(float(implied_a_m), 3),
                "slam_lowhigh": round(fa["lowhigh"],2),
                "slam_centroid_hz": fa["centroid_hz"],
                "pinch_centroid_hz": fb["centroid_hz"],
                "n_neighbors": comparable,
                "ring_snr_db": round(ring_snr, 1),
                "ring_f0_hz": round(f0, 1),
                "grade": grade,
            })
    # de-duplicate overlapping pairs sharing the slam, keep best-graded
    out.sort(key=lambda d: (d["t_slam_s"], d["grade"]))
    dedup, used = [], set()
    for d in out:
        key = round(d["t_slam_s"], 2)
        if key in used:
            continue
        used.add(key)
        dedup.append(d)
    return dedup


def main():
    import argparse
    import soundfile as sf
    ap = argparse.ArgumentParser()
    ap.add_argument("wav")
    args = ap.parse_args()
    from spectral_analysis import load_pressure
    p, fs = load_pressure(args.wav)
    p = p * 1e-6  # uPa -> Pa
    for d in detect(p, fs):
        print(d)


if __name__ == "__main__":
    main()
