#!/usr/bin/env python3
"""Doppler-drift tracker (discriminants D1/D3): the kinematic core.

A tone (or PRF comb line) from a source transiting past the hydrophone
sweeps a monotonic S-curve through closest point of approach (CPA):
up-shifted plateau on approach, inflection at CPA, down-shifted plateau on
recession. For a straight constant-velocity pass,

    f_obs(t) = f0 * c / (c + v_r(t)),   v_r(t) = v^2 (t - t_c) / R(t),
    R(t) = sqrt(R_cpa^2 + v^2 (t - t_c)^2),   c = 1480 m/s

Fitting the ridge gives (f0, v, R_cpa, t_c). The total fractional swing
gives speed v; the CPA slope df/dt = -f0 v^2 / (c R_cpa) gives range.
Worked: f0=10 kHz, R_cpa=1 km -> 4.5 Hz/s at 50 kn, 17.9 at 100 kn,
286 at 400 kn.

Discriminator vs whale whistles / drifting ship lines: the S-curve is a
specific plateau->inflection->plateau shape; a good fit with a significant,
monotonic swing and a physically bounded v is required. The single most
decisive check (applied at M2/M3) is that ALL comb harmonics share one
fractional S-curve; a lone tone here is graded, not confirmed.

Library + CLI: detect(p, fs, band) -> best S-curve fit or None.
"""

import numpy as np
from scipy import signal, optimize

C_WATER = 1480.0
KN = 0.514444
NPERSEG = 8192
V_MAX = 500 * KN                  # physical search ceiling (~500 kn)


def scurve(t, f0, v, r_cpa, t_c):
    dt = t - t_c
    R = np.sqrt(r_cpa ** 2 + (v * dt) ** 2)
    v_r = (v ** 2 * dt) / R
    return f0 * C_WATER / (C_WATER + v_r)


def track_ridge(p, fs, band, min_prom_db=6.0):
    """Per-frame narrowband peak in `band` -> (t, f, snr) ridge points."""
    f, t, sxx = signal.spectrogram(p, fs=fs, window="hann", nperseg=NPERSEG,
                                   noverlap=NPERSEG // 2)
    sel = (f >= band[0]) & (f <= band[1])
    fb, S = f[sel], sxx[sel, :]
    tt, ff, ss = [], [], []
    for k in range(S.shape[1]):
        col = 10 * np.log10(S[:, k] + 1e-30)
        base = np.median(col)
        i = int(np.argmax(col))
        if col[i] - base >= min_prom_db:
            # parabolic sub-bin refinement
            if 0 < i < len(col) - 1:
                d = col[i - 1] - 2 * col[i] + col[i + 1]
                off = 0.5 * (col[i - 1] - col[i + 1]) / d if d else 0.0
            else:
                off = 0.0
            df = fb[1] - fb[0]
            tt.append(t[k]); ff.append(fb[i] + off * df)
            ss.append(col[i] - base)
    return np.array(tt), np.array(ff), np.array(ss)


def fit_scurve(t, f, snr):
    if len(t) < 12:
        return None
    f0_0 = float(np.median(f))
    span = t[-1] - t[0]
    best = None
    for v0 in (30 * KN, 80 * KN, 200 * KN):
        for tc0 in (t[0] + 0.35 * span, t[0] + 0.5 * span, t[0] + 0.65 * span):
            try:
                popt, _ = optimize.curve_fit(
                    scurve, t, f, p0=[f0_0, v0, 500.0, tc0],
                    bounds=([f0_0 * 0.8, 1.0, 20.0, t[0] - span],
                            [f0_0 * 1.2, V_MAX, 20000.0, t[-1] + span]),
                    sigma=1.0 / np.sqrt(np.maximum(snr, 1.0)), maxfev=8000)
                resid = float(np.sqrt(np.mean((scurve(t, *popt) - f) ** 2)))
                if best is None or resid < best["resid"]:
                    best = {"f0": popt[0], "v_ms": popt[1],
                            "r_cpa_m": popt[2], "t_c": popt[3],
                            "resid_hz": resid}
            except Exception:
                continue
    return best


def detect(p, fs, band, min_swing_frac=3e-3, max_resid_frac=2e-3):
    """Return an S-curve fit if a tone in `band` follows a Doppler transit."""
    t, f, snr = track_ridge(p, fs, band)
    fit = fit_scurve(t, f, snr) if len(t) else None
    if not fit:
        return None
    # observed fractional swing and monotonic-through-inflection shape
    swing = (f.max() - f.min()) / fit["f0"]
    slope_cpa = -fit["f0"] * fit["v_ms"] ** 2 / (C_WATER * fit["r_cpa_m"])
    good = (swing >= min_swing_frac
            and fit["resid_hz"] / fit["f0"] <= max_resid_frac
            and fit["v_ms"] < V_MAX)
    fit.update({"v_kn": round(fit["v_ms"] / KN, 1),
                "swing_frac": round(float(swing), 5),
                "cpa_slope_hz_s": round(slope_cpa, 2),
                "n_ridge": int(len(t)),
                "is_transit": bool(good)})
    return fit


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("wav")
    ap.add_argument("--lo", type=float, default=25000.0)
    ap.add_argument("--hi", type=float, default=45000.0)
    args = ap.parse_args()
    from spectral_analysis import load_pressure
    p, fs = load_pressure(args.wav)
    print(detect(p * 1e-6, fs, (args.lo, args.hi)))


if __name__ == "__main__":
    main()
