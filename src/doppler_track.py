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


def track_ridge(p, fs, band, min_prom_db=8.0, max_jump_hz=400.0):
    """Continuity-constrained narrowband ridge in `band`.

    A real transiting tone is a SINGLE continuous, prominent ridge. Taking
    the per-frame global argmax (old behaviour) tracks broadband noise when
    no tone is present -- huge apparent swing, garbage fit. Instead: seed at
    the most prominent frame, then grow forward/back frame-to-frame within
    max_jump_hz, keeping only frames whose local peak clears min_prom_db.
    Returns (t, f, snr) for the tracked ridge plus a tonality score = the
    fraction of the band's time frames that the ridge actually spans."""
    f, t, sxx = signal.spectrogram(p, fs=fs, window="hann", nperseg=NPERSEG,
                                   noverlap=NPERSEG // 2)
    sel = (f >= band[0]) & (f <= band[1])
    fb, S = f[sel], sxx[sel, :]
    dbf = 10 * np.log10(S + 1e-30)
    base = np.median(dbf, axis=0)                    # per-frame noise floor
    prom = dbf - base                                # prominence map
    nfr = dbf.shape[1]
    df = fb[1] - fb[0]
    jump = max(1, int(round(max_jump_hz / df)))
    # seed at the single most prominent (frame, bin)
    seed_bin, seed_fr = np.unravel_index(np.argmax(prom), prom.shape)
    if prom[seed_bin, seed_fr] < min_prom_db:
        return np.array([]), np.array([]), np.array([]), 0.0
    ridge = {seed_fr: seed_bin}

    def grow(step):
        cur = seed_bin
        k = seed_fr + step
        while 0 <= k < nfr:
            lo, hi = max(0, cur - jump), min(len(fb), cur + jump + 1)
            j = lo + int(np.argmax(prom[lo:hi, k]))
            if prom[j, k] < min_prom_db:
                break                                # ridge ended
            ridge[k] = j
            cur = k = k
            cur = j
            k += step
    grow(+1); grow(-1)
    ks = sorted(ridge)
    tt, ff, ss = [], [], []
    for k in ks:
        j = ridge[k]
        if 0 < j < len(fb) - 1:
            d = dbf[j - 1, k] - 2 * dbf[j, k] + dbf[j + 1, k]
            off = 0.5 * (dbf[j - 1, k] - dbf[j + 1, k]) / d if d else 0.0
        else:
            off = 0.0
        tt.append(t[k]); ff.append(fb[j] + off * df); ss.append(prom[j, k])
    tonality = len(ks) / nfr
    return np.array(tt), np.array(ff), np.array(ss), float(tonality)


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


def detect(p, fs, band, min_swing_frac=3e-3, max_resid_frac=2e-3,
           min_tonality=0.25):
    """Return an S-curve fit if a tone in `band` follows a Doppler transit."""
    t, f, snr, tonality = track_ridge(p, fs, band)
    # a real transit is a continuous tone spanning much of the window; a
    # short/broken ridge is noise and is not fit at all
    if len(t) == 0 or tonality < min_tonality:
        return None
    fit = fit_scurve(t, f, snr)
    if not fit:
        return None
    # observed fractional swing and monotonic-through-inflection shape
    swing = (f.max() - f.min()) / fit["f0"]
    slope_cpa = -fit["f0"] * fit["v_ms"] ** 2 / (C_WATER * fit["r_cpa_m"])
    good = (swing >= min_swing_frac
            and fit["resid_hz"] / fit["f0"] <= max_resid_frac
            and fit["v_ms"] < V_MAX)
    fit.update({"v_kn": round(fit["v_ms"] / KN, 1),
                "tonality": round(tonality, 3),
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
