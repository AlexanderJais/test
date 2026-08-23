#!/usr/bin/env python3
"""Fine-grain onset shape of the CI01 snap-rate rise at the Omaha event.

Is the rise a STEP (rise time of seconds -> a discrete trigger at an instant)
or a RAMP (rise time of many minutes -> a gradual driver / behavioural drift)?
And where does the change point sit relative to the event arrival t0?

Method (CI01 event file, 48 kHz):
  * fixed detection threshold from a clean pre-event baseline (median+8*MAD of
    the 2-20 kHz envelope) -- so a genuine RATE change is not normalised away
    by a per-bin adaptive threshold;
  * snap rate in 5 s bins across t0 +/- 50 min;
  * CUSUM change-point (independent of any shape assumption);
  * fit three shapes and compare by BIC:
      step  r = a + b*H(t-tau)                       (rise time 0)
      ramp  r = a, linear a->a+b over [t1,t2], then a+b   (rise time t2-t1)
      tanh  r = a + b*0.5*(1+tanh((t-tau)/s))         (10-90% rise = 2.197*s)
  * report the 10-90% rise time and the onset time relative to t0.
Overlays the per-bin adaptive metric as a robustness check.
"""
import datetime
import numpy as np
import soundfile as sf
from scipy import signal, optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EVENT_FILE = "data/omaha/SanctSound_CI01_02_671883305_190716033643.flac"
FSTART = datetime.datetime(2019, 7, 16, 3, 36, 43)
EVENT_UTC = datetime.datetime(2019, 7, 16, 6, 0, 0)
CI01 = (34.0438, -120.0811)
EVENT_LATLON = (32.4894, -119.3647)
C_KM_S = 1.487
SNAP_BAND = (2000.0, 20000.0)
WIN_S = 5.0                 # fine bin
ZOOM_MIN = 50              # +/- around arrival
BASE_LO_MIN, BASE_HI_MIN = -45, -40   # pre-event baseline for the fixed threshold


def haversine_km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, [a[0], a[1], b[0], b[1]])
    return 2 * 6371 * np.arcsin(np.sqrt(
        np.sin((la2 - la1) / 2) ** 2 +
        np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2))


def main():
    info = sf.info(EVENT_FILE)
    fs = info.samplerate
    sos = signal.butter(4, SNAP_BAND, btype="bandpass", fs=fs, output="sos")
    rng = haversine_km(CI01, EVENT_LATLON)
    t0_s = (EVENT_UTC - FSTART).total_seconds() + rng / C_KM_S
    t0_min = t0_s / 60.0

    def env_of(a_s, b_s):
        with sf.SoundFile(EVENT_FILE) as fh:
            fh.seek(int(a_s * fs))
            x = fh.read(int((b_s - a_s) * fs), dtype="float64")
        if x.ndim > 1:
            x = x[:, 0]
        return np.abs(signal.hilbert(signal.sosfiltfilt(sos, x)))

    # fixed threshold from clean baseline
    be = env_of(t0_s + BASE_LO_MIN * 60, t0_s + BASE_HI_MIN * 60)
    med = np.median(be)
    mad = np.median(np.abs(be - med)) + 1e-30
    thr = med + 8 * mad
    dist = int(0.01 * fs)

    lo_s = max(0.0, t0_s - ZOOM_MIN * 60)
    hi_s = min(info.frames / fs, t0_s + ZOOM_MIN * 60)
    n = int(WIN_S * fs)
    offs = np.arange(int(lo_s * fs), int(hi_s * fs) - n, n)
    tmin = (offs / fs) / 60.0
    rate_fixed = np.empty(len(offs))
    rate_adapt = np.empty(len(offs))
    with sf.SoundFile(EVENT_FILE) as fh:
        for i, off in enumerate(offs):
            fh.seek(int(off))
            x = fh.read(n, dtype="float64")
            if x.ndim > 1:
                x = x[:, 0]
            env = np.abs(signal.hilbert(signal.sosfiltfilt(sos, x)))
            rate_fixed[i] = len(signal.find_peaks(env, height=thr, distance=dist)[0])
            m2 = np.median(env)
            a2 = m2 + 8 * (np.median(np.abs(env - m2)) + 1e-30)
            rate_adapt[i] = len(signal.find_peaks(env, height=a2, distance=dist)[0])
    rate_fixed *= 60.0 / WIN_S      # per minute
    rate_adapt *= 60.0 / WIN_S

    t = tmin - t0_min               # minutes relative to event arrival
    r = rate_fixed

    # --- CUSUM change point (shape-agnostic) ---
    mu = r.mean()
    cs = np.cumsum(r - mu)
    tau_cusum = t[np.argmax(np.abs(cs - cs[0]))]

    # --- model fits ---
    def sse(res):
        return np.sum(res ** 2)

    # step: optimal tau by scan
    best = (np.inf, None)
    for k in range(5, len(t) - 5):
        a = r[:k].mean(); b = r[k:].mean()
        s = sse(r[:k] - a) + sse(r[k:] - b)
        if s < best[0]:
            best = (s, (t[k], a, b - a))
    sse_step, (tau_step, a_step, b_step) = best

    # tanh
    def ftanh(tt, a, b, tau, s):
        return a + b * 0.5 * (1 + np.tanh((tt - tau) / s))
    try:
        p0 = [r[t < -10].mean(), r[t > 10].mean() - r[t < -10].mean(), tau_cusum, 5.0]
        popt, _ = optimize.curve_fit(ftanh, t, r, p0=p0, maxfev=20000,
                                     bounds=([-np.inf, -np.inf, t.min(), 0.05],
                                             [np.inf, np.inf, t.max(), 60]))
        sse_tanh = sse(r - ftanh(t, *popt))
        rise_1090 = 2.197 * popt[3]     # minutes
        tau_tanh = popt[2]
    except Exception as e:
        popt = None; sse_tanh = np.nan; rise_1090 = np.nan; tau_tanh = np.nan
        print("tanh fit failed:", e)

    # ramp (piecewise-linear)
    def framp(tt, a, b, t1, t2):
        w = np.clip((tt - t1) / max(t2 - t1, 1e-6), 0, 1)
        return a + b * w
    try:
        pr0 = [a_step, b_step, tau_step - 5, tau_step + 5]
        propt, _ = optimize.curve_fit(framp, t, r, p0=pr0, maxfev=20000,
                                      bounds=([-np.inf, -np.inf, t.min(), t.min()],
                                              [np.inf, np.inf, t.max(), t.max()]))
        sse_ramp = sse(r - framp(t, *propt))
        ramp_width = propt[3] - propt[2]
    except Exception as e:
        propt = None; sse_ramp = np.nan; ramp_width = np.nan
        print("ramp fit failed:", e)

    N = len(r)
    def bic(sse_, k):
        return N * np.log(sse_ / N) + k * np.log(N)
    bic_step = bic(sse_step, 3)
    bic_tanh = bic(sse_tanh, 4) if np.isfinite(sse_tanh) else np.nan
    bic_ramp = bic(sse_ramp, 4) if np.isfinite(sse_ramp) else np.nan

    print(f"CI01 event arrival t0 = {t0_min:.1f} min; baseline snap rate "
          f"~{a_step:.0f}/min, jump ~{b_step:+.0f}/min ({100*b_step/a_step:+.1f}%)")
    print(f"\nchange-point (CUSUM, shape-agnostic): {tau_cusum:+.1f} min rel. to t0")
    print(f"step  fit: tau={tau_step:+.1f} min           SSE={sse_step:.3e}  BIC={bic_step:.1f}")
    print(f"ramp  fit: window=[{propt[2]:+.1f},{propt[3]:+.1f}] min, width={ramp_width:.1f} min "
          f"SSE={sse_ramp:.3e}  BIC={bic_ramp:.1f}")
    print(f"tanh  fit: mid tau={tau_tanh:+.1f} min, 10-90% rise={rise_1090:.1f} min "
          f"SSE={sse_tanh:.3e}  BIC={bic_tanh:.1f}")
    bics = {"step": bic_step, "ramp": bic_ramp, "tanh": bic_tanh}
    best_model = min((k for k in bics if np.isfinite(bics[k])), key=lambda k: bics[k])
    print(f"\nbest by BIC: {best_model.upper()}  "
          f"(step-best means instantaneous; ramp/tanh best means gradual)")
    if np.isfinite(rise_1090):
        verdict = ("STEP-like (discrete trigger)" if rise_1090 < 2
                   else "RAMP-like (gradual driver)" if rise_1090 > 8
                   else "intermediate")
        print(f"10-90% rise time = {rise_1090:.1f} min  ->  {verdict}")

    # figure
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.plot(t, rate_fixed, color="black", lw=1.0, label="snap rate (fixed thr, 5 s)")
    ax.plot(t, rate_adapt, color="0.6", lw=0.7, label="snap rate (per-bin adaptive)")
    tt = np.linspace(t.min(), t.max(), 800)
    if popt is not None:
        ax.plot(tt, ftanh(tt, *popt), color="tab:red", lw=2,
                label=f"tanh fit (10-90% = {rise_1090:.1f} min)")
    ax.axvline(0, color="tab:blue", ls="--", lw=2, label="event arrival t0")
    ax.axvline(tau_cusum, color="tab:green", ls=":", lw=2,
               label=f"CUSUM change pt ({tau_cusum:+.1f} min)")
    ax.set_xlabel("minutes relative to event arrival")
    ax.set_ylabel("snap rate / min")
    ax.set_title("CI01 snap-rate onset, fine-grain (5 s bins)")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig("figures/fig_onset_timing.png", dpi=115)
    print("\nwrote figures/fig_onset_timing.png")


if __name__ == "__main__":
    main()
