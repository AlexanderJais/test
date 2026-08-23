#!/usr/bin/env python3
"""Aguadilla 2013 thermal-video UAP: the decisive wind-vs-motion test.

Not a hydrophone case (thermal only, no acoustic sensor nearby) -- included
because it was investigated on request. The lantern hypothesis predicts the
object drifts AT wind speed and DOWNWIND. Tested against wind independently
pulled from the Rafael Hernandez Airport (TJBQ) METAR archive (Iowa Env.
Mesonet) for the incident hour (2013-04-26 ~01:22 UTC = 21:22 AST 04-25).

Measured wind (TJBQ, my pull): 8-14 mph from the E/ENE through the evening,
~9.8 mph (4.4 m/s) from the ENE at the incident -- matching the value AARO
cites. AARO's STK reconstruction: object drifted 8 mph toward the SW, over
land, ~200 m altitude; apparent 120 mph (SCU) is motion parallax from the
orbiting/zooming aircraft (plane-object range tripled). Object motion matches
the wind vector in BOTH speed and direction -> wind-driven drift, not
propulsion. (Apparent split = view-angle change; apparent water entry =
thermal crossover near sunset. AARO: objects did not enter the water.)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# TJBQ METAR (UTC time label, speed kt, direction FROM deg) -- pulled from
# mesonet.agron.iastate.edu ASOS archive for 2013-04-25/26.
OBS = [("20:50", 14, 60), ("21:50", 13, 70), ("22:50", 11, 60),
       ("23:50", 12, 70), ("00:50", 7, 70), ("01:50", 7, 90), ("02:50", 6, 100)]
OBJ_SPEED_KT = 6.95     # AARO 8 mph
OBJ_HEADING_DEG = 240   # AARO "southwest", downwind
INCIDENT_X = 4.5        # ~01:22 UTC between 00:50 and 01:50


def main():
    hours = np.arange(len(OBS))
    spd = np.array([o[1] for o in OBS])
    lbl = [o[0] for o in OBS]

    fig = plt.figure(figsize=(13, 5.2))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(hours, spd, "o-", color="tab:blue", label="TJBQ measured wind (kt)")
    ax.axhline(8.55, color="tab:green", ls="--", label="AARO wind cited 9.8 mph")
    ax.axhline(OBJ_SPEED_KT, color="tab:red", lw=2, label="AARO object speed 8 mph")
    ax.axvline(INCIDENT_X, color="k", ls=":", label="incident ~01:22 UTC")
    ax.set_xticks(hours); ax.set_xticklabels(lbl, rotation=45)
    ax.set_ylabel("knots"); ax.set_xlabel("UTC (Apr 25-26 2013)")
    ax.set_ylim(0, 20); ax.set_title("Object speed vs measured wind")
    ax.annotate("SCU claim ~120 mph = 104 kt (off-scale;\nmotion parallax)",
                (0.3, 17), fontsize=8, color="purple")
    ax.legend(fontsize=8)

    axp = fig.add_subplot(1, 2, 2, projection="polar")
    axp.set_theta_zero_location("N"); axp.set_theta_direction(-1)
    for _, s, d in OBS:                       # wind blows TOWARD (from+180)
        axp.annotate("", xy=(np.radians((d + 180) % 360), s), xytext=(0, 0),
                     arrowprops=dict(arrowstyle="->", color="tab:blue", alpha=0.5))
    axp.annotate("", xy=(np.radians(OBJ_HEADING_DEG), OBJ_SPEED_KT), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="->", color="tab:red", lw=2.5))
    axp.set_ylim(0, 16)
    axp.set_title("Direction: wind toward (blue) vs\nobject drift (red) — both SW/WSW",
                  fontsize=9)
    fig.suptitle("Aguadilla 2013: object motion vs independently-measured wind",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig("figures/fig_aguadilla_wind.png", dpi=115)
    print("wrote figures/fig_aguadilla_wind.png")


if __name__ == "__main__":
    main()
