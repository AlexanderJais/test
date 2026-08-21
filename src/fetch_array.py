#!/usr/bin/env python3
"""Fetch OOI Regional Cabled Array low-frequency hydrophone data (an
actual underwater *array*) from EarthScope/IRIS FDSN web services.

Network OO, channel HDH: five bottom-mounted hydrophones at 200 Hz with a
common GPS-disciplined timebase — the property a single hydrophone lacks
and localization requires:

  Axial Seamount cluster : AXBA1, AXCC1, AXEC2   (3-26 km baselines)
  Hydrate Ridge pair     : HYS14, HYSB1          (~21 km baseline)
  total aperture         : ~380 km

Windows fetched:

* EVENT   — 11 min from the origin time of the USGS-reviewed M5.0
  earthquake of 2023-01-11 10:17:18 UTC at 43.9646 N, 128.7389 W
  (Blanco region). Its water-borne T-phase is a loud impulsive source
  with independently known position: ground truth for validating the
  TDOA localization pipeline.
* AMBIENT — two 10-min comparison windows for signature screening.

Instrument response is removed (output in Pa) and each trace is saved to
one .npz per window with station coordinates and start times, so the
downstream scripts need only numpy.

Usage:
    python3 src/fetch_array.py [--out data/array]
"""

import argparse
import os

import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

FDSN_BASE = "https://service.earthscope.org"
NET, CHA = "OO", "HDH"
STATIONS = ["AXBA1", "AXCC1", "AXEC2", "HYS14", "HYSB1"]

EVENT = {                     # USGS us7000j3ld (reviewed)
    "origin": "2023-01-11T10:17:18.204",
    "lat": 43.9646, "lon": -128.7389, "mag": 5.0,
}

WINDOWS = {
    # T-phase arrives ~150-200 s after origin at these ranges; include
    # margin on both sides.
    "event":    (UTCDateTime(EVENT["origin"]) - 60, 660),
    # Quiet-day comparisons for the signature screening (same winter
    # season -> fin whale calls likely; different day/hour).
    "ambient1": (UTCDateTime("2023-01-11T06:00:00"), 600),
    "ambient2": (UTCDateTime("2023-01-14T18:00:00"), 600),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/array")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    client = Client(FDSN_BASE)
    inv = client.get_stations(network=NET, channel=CHA, level="response")
    coords = {}
    for net in inv:
        for sta in net:
            coords[sta.code] = (sta.latitude, sta.longitude,
                                float(sta.elevation))

    for wname, (t0, dur) in WINDOWS.items():
        out = os.path.join(args.out, f"{wname}.npz")
        if os.path.exists(out):
            print(f"cached: {out}")
            continue
        payload = {"t0": str(t0), "duration": float(dur)}
        for code in STATIONS:
            try:
                st = client.get_waveforms(NET, code, "*", CHA,
                                          t0, t0 + dur)
            except Exception as e:
                print(f"  {wname}/{code}: no data ({e})")
                continue
            st.merge(fill_value="interpolate")
            tr = st[0]
            tr.detrend("linear")
            tr.remove_response(inventory=inv, output="DEF",
                               water_level=60,
                               pre_filt=(0.2, 0.5, 90.0, 99.0))
            payload[f"{code}_x"] = tr.data.astype(np.float32)
            payload[f"{code}_fs"] = float(tr.stats.sampling_rate)
            payload[f"{code}_start"] = float(tr.stats.starttime - t0)
            payload[f"{code}_lat"], payload[f"{code}_lon"] = \
                coords[code][0], coords[code][1]
            payload[f"{code}_elev"] = coords[code][2]
            print(f"  {wname}/{code}: {len(tr.data)} samples @ "
                  f"{tr.stats.sampling_rate:g} Hz, "
                  f"rms={np.std(tr.data):.3f} Pa")
        np.savez_compressed(out, **payload)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
