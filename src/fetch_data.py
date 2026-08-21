#!/usr/bin/env python3
"""Fetch short slices of MBARI Pacific Sound 256 kHz hydrophone data.

The full-resolution archive lives in public S3 buckets named
``pacific-sound-256khz-YYYY`` (10-minute WAV files, ~460 MB each, 24-bit
PCM at 256,000 samples/s, single channel).  Downloading whole files is
wasteful for spectral survey work, so this script:

1. lists the bucket with an anonymous HTTPS request to find the file
   covering a requested timestamp,
2. range-GETs the RIFF header plus the first N seconds of PCM data,
3. patches the RIFF/data chunk sizes so the result is a valid,
   self-contained WAV slice.

No AWS credentials or SDK required — plain HTTPS byte-range requests.

Usage:
    python3 src/fetch_data.py [--out DIR] [--seconds 60]
"""

import argparse
import os
import re
import struct
import sys
import urllib.request

BUCKET_FMT = "https://pacific-sound-256khz-{year}.s3.amazonaws.com"

# Survey slices: spread across seasons and times of day (UTC) in 2022 so the
# narrowband-line persistence analysis can separate permanent instrument
# lines from transient tonal sources (ships, biologics, echosounders).
DEFAULT_SLICES = [
    ("2022", "01", "20220101_00"),  # winter, night (UTC 00 ≈ 16:00 local PST)
    ("2022", "03", "20220315_12"),  # spring, UTC noon (04:00 local, true night)
    ("2022", "05", "20220510_06"),  # late spring, UTC 06 (22:00 local)
    ("2022", "07", "20220720_18"),  # summer, UTC 18 (11:00 local, daytime)
    ("2022", "09", "20220905_03"),  # early fall, UTC 03 (20:00 local)
    ("2022", "11", "20221125_21"),  # late fall, UTC 21 (13:00 local, daytime)
]


def http_get(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def list_keys(year, month, name_prefix):
    """List bucket keys under MM/ starting with MARS_<name_prefix>."""
    base = BUCKET_FMT.format(year=year)
    prefix = f"{month}/MARS_{name_prefix}"
    xml = http_get(f"{base}/?list-type=2&max-keys=20&prefix={prefix}").decode()
    return re.findall(r"<Key>([^<]+)</Key>", xml)


def parse_riff_header(blob):
    """Return (data_chunk_offset, fs, channels, bits) from a RIFF blob."""
    assert blob[:4] == b"RIFF" and blob[8:12] == b"WAVE", "not a WAV file"
    pos, fs, ch, bits = 12, None, None, None
    while pos + 8 <= len(blob):
        cid = blob[pos:pos + 4]
        csize = struct.unpack("<I", blob[pos + 4:pos + 8])[0]
        if cid == b"fmt ":
            ch, fs = struct.unpack("<HI", blob[pos + 10:pos + 16])
            bits = struct.unpack("<H", blob[pos + 22:pos + 24])[0]
        elif cid == b"data":
            return pos + 8, fs, ch, bits
        pos += 8 + csize + (csize & 1)
    raise ValueError("no data chunk found in first block")


def fetch_slice(year, month, name_prefix, seconds, out_dir):
    keys = list_keys(year, month, name_prefix)
    if not keys:
        print(f"  !! no file found for {name_prefix}", file=sys.stderr)
        return None
    key = keys[0]
    base = BUCKET_FMT.format(year=year)
    url = f"{base}/{key}"

    head = http_get(url, {"Range": "bytes=0-4095"})
    data_off, fs, ch, bits = parse_riff_header(head)
    bytes_per_s = fs * ch * (bits // 8)
    need = data_off + seconds * bytes_per_s

    out = os.path.join(out_dir, os.path.basename(key).replace(
        ".wav", f"_first{seconds}s.wav"))
    if os.path.exists(out) and os.path.getsize(out) == need:
        print(f"  cached: {out}")
        return out

    blob = bytearray(http_get(url, {"Range": f"bytes=0-{need - 1}"}))
    # Patch RIFF size and data-chunk size so the slice is a valid WAV.
    blob[4:8] = struct.pack("<I", need - 8)
    blob[data_off - 4:data_off] = struct.pack("<I", need - data_off)
    with open(out, "wb") as f:
        f.write(blob)
    print(f"  fetched {key}: fs={fs} Hz, {bits}-bit, {ch} ch -> {out} "
          f"({need / 1e6:.1f} MB)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.environ.get("PACIFIC_SOUND_DIR",
                                                    "data/slices"))
    ap.add_argument("--seconds", type=int, default=60)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    paths = []
    for year, month, name_prefix in DEFAULT_SLICES:
        print(f"slice {name_prefix} ...")
        p = fetch_slice(year, month, name_prefix, args.seconds, args.out)
        if p:
            paths.append(p)
    print(f"\n{len(paths)} slices ready in {args.out}")


if __name__ == "__main__":
    main()
