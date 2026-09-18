#!/usr/bin/env python3
"""
Check the detector against the hand-drawn sketches supplied by the user.

Their magenta outlines are loose indications of the 3-D shape, not exact
boundaries, and they only drew the marks they could identify - so they are used
to test DETECTION (did we find the mark, and roughly its extent), never to tune
the shape morphology.  Fitting the outlines directly produces blobs.

usage:  python3 scripts/validate_user.py [outdir]
"""
import cv2, numpy as np, os, sys, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engravings import surfaces, detect, ROOT, C_SHADOW, C_LIGHT, FD
OUT = sys.argv[1] if len(sys.argv) > 1 else ROOT
SRC = os.path.join(ROOT, 'source')
WHT = (252, 252, 252)

def main():
    clean, enh, egg, valid = surfaces()
    SD, HL, marks, cl = detect(clean, valid)
    uo = json.load(open(os.path.join(SRC, 'user_outlines.json')))
    canvas = np.zeros(valid.shape, np.uint8)
    for pts in sum(uo.values(), []):
        for x, y in pts:
            if 0 <= y < canvas.shape[0] and 0 <= x < canvas.shape[1]:
                canvas[y, x] = 255
    ct, _ = cv2.findContours(cv2.dilate(canvas, np.ones((3, 3), np.uint8)),
                             cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    filled = np.zeros_like(canvas); cv2.drawContours(filled, ct, -1, 255, -1)
    uf = filled > 0
    zone = (cv2.dilate(uf.astype(np.uint8), np.ones((61, 61), np.uint8)) > 0) & (valid > 0)
    det = marks > 0
    hit = (det & uf & zone).sum()
    print('within the sketched zone:')
    print('  recall of the sketched area : %.3f' % (hit/max((uf & zone).sum(), 1)))
    print('  of what we detect there, inside a sketch: %.3f' % (hit/max((det & zone).sum(), 1)))
    found = 0
    for gid, polys in uo.items():
        m = np.zeros(valid.shape, bool)
        for pts in polys:
            for x, y in pts:
                if 0 <= y < m.shape[0] and 0 <= x < m.shape[1]:
                    m[max(0, y-30):y+30, max(0, x-30):x+30] = True
        if (det & m).sum() > 120: found += 1
    print('  sketched marks with a detection: %d / %d' % (found, len(uo)))

    tiles = []
    ids = sorted(uo, key=lambda k: int(k))
    for gid in ids:
        pts = np.concatenate(uo[gid]); cx, cy = int(pts[:, 0].mean()), int(pts[:, 1].mean())
        R = int(max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))/2)+30; T = 230
        def cut(im, f=0):
            o = np.full((2*R, 2*R)+((3,) if im.ndim == 3 else ()), f, im.dtype)
            sy0, sx0 = max(0, cy-R), max(0, cx-R)
            sy1, sx1 = min(im.shape[0], cy+R), min(im.shape[1], cx+R)
            o[sy0-(cy-R):sy1-(cy-R), sx0-(cx-R):sx1-(cx-R)] = im[sy0:sy1, sx0:sx1]
            return o
        ph = cv2.resize(cut(enh), (T, T), interpolation=cv2.INTER_CUBIC)
        sd = cv2.resize(cut(SD*255), (T, T), interpolation=cv2.INTER_CUBIC) > 110
        hl = cv2.resize(cut(HL*255), (T, T), interpolation=cv2.INTER_CUBIC) > 110
        fg = cv2.resize(cut(marks*255), (T, T), interpolation=cv2.INTER_CUBIC) > 110
        ov = ph.copy()
        ov[hl] = (0.55*ov[hl] + 0.45*np.array(C_LIGHT)).astype(np.uint8)
        ov[sd] = (0.55*ov[sd] + 0.45*np.array(C_SHADOW)).astype(np.uint8)
        ov[fg] = (0.25*ov[fg] + 0.75*np.array((70, 210, 120))).astype(np.uint8)
        us = ph.copy()
        um = cv2.resize(cut(canvas), (T, T), interpolation=cv2.INTER_NEAREST) > 60
        us[um] = (200, 0, 220)
        row = np.hstack([ph, ov, us])
        cv2.rectangle(row, (0, 0), (row.shape[1]-1, T-1), (150, 150, 156), 1)
        cv2.putText(row, '#'+gid, (6, 20), FD, 0.55, WHT, 2, cv2.LINE_AA)
        tiles.append(row)
    cols = 3; TWt, THt = tiles[0].shape[1], tiles[0].shape[0]
    rows = (len(tiles)+cols-1)//cols
    sheet = np.full((rows*(THt+12)+104, cols*(TWt+12)+12, 3), 244, np.uint8)
    cv2.rectangle(sheet, (0, 0), (sheet.shape[1], 92), (30, 30, 34), -1)
    cv2.putText(sheet, 'CHECK AGAINST THE HAND SKETCHES', (24, 40), FD, 0.85, WHT, 2, cv2.LINE_AA)
    cv2.putText(sheet, 'per mark: enhanced photo  |  detected relief (red = shadow, blue = lit side)  |  '
                'the hand sketch, in magenta', (24, 70), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    for i, t in enumerate(tiles):
        r, cc = divmod(i, cols)
        sheet[104+r*(THt+12):104+r*(THt+12)+THt, 12+cc*(TWt+12):12+cc*(TWt+12)+TWt] = t
    cv2.imwrite(os.path.join(OUT, '05_check_vs_sketches.png'), sheet)
    print('wrote 05_check_vs_sketches.png')

if __name__ == '__main__':
    main()
