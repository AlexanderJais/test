#!/usr/bin/env python3
"""
Rebuild the "egg craft / hieratic numerals" figures from the two source images.

  source/annotated.webp       the annotated image (chart + 13 red arrows), 1456x1941
  source/post_screenshot.webp the original post; the photo panel is crop x>=463, y>=116

Stages
  1. read the 13 red arrows and trace each one back to its cell in the numeral chart
  2. place the annotated frame onto the photo (transform below, see README)
  3. erase the arrows from the annotated image, restoring surface texture
  4. enhance the hull and render the three panels + positions.csv

usage:  python3 scripts/pipeline.py [outdir]
"""
import cv2, numpy as np, os, sys, math, json, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = sys.argv[1] if len(sys.argv) > 1 else ROOT
SRC  = os.path.join(ROOT, 'source')
os.makedirs(OUT, exist_ok=True)
W, H = 1456, 1941
GOLD, INK, WHT = (38, 196, 255), (24, 24, 28), (252, 252, 252)
FD = cv2.FONT_HERSHEY_DUPLEX

# transform annotated -> photo panel: similarity, recovered by multi-scale template
# matching on the cave-rock texture (17 correspondences, 1.8 px rms).
S_, ROT_, TX_, TY_ = 0.25889, -2.318, 519.29, 389.93
_c, _s = math.cos(math.radians(ROT_)), math.sin(math.radians(ROT_))
A = np.array([[S_*_c, -S_*_s, TX_], [S_*_s, S_*_c, TY_]])

# chart geometry inside the annotated image (detected from its rules)
COLS = {1: (42, 155), 10: (165, 277), 100: (287, 429), 1000: (439, 581)}
HY   = [565, 653, 741, 830, 918, 1007, 1095, 1183, 1272, 1360]
CHART_BOX = (33, 555, 633, 1416)


def red_mask(img, r=100, d=35):
    B, G, R = (img[:, :, i].astype(int) for i in range(3))
    return ((R > r) & (R - G > d) & (R - B > d)).astype(np.uint8)


def find_arrows(a):
    """Locate arrowheads, trace each shaft back, return [(head, tail, value)].

    For each head the shaft direction is the ray (0.05 deg steps) collecting the
    most red pixels within 2.2 px; the tail is the far end of the contiguous red
    run along it.  Global support + contiguity makes this stable where arrows
    cross, which a step-along-the-ray walk is not.
    """
    red = red_mask(a)
    op = cv2.morphologyEx(red, cv2.MORPH_OPEN,
                          cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9)))
    n, _, st, ce = cv2.connectedComponentsWithStats(op, 8)
    heads = sorted([(float(ce[i][0]), float(ce[i][1]))
                    for i in range(1, n) if st[i, cv2.CC_STAT_AREA] >= 30],
                   key=lambda p: p[1])
    ys, xs = np.nonzero(red)
    xs = xs.astype(float); ys = ys.astype(float)

    def cell(x, y):
        c = [k for k, (x0, x1) in COLS.items() if x0 - 6 <= x <= x1 + 6]
        r = [i for i in range(9) if HY[i] - 5 <= y <= HY[i+1] + 5]
        return (r[0] + 1) * c[0] if c and r else None

    out = []
    for h in heads:
        dx, dy = xs - h[0], ys - h[1]
        best = None
        for th in np.arange(0, 2*np.pi, np.pi/3600):
            c_, s_ = math.cos(th), math.sin(th)
            t = dx*c_ + dy*s_; p = np.abs(-dx*s_ + dy*c_)
            k = int(((t > 20) & (t < 1500) & (p < 2.2)).sum())
            if best is None or k > best[0]:
                best = (k, th)
        th = best[1]; c_, s_ = math.cos(th), math.sin(th)
        t = dx*c_ + dy*s_; p = np.abs(-dx*s_ + dy*c_)
        tt = np.sort(t[(t > 20) & (p < 2.2)])
        br = np.nonzero(np.diff(tt) > 25)[0]
        e = tt[br[0]] if len(br) else tt[-1]
        tail = (h[0] + c_*e, h[1] + s_*e)
        out.append((h, tail, cell(*tail)))
    return out


def extract_glyphs(a, values):
    """Crop each hieratic sign out of the chart, dropping label text and arrows."""
    g = {}
    for v in sorted(set(values)):
        mag = 10 ** (len(str(v)) - 1); d = v // mag
        x0, x1 = COLS[mag]; y0, y1 = HY[d-1], HY[d]
        cell = a[y0+4:y1-4, x0+4:x1-4].copy()
        B, G, R = (cell[:, :, i].astype(int) for i in range(3))
        red = ((R-G > 25) & (R-B > 25)) | ((R > 150) & (G < 210) & (B < 210) & (R-G > 12))
        cell[cv2.dilate(red.astype(np.uint8), np.ones((3, 3), np.uint8)) > 0] = 255
        ink = cv2.morphologyEx((cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY) < 125).astype(np.uint8),
                               cv2.MORPH_OPEN, np.ones((2, 2), np.uint8))
        prof = ink.sum(axis=0) == 0                       # widest gap splits label|sign
        runs, s = [], None
        for i, e in enumerate(prof):
            if e and s is None: s = i
            if not e and s is not None: runs.append((i-s, s, i)); s = None
        if s is not None: runs.append((len(prof)-s, s, len(prof)))
        inner = [r for r in runs if r[1] > 4 and r[2] < len(prof)-4]
        ink[:, :max(inner)[2] if inner else 0] = 0
        n, _, st, _ = cv2.connectedComponentsWithStats(ink, 8)
        keep = [i for i in range(1, n) if st[i, cv2.CC_STAT_AREA] >= 12] or list(range(1, n))
        p = 5
        ys0 = max(0, min(st[i, 1] for i in keep) - p); ys1 = max(st[i, 1]+st[i, 3] for i in keep) + p
        xs0 = max(0, min(st[i, 0] for i in keep) - p); xs1 = max(st[i, 0]+st[i, 2] for i in keep) + p
        crop = cell[ys0:ys1, xs0:xs1]
        g[v] = cv2.cvtColor(cv2.cvtColor(crop, cv2.COLOR_BGR2GRAY), cv2.COLOR_GRAY2BGR)
    return g


def erase_arrows(a):
    """Remove the red strokes: unbiased local mean + high-frequency texture donor."""
    af = a.astype(np.float32)
    red = cv2.dilate(red_mask(a, 95, 28), np.ones((5, 5), np.uint8))
    keep = (1 - red).astype(np.float32)
    lo = (cv2.GaussianBlur(af*keep[:, :, None], (0, 0), 6.0) /
          (cv2.GaussianBlur(keep, (0, 0), 6.0)[:, :, None] + 1e-6))
    hf = af - cv2.GaussianBlur(af, (0, 0), 3.0)
    need = red > 0
    filled = np.zeros_like(hf); done = np.zeros(red.shape, bool)
    for dy, dx in ([(0, d) for d in (26, -26, 40, -40, 58, -58, 80, -80)] +
                   [(d, 0) for d in (26, -26, 40, -40, 58, -58)]):
        take = need & ~done & np.roll(~need, (dy, dx), axis=(0, 1))
        if not take.any():
            continue
        filled[take] = np.roll(hf, (dy, dx), axis=(0, 1))[take]; done |= take
    rem = need & ~done
    if rem.any():
        filled[rem] = np.random.default_rng(0).normal(
            0, float(hf[~need].std()), size=(int(rem.sum()), 3))
    soft = cv2.GaussianBlur(red.astype(np.float32), (0, 0), 1.6)[:, :, None]
    return np.clip(af*(1-soft) + (lo+filled)*soft, 0, 255).astype(np.uint8)


def hull_ellipse(photo):
    """Fit the hull outline by radial edge search + RANSAC."""
    g = cv2.GaussianBlur(cv2.cvtColor(photo, cv2.COLOR_BGR2GRAY), (0, 0), 3).astype(np.float32)
    hh, ww = g.shape; c, pts = (490, 700), []
    for k in range(720):
        th = 2*np.pi*k/720; d = np.array([math.cos(th), math.sin(th)])
        rs = np.arange(60, 520, 1.0)
        xs, ys = c[0]+d[0]*rs, c[1]+d[1]*rs
        ok = (xs > 2) & (xs < ww-3) & (ys > 2) & (ys < hh-3)
        if ok.sum() < 40 or not ok[:40].all():
            continue
        last = int(np.argmax(~ok)) if (~ok).any() else len(ok)
        if last < 60:
            continue
        xs, ys = xs[:last], ys[:last]
        prof = cv2.remap(g, xs.astype(np.float32).reshape(-1, 1),
                         ys.astype(np.float32).reshape(-1, 1), cv2.INTER_LINEAR).ravel()
        d1 = np.gradient(prof); i = int(np.argmin(d1))
        if 8 <= i <= len(prof)-6 and -d1[i] >= 1.2:
            pts.append((xs[i], ys[i]))
    XY = np.array(pts, np.float32)

    def dist(e):
        (cx, cy), (MA, ma), ang = e; t = math.radians(ang)
        R = np.array([[math.cos(t), math.sin(t)], [-math.sin(t), math.cos(t)]])
        q = (XY - [cx, cy]) @ R.T
        return (np.sqrt((q[:, 0]/(MA/2))**2 + (q[:, 1]/(ma/2))**2) - 1) * ((MA+ma)/4)
    rng = np.random.default_rng(0); best = None
    for _ in range(4000):
        try: e = cv2.fitEllipse(XY[rng.choice(len(XY), 6, replace=False)])
        except cv2.error: continue
        if not (20 < e[1][0] < 3000 and 20 < e[1][1] < 3000):
            continue
        n = int((np.abs(dist(e)) < 3.0).sum())
        if best is None or n > best[0]: best = (n, e)
    e = best[1]
    for _ in range(12):
        k = np.abs(dist(e)) < 4.5
        if k.sum() < 10: break
        e = cv2.fitEllipse(XY[k])
    return e


def fit_text(t, mw, base=1.0, th=2):
    s = base
    while s > 0.2:
        (tw, tht), _ = cv2.getTextSize(t, FD, s, th)
        if tw <= mw: return s, tw, tht
        s -= 0.02
    return s, 0, 0


def glyph_card(gl, v, cw, ch, bar=True, bw=3):
    card = np.full((ch, cw, 3), 252, np.uint8)
    bh = int(ch*0.26) if bar else 0
    ah, aw = ch-bh-14, cw-20
    s = min(aw/gl.shape[1], ah/gl.shape[0])
    r = cv2.resize(gl, None, fx=s, fy=s, interpolation=cv2.INTER_AREA)
    y0 = 7 + (ah-r.shape[0])//2; x0 = (cw-r.shape[1])//2
    card[y0:y0+r.shape[0], x0:x0+r.shape[1]] = r
    if bar:
        cv2.rectangle(card, (0, ch-bh), (cw, ch), (44, 44, 48), -1)
        fs, tw, tht = fit_text(str(v), cw-14, 0.75, 2)
        cv2.putText(card, str(v), ((cw-tw)//2, ch-(bh-tht)//2-2), FD, fs, WHT, 2, cv2.LINE_AA)
    cv2.rectangle(card, (0, 0), (cw-1, ch-1), GOLD, bw)
    return card


def main():
    a = cv2.imread(os.path.join(SRC, 'annotated.webp'))
    photo = cv2.imread(os.path.join(SRC, 'post_screenshot.webp'))[116:, 463:]
    arrows = find_arrows(a)
    pts1 = np.array([h for h, _, _ in arrows], np.float32)
    vals = [v for _, _, v in arrows]
    print('arrows:', list(zip(range(1, 14), vals)))
    glyphs = extract_glyphs(a, vals)
    ell = hull_ellipse(photo); (ecx, ecy), (eMA, ema), eang = ell

    clean = erase_arrows(a)
    egg = (cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY) > 0).astype(np.uint8)
    _, egg = cv2.threshold(cv2.GaussianBlur(cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY), (9, 9), 0),
                           0, 1, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    egg = cv2.morphologyEx(egg, cv2.MORPH_CLOSE, np.ones((35, 35), np.uint8))
    n, lab, st, _ = cv2.connectedComponentsWithStats(egg, 8)
    egg = (lab == 1 + int(np.argmax(st[1:, cv2.CC_STAT_AREA]))).astype(np.float32)
    egg = cv2.GaussianBlur(cv2.erode(egg, np.ones((9, 9), np.uint8)), (0, 0), 12)

    lb = cv2.cvtColor(clean, cv2.COLOR_BGR2LAB); L = lb[:, :, 0].astype(np.float32)
    D = cv2.GaussianBlur(L, (0, 0), 2.6) - cv2.GaussianBlur(L, (0, 0), 34.)
    L2 = cv2.createCLAHE(1.1, (16, 16)).apply(L.astype(np.uint8)).astype(np.float32)
    o = lb.copy(); o[:, :, 0] = np.clip(L2 + (1.0 + 2.4*egg)*D, 0, 255).astype(np.uint8)
    for ch in (1, 2):
        o[:, :, ch] = np.clip(128 + (o[:, :, ch].astype(np.float32)-128)*0.30, 0, 255).astype(np.uint8)
    enh = cv2.cvtColor(o, cv2.COLOR_LAB2BGR)

    img = enh.copy()
    # ---- legend panel over the old chart footprint ----
    x0, y0, x1, y1 = CHART_BOX; pw, ph = x1-x0, y1-y0
    pan = np.full((ph, pw, 3), 248, np.uint8)
    cv2.rectangle(pan, (0, 0), (pw, 96), (30, 30, 34), -1)
    cv2.putText(pan, 'HIERATIC NUMERALS', (20, 44), FD, 1.0, WHT, 2, cv2.LINE_AA)
    cv2.putText(pan, '13 marked positions on the hull, keyed to the rings',
                (20, 76), FD, 0.5, (178, 198, 220), 1, cv2.LINE_AA)
    cw_, chh = pw//2, (ph-96-30)//7
    for i in range(13):
        cx0 = (i//7)*cw_; cy0 = 96 + (i % 7)*chh
        cv2.rectangle(pan, (cx0+7, cy0+4), (cx0+cw_-7, cy0+chh-4), (224, 224, 228), 1)
        cv2.circle(pan, (cx0+34, cy0+chh//2), 19, (44, 44, 48), -1, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(i+1), 26, 0.62, 1)
        cv2.putText(pan, str(i+1), (cx0+34-tw//2, cy0+chh//2+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
        gh = chh-22; gc = glyph_card(glyphs[vals[i]], vals[i], int(gh*0.92), gh, False, 2)
        gx, gy = cx0+64, cy0+11
        pan[gy:gy+gc.shape[0], gx:gx+gc.shape[1]] = gc
        fs, tw, tht = fit_text(str(vals[i]), cw_-(gx-cx0)-gc.shape[1]-26, 0.95, 2)
        cv2.putText(pan, str(vals[i]), (gx+gc.shape[1]+18, cy0+chh//2+tht//2), FD, fs, (30, 30, 34), 2, cv2.LINE_AA)
    cv2.rectangle(pan, (0, ph-30), (pw, ph), (236, 236, 240), -1)
    cv2.putText(pan, 'glyph shapes taken from the source annotation chart',
                (14, ph-10), FD, 0.44, (88, 88, 94), 1, cv2.LINE_AA)
    cv2.rectangle(pan, (0, 0), (pw-1, ph-1), (58, 58, 64), 2)
    img[y0:y1, x0:x1] = pan

    # ---- rings + label cards ----
    RING, CW, CH = 58, 104, 132
    placed, pos = [], {}
    for i in np.argsort(-pts1[:, 1]):
        px, py = pts1[i]; best = None
        for ang in range(0, 360, 10):
            for rad in (128, 164, 200, 240, 285):
                lx = px + rad*math.cos(math.radians(ang)); ly = py + rad*math.sin(math.radians(ang))
                if lx-CW/2 < 14 or ly-CH/2 < 14 or lx+CW/2 > W-14 or ly+CH/2 > H-14:
                    continue
                pen = rad*0.30
                if ly-CH/2 < 104: pen += 1500
                if (lx-CW/2 < x1+14 and lx+CW/2 > x0-14 and ly-CH/2 < y1+14 and ly+CH/2 > y0-14):
                    pen += 2500
                for ox, oy in placed:
                    if abs(ox-lx) < CW+22 and abs(oy-ly) < CH+22: pen += 1200
                for qx, qy in pts1:
                    if abs(qx-lx) < CW/2+RING+10 and abs(qy-ly) < CH/2+RING+10: pen += 800
                if best is None or pen < best[0]: best = (pen, lx, ly)
        pos[i] = best[1:]; placed.append(best[1:])
    for i in range(13):
        px, py = pts1[i]; lx, ly = pos[i]
        d = np.array([lx-px, ly-py]); d /= np.linalg.norm(d)
        s = (px+d[0]*(RING+3), py+d[1]*(RING+3)); e = (lx-d[0]*CH*0.40, ly-d[1]*CH*0.40)
        cv2.line(img, tuple(np.int32(s)), tuple(np.int32(e)), INK, 6, cv2.LINE_AA)
        cv2.line(img, tuple(np.int32(s)), tuple(np.int32(e)), GOLD, 2, cv2.LINE_AA)
    for px, py in pts1:
        cv2.circle(img, (int(px), int(py)), RING+2, INK, 6, cv2.LINE_AA)
        cv2.circle(img, (int(px), int(py)), RING, GOLD, 3, cv2.LINE_AA)
    for i in range(13):
        lx, ly = pos[i]; cx0, cy0 = int(lx-CW/2), int(ly-CH/2)
        sh = img[cy0+7:cy0+CH+7, cx0+7:cx0+CW+7]
        if sh.shape[:2] == (CH, CW):
            img[cy0+7:cy0+CH+7, cx0+7:cx0+CW+7] = (sh*0.40).astype(np.uint8)
        img[cy0:cy0+CH, cx0:cx0+CW] = glyph_card(glyphs[vals[i]], vals[i], CW, CH)
        cv2.circle(img, (cx0+1, cy0+1), 17, INK, -1, cv2.LINE_AA)
        cv2.circle(img, (cx0+1, cy0+1), 17, GOLD, 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(i+1), 24, 0.55, 1)
        cv2.putText(img, str(i+1), (cx0+1-tw//2, cy0+1+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
    img[0:100, :] = (img[0:100, :]*0.28).astype(np.uint8)
    cv2.putText(img, '"EGG" CRAFT  -  HIERATIC NUMERALS HIGHLIGHTED', (26, 46), FD, 0.95, WHT, 2, cv2.LINE_AA)
    cv2.putText(img, 'surface enhanced (band-pass + local contrast) - annotation arrows digitally '
                'removed - rings mark the 13 annotated positions', (26, 78), FD, 0.45, (188, 204, 222), 1, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, '01_ship_numerals_highlighted.png'), img)

    # ---- 02 detail sheet ----
    TW, TH, GW, cols, R = 250, 250, 150, 4, 125
    rows_ = math.ceil(13/cols)
    sheet = np.full((rows_*(TH+34)+96, cols*(TW+GW+14)+14, 3), 244, np.uint8)
    cv2.rectangle(sheet, (0, 0), (sheet.shape[1], 86), (30, 30, 34), -1)
    cv2.putText(sheet, 'THE 13 MARKED SITES  -  enhanced hull surface vs. the claimed hieratic sign',
                (24, 40), FD, 0.82, WHT, 2, cv2.LINE_AA)
    cv2.putText(sheet, 'left: band-pass enhanced crop, 250 px = ~64 px of the source photo   |   '
                'right: reference glyph + value', (24, 68), FD, 0.46, (180, 200, 220), 1, cv2.LINE_AA)
    for i in range(13):
        r, cc = divmod(i, cols); ox = 14+cc*(TW+GW+14); oy = 96+r*(TH+34)
        x, y = int(pts1[i][0]), int(pts1[i][1])
        pad = np.zeros((2*R, 2*R, 3), np.uint8)
        ys0, xs0 = max(0, y-R), max(0, x-R)
        ys1, xs1 = min(enh.shape[0], y+R), min(enh.shape[1], x+R)
        pad[ys0-(y-R):ys1-(y-R), xs0-(x-R):xs1-(x-R)] = enh[ys0:ys1, xs0:xs1]
        crop = cv2.resize(pad, (TW, TH), interpolation=cv2.INTER_CUBIC)
        cv2.circle(crop, (TW//2, TH//2), int(58*TW/(2*R)), GOLD, 2, cv2.LINE_AA)
        sheet[oy:oy+TH, ox:ox+TW] = crop
        cv2.rectangle(sheet, (ox, oy), (ox+TW+GW, oy+TH), (120, 120, 126), 1)
        gl = glyphs[vals[i]]; box = np.full((TH, GW, 3), 252, np.uint8)
        s = min((GW-30)/gl.shape[1], (TH-86)/gl.shape[0])
        g2 = cv2.resize(gl, None, fx=s, fy=s, interpolation=cv2.INTER_AREA)
        yy = 24+((TH-70-g2.shape[0])//2)
        box[yy:yy+g2.shape[0], (GW-g2.shape[1])//2:(GW-g2.shape[1])//2+g2.shape[1]] = g2
        cv2.rectangle(box, (0, TH-46), (GW, TH), (44, 44, 48), -1)
        fs, tw, tht = fit_text(str(vals[i]), GW-16, 0.9, 2)
        cv2.putText(box, str(vals[i]), ((GW-tw)//2, TH-14), FD, fs, WHT, 2, cv2.LINE_AA)
        sheet[oy:oy+TH, ox+TW:ox+TW+GW] = box
        cv2.circle(sheet, (ox+20, oy+20), 17, INK, -1, cv2.LINE_AA)
        cv2.circle(sheet, (ox+20, oy+20), 17, GOLD, 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(i+1), 24, 0.6, 1)
        cv2.putText(sheet, str(i+1), (ox+20-tw//2, oy+20+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, '02_site_details.png'), sheet)

    # ---- 03 context ----
    P = cv2.transform(pts1.reshape(-1, 1, 2), A).reshape(-1, 2)
    g = cv2.cvtColor(photo, cv2.COLOR_BGR2GRAY)
    den = cv2.fastNlMeansDenoising(g, None, 10, 7, 21).astype(np.float32)
    D2 = cv2.GaussianBlur(den, (0, 0), 1.6) - cv2.GaussianBlur(den, (0, 0), 17.)
    m = np.zeros(g.shape, np.float32)
    cv2.ellipse(m, ((ecx, ecy), (eMA, ema), eang), 1., -1)
    m = cv2.GaussianBlur(m, (0, 0), 9)
    lb2 = cv2.cvtColor(photo, cv2.COLOR_BGR2LAB)
    lb2[:, :, 0] = np.clip(cv2.createCLAHE(1.1, (16, 16)).apply(lb2[:, :, 0]).astype(np.float32)
                           + 3.6*D2*m, 0, 255).astype(np.uint8)
    ctx = cv2.resize(cv2.cvtColor(lb2, cv2.COLOR_LAB2BGR), None, fx=1.5, fy=1.5, interpolation=cv2.INTER_CUBIC)
    S = 1.5
    cv2.ellipse(ctx, ((ecx*S, ecy*S), (eMA*S, ema*S), eang), (120, 200, 120), 2, cv2.LINE_AA)
    fr = cv2.transform(np.array([[[0, 0]], [[W, 0]], [[W, H]], [[0, H]]], np.float32), A).reshape(-1, 2)*S
    cv2.polylines(ctx, [fr.astype(np.int32)], True, (255, 170, 60), 3, cv2.LINE_AA)
    cv2.putText(ctx, 'zoom frame of panel 1', (int(fr[:, 0].min())+10, int(fr[:, 1].min())-12),
                FD, 0.6, (255, 170, 60), 2, cv2.LINE_AA)
    for i, (x, y) in enumerate(P):
        px, py = int(x*S), int(y*S)
        cv2.circle(ctx, (px, py), 19, INK, -1, cv2.LINE_AA)
        cv2.circle(ctx, (px, py), 19, GOLD, 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(i+1), 26, 0.6, 1)
        cv2.putText(ctx, str(i+1), (px-tw//2, py+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
    hdr = np.full((96, ctx.shape[1], 3), (30, 30, 34), np.uint8)
    cv2.putText(hdr, 'CONTEXT  -  where the 13 marks sit on the whole craft', (24, 42), FD, 0.85, WHT, 2, cv2.LINE_AA)
    cv2.putText(hdr, 'green = fitted hull outline   |   orange = the zoomed region shown in panel 1   |   '
                'all 13 lie on the upper-right quadrant', (24, 72), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, '03_context_whole_craft.png'), np.vstack([hdr, ctx]))

    # ---- 04 data ----
    t = math.radians(eang)
    emin = np.array([math.cos(t), math.sin(t)]); emaj = np.array([-math.sin(t), math.cos(t)])
    sa, sb = ema/2.0, eMA/2.0
    recs = []
    for i, (x, y) in enumerate(P):
        d = np.array([x-ecx, y-ecy]); u = float(d@emaj)/sa; w = float(d@emin)/sb
        recs.append(dict(idx=i+1, value=vals[i],
                         img1_x=round(float(pts1[i][0]), 1), img1_y=round(float(pts1[i][1]), 1),
                         photo_x=round(float(x), 1), photo_y=round(float(y), 1),
                         screenshot_x=round(float(x)+463, 1), screenshot_y=round(float(y)+116, 1),
                         u_major=round(u, 4), v_minor=round(w, 4),
                         r_norm=round(math.hypot(u, w), 4),
                         theta_deg=round(math.degrees(math.atan2(w, u)), 2)))
    with open(os.path.join(OUT, 'positions.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(recs[0].keys())); wr.writeheader(); wr.writerows(recs)
    json.dump(dict(
        transform_img1_to_photo=dict(matrix=A.tolist(), scale=S_, rotation_deg=ROT_,
                                     translation=[TX_, TY_], residual_rms_px=1.8,
                                     method='multi-scale template matching on rock texture, 17 correspondences'),
        hull_ellipse_photo_px=dict(center=[round(float(ecx), 1), round(float(ecy), 1)],
                                   semi_major=round(float(sa), 1), semi_minor=round(float(sb), 1),
                                   angle_deg=round(float(eang), 2)),
        values=vals, sum=int(sum(vals)),
        notes=['ring radius 58 px in panel-1 frame = 15 px in the source photo',
               'theta_deg measured from the hull major axis, image y down',
               'r_norm = 1.0 is the fitted hull rim']),
        open(os.path.join(OUT, 'analysis_meta.json'), 'w'), indent=2)
    print('wrote panels + positions.csv to', OUT)


if __name__ == '__main__':
    main()
