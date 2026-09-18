#!/usr/bin/env python3
"""
Detect and trace the engravings on the hull, without assuming any script.

Method
  1. arrows erased and the hull enhanced (shared with pipeline.py)
  2. Hessian ridge filter (Frangi, dark-line polarity) at sigma 2.5/3.5/4.5 -
     the grooves are 4-10 px wide in this frame, the residual scan-line noise
     is 1-2 px, so scale selection separates them
  3. thresholds are taken RELATIVE to a local ridge-response floor, because the
     striping raises that floor unevenly across the hull
  4. hysteresis -> geodesic reconstruction -> oriented closing, so each groove
     comes out as one connected stroke rather than a dashed line
  5. strokes grouped into mark groups, ranked by contrast against the local floor

usage:  python3 scripts/engravings.py [outdir]
"""
import cv2, numpy as np, os, sys, math, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline import (erase_arrows, hull_ellipse, A, CHART_BOX, fit_text, FD,
                      ROOT, W, H)

OUT = sys.argv[1] if len(sys.argv) > 1 else ROOT
SRC = os.path.join(ROOT, 'source')
os.makedirs(OUT, exist_ok=True)
SCALE = 0.25889                      # this frame -> source-photo pixels
WHT = (252, 252, 252)
TIER = {'strong': (60, 190, 255), 'moderate': (110, 175, 110), 'weak': (150, 150, 150)}


def frangi_dark(I, scales=(2.5, 3.5, 4.5), beta=0.5):
    """Ridge response for dark line structures on a brighter background."""
    out = np.zeros_like(I, np.float32)
    for s in scales:
        L = cv2.GaussianBlur(I, (0, 0), s)
        Lxx = cv2.Sobel(L, cv2.CV_32F, 2, 0, ksize=5)*(s**2)
        Lyy = cv2.Sobel(L, cv2.CV_32F, 0, 2, ksize=5)*(s**2)
        Lxy = cv2.Sobel(L, cv2.CV_32F, 1, 1, ksize=5)*(s**2)
        tmp = np.sqrt(np.maximum((Lxx-Lyy)**2 + 4*Lxy**2, 0))
        l1, l2 = 0.5*(Lxx+Lyy+tmp), 0.5*(Lxx+Lyy-tmp)
        sw = np.abs(l1) > np.abs(l2)
        lo, hi = np.where(sw, l2, l1), np.where(sw, l1, l2)
        S = np.sqrt(lo**2 + hi**2); c = 0.5*S.max()
        Rb = np.abs(lo)/(np.abs(hi)+1e-6)
        Vv = np.exp(-Rb**2/(2*beta**2))*(1-np.exp(-S**2/(2*c**2+1e-9)))
        Vv[hi <= 0] = 0                       # dark grooves only
        out = np.maximum(out, Vv)
    return out


def tier(r):
    return 'strong' if r >= 8 else ('moderate' if r >= 5 else 'weak')


def surfaces():
    a = cv2.imread(os.path.join(SRC, 'annotated.webp'))
    clean = erase_arrows(a)
    g = cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY)
    _, egg = cv2.threshold(cv2.GaussianBlur(g, (9, 9), 0), 0, 1,
                           cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    egg = cv2.morphologyEx(egg, cv2.MORPH_CLOSE, np.ones((35, 35), np.uint8))
    n, lab, st, _ = cv2.connectedComponentsWithStats(egg, 8)
    egg = (lab == 1 + int(np.argmax(st[1:, cv2.CC_STAT_AREA]))).astype(np.uint8)
    soft = cv2.GaussianBlur(cv2.erode(egg, np.ones((9, 9), np.uint8)), (0, 0), 12).astype(np.float32)

    lb = cv2.cvtColor(clean, cv2.COLOR_BGR2LAB); L = lb[:, :, 0].astype(np.float32)
    D = cv2.GaussianBlur(L, (0, 0), 2.6) - cv2.GaussianBlur(L, (0, 0), 34.)
    L2 = cv2.createCLAHE(1.1, (16, 16)).apply(L.astype(np.uint8)).astype(np.float32)
    o = lb.copy(); o[:, :, 0] = np.clip(L2 + (1.0+2.4*soft)*D, 0, 255).astype(np.uint8)
    for ch in (1, 2):
        o[:, :, ch] = np.clip(128+(o[:, :, ch].astype(np.float32)-128)*0.30, 0, 255).astype(np.uint8)
    enh = cv2.cvtColor(o, cv2.COLOR_LAB2BGR)

    valid = cv2.erode(egg, np.ones((3, 3), np.uint8), iterations=9)
    x0, y0, x1, y1 = CHART_BOX2
    valid[y0:y1, x0:x1] = 0               # the source chart hides this patch
    return clean, enh, egg, valid


CHART_BOX2 = (22, 545, 645, 1428)         # chart footprint, with margin


def detect(clean, valid):
    g = cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY).astype(np.float32)
    pre = cv2.GaussianBlur(g, (0, 0), 1.8)
    V = frangi_dark(pre - cv2.GaussianBlur(pre, (0, 0), 26)) * valid
    vm = V.copy(); vm[valid == 0] = 0; w = valid.astype(np.float32)
    bg = cv2.GaussianBlur(vm, (0, 0), 55)/(cv2.GaussianBlur(w, (0, 0), 55)+1e-6)
    rel = V/(bg+1e-6); rel[valid == 0] = 0

    hiA = np.percentile(V[valid > 0], 99.3); loA = np.percentile(V[valid > 0], 97.5)
    seed = ((rel > 4.2) & (V > hiA) & (valid > 0)).astype(np.uint8)
    grow = cv2.morphologyEx(((rel > 2.6) & (V > loA) & (valid > 0)).astype(np.uint8),
                            cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))
    n, lab, st, _ = cv2.connectedComponentsWithStats(grow, 8)
    has = np.zeros(n, bool)
    for l in np.unique(lab[seed > 0]):
        if l: has[l] = True
    core = np.zeros_like(grow)
    for i in range(1, n):
        if not has[i]: continue
        if st[i, cv2.CC_STAT_AREA] < 55 or max(st[i, cv2.CC_STAT_WIDTH], st[i, cv2.CC_STAT_HEIGHT]) < 15:
            continue
        ys, xs = np.nonzero(lab == i)
        ev = np.linalg.eigvalsh(np.cov(np.stack([xs-xs.mean(), ys-ys.mean()])))
        if ev[0] < 1e-6 or ev[1]/max(ev[0], 1e-6) < 4.0:       # must be elongated
            continue
        core[lab == i] = 1

    perm = ((rel > 2.2) & (V > np.percentile(V[valid > 0], 97.0)) & (valid > 0)).astype(np.uint8)
    m = core.copy(); k = np.ones((3, 3), np.uint8)
    for _ in range(400):                                        # geodesic reconstruction
        d = cv2.dilate(m, k) & perm
        if (d == m).all(): break
        m = d
    br = m.copy()                                               # bridge gaps along the stroke
    for ang in range(0, 180, 15):
        L = 9
        ker = np.zeros((L, L), np.uint8); cv2.line(ker, (0, L//2), (L-1, L//2), 1, 1)
        M = cv2.getRotationMatrix2D((L/2-0.5, L/2-0.5), ang, 1.0)
        kr = (cv2.warpAffine(ker*255, M, (L, L)) > 90).astype(np.uint8)
        if kr.sum() >= 3:
            br |= cv2.morphologyEx(m, cv2.MORPH_CLOSE, kr)
    s = cv2.morphologyEx(br, cv2.MORPH_OPEN, np.ones((2, 2), np.uint8))
    n2, l2, _, _ = cv2.connectedComponentsWithStats(s, 8)
    ok = set(np.unique(l2[core > 0])); ok.discard(0)
    strokes = np.isin(l2, list(ok)).astype(np.uint8); strokes[valid == 0] = 0

    gl = cv2.dilate(strokes, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35, 35)))
    cn, clab, _, _ = cv2.connectedComponentsWithStats(gl, 8)
    cl = []
    for i in range(1, cn):
        mm = (clab == i) & (strokes > 0)
        if mm.sum() < 200: continue
        ys, xs = np.nonzero(mm)
        cl.append(dict(px=int(mm.sum()), n_strokes=int(len(set(l2[mm].ravel()))),
                       x0=int(xs.min()), x1=int(xs.max()), y0=int(ys.min()), y1=int(ys.max()),
                       cx=float(xs.mean()), cy=float(ys.mean()),
                       V=float(V[mm].mean()), rel=float(rel[mm].mean())))
    cl.sort(key=lambda c: -c['rel'])
    for i, c in enumerate(cl): c['id'] = i+1
    return V, rel, strokes, cl


def plural(n, w):
    return '%d %s%s' % (n, w, '' if n == 1 else 's')


def main():
    clean, enh, egg, valid = surfaces()
    V, rel, strokes, cl = detect(clean, valid)
    print('%d mark groups, %d stroke px' % (len(cl), int(strokes.sum())))
    ell = hull_ellipse(cv2.imread(os.path.join(SRC, 'post_screenshot.webp'))[116:, 463:])
    (ecx, ecy), (eMA, ema), eang = ell
    x0, y0, x1, y1 = CHART_BOX2

    # ---------- 01 highlighted ----------
    img = enh.copy()
    col = np.zeros_like(img)
    for c in cl:
        m = np.zeros((H, W), np.uint8); m[c['y0']-2:c['y1']+3, c['x0']-2:c['x1']+3] = 1
        col[(strokes > 0) & (m > 0)] = TIER[tier(c['rel'])]
    halo = cv2.GaussianBlur(cv2.dilate(strokes, np.ones((7, 7), np.uint8)).astype(np.float32), (0, 0), 5)[:, :, None]
    img = np.clip(img.astype(np.float32)*(1-0.45*halo) + col.astype(np.float32)*0.95*halo, 0, 255).astype(np.uint8)
    img[strokes > 0] = col[strokes > 0]
    cv2.rectangle(img, (x0, y0), (x1, y1), (42, 42, 48), -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), (92, 92, 102), 2)
    for t, yy, sz, th in (('SURFACE NOT VISIBLE', (y0+y1)//2-18, 0.78, 2),
                          ('(covered by the chart in the source image)', (y0+y1)//2+18, 0.5, 1)):
        fs, tw, _ = fit_text(t, x1-x0-60, sz, th)
        cv2.putText(img, t, (x0+((x1-x0)-tw)//2, yy), FD, fs, (150, 150, 160), th, cv2.LINE_AA)
    for c in cl:
        cx, cy = int(c['cx']), int(c['cy'])
        r = int(max(c['x1']-c['x0'], c['y1']-c['y0'])/2)+26
        cv2.circle(img, (cx, cy), r, TIER[tier(c['rel'])], 2, cv2.LINE_AA)
        bx, by = cx+int(r*0.72), cy-int(r*0.72)
        cv2.circle(img, (bx, by), 17, (24, 24, 28), -1, cv2.LINE_AA)
        cv2.circle(img, (bx, by), 17, TIER[tier(c['rel'])], 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(c['id']), 24, 0.55, 1)
        cv2.putText(img, str(c['id']), (bx-tw//2, by+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
    img[0:104, :] = (img[0:104, :]*0.26).astype(np.uint8)
    cv2.putText(img, 'ENGRAVINGS ON THE HULL  -  detected and traced', (26, 46), FD, 0.95, WHT, 2, cv2.LINE_AA)
    cv2.putText(img, '%d mark groups from a Hessian ridge detector on the raking-light relief.  '
                'No script assumed.' % len(cl), (26, 76), FD, 0.45, (188, 204, 222), 1, cv2.LINE_AA)
    lx = 26
    for t, k in (('strong', 'strong'), ('moderate', 'moderate'), ('weak / possible artefact', 'weak')):
        cv2.circle(img, (lx+8, 94), 7, TIER[k], -1, cv2.LINE_AA)
        cv2.putText(img, t, (lx+22, 99), FD, 0.42, (200, 210, 225), 1, cv2.LINE_AA)
        lx += int(cv2.getTextSize(t, FD, 0.42, 1)[0][0])+58
    cv2.imwrite(os.path.join(OUT, '01_engravings_highlighted.png'), img)

    # ---------- 02 plate ----------
    TW = TH = 250; PAD = 18; tiles = []
    for c in cl:
        bw, bh = c['x1']-c['x0'], c['y1']-c['y0']
        R = int(max(bw, bh)/2)+PAD; cx, cy = (c['x0']+c['x1'])//2, (c['y0']+c['y1'])//2
        def crop(im, fill):
            o = np.full((2*R, 2*R)+((3,) if im.ndim == 3 else ()), fill, im.dtype)
            sy0, sx0 = max(0, cy-R), max(0, cx-R)
            sy1, sx1 = min(im.shape[0], cy+R), min(im.shape[1], cx+R)
            o[sy0-(cy-R):sy1-(cy-R), sx0-(cx-R):sx1-(cx-R)] = im[sy0:sy1, sx0:sx1]
            return o
        ph = cv2.resize(crop(enh, 0), (TW, TH), interpolation=cv2.INTER_CUBIC)
        tr = cv2.resize(crop((strokes*255).astype(np.uint8), 0), (TW, TH), interpolation=cv2.INTER_CUBIC)
        tr = cv2.cvtColor(255-((tr > 110)*255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
        t = np.hstack([ph, tr])
        cv2.line(t, (TW, 0), (TW, TH), (150, 150, 156), 1)
        cv2.rectangle(t, (0, 0), (t.shape[1]-1, TH-1), (150, 150, 156), 1)
        hd = np.full((34, t.shape[1], 3), 248, np.uint8)
        cv2.circle(hd, (19, 17), 13, (34, 34, 38), -1, cv2.LINE_AA)
        cv2.circle(hd, (19, 17), 13, TIER[tier(c['rel'])], 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(c['id']), 20, 0.5, 1)
        cv2.putText(hd, str(c['id']), (19-tw//2, 17+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
        cv2.putText(hd, '%s   %s   %.0f x %.0f px on the photo' % (
            tier(c['rel']), plural(c['n_strokes'], 'stroke'), bw*SCALE, bh*SCALE),
            (40, 23), FD, 0.42, (60, 60, 66), 1, cv2.LINE_AA)
        tiles.append(np.vstack([hd, t]))
    cols = 3; TWt, THt = tiles[0].shape[1], tiles[0].shape[0]
    rows = math.ceil(len(tiles)/cols)
    plate = np.full((rows*(THt+14)+112, cols*(TWt+14)+14, 3), 244, np.uint8)
    cv2.rectangle(plate, (0, 0), (plate.shape[1], 98), (30, 30, 34), -1)
    cv2.putText(plate, 'ENGRAVING PLATE  -  %d mark groups, photo vs. traced shape' % len(cl),
                (24, 42), FD, 0.85, WHT, 2, cv2.LINE_AA)
    cv2.putText(plate, 'left: enhanced surface    right: extracted shape.    Sizes are as the mark '
                'appears on the source photograph.', (24, 72), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    for i, t in enumerate(tiles):
        r, cc = divmod(i, cols)
        plate[112+r*(THt+14):112+r*(THt+14)+THt, 14+cc*(TWt+14):14+cc*(TWt+14)+TWt] = t
    cv2.imwrite(os.path.join(OUT, '02_engraving_plate.png'), plate)

    # ---------- 03 facsimile ----------
    fac = np.full((H, W, 3), 252, np.uint8)
    Ai = np.linalg.inv(np.vstack([A, [0, 0, 1]]))[:2]
    ctr = cv2.transform(np.array([[[ecx, ecy]]], np.float32), Ai).reshape(2)
    cv2.ellipse(fac, ((float(ctr[0]), float(ctr[1])), (eMA/SCALE, ema/SCALE), eang+2.318),
                (206, 206, 212), 3, cv2.LINE_AA)
    cv2.rectangle(fac, (x0, y0), (x1, y1), (233, 233, 237), -1)
    cv2.rectangle(fac, (x0, y0), (x1, y1), (206, 206, 212), 2)
    fs, tw, _ = fit_text('surface not visible', x1-x0-80, 0.7, 2)
    cv2.putText(fac, 'surface not visible', (x0+((x1-x0)-tw)//2, (y0+y1)//2), FD, fs, (170, 170, 176), 2, cv2.LINE_AA)
    fac[strokes > 0] = (20, 20, 24)
    for c in cl:
        cx, cy = int(c['cx']), int(c['cy'])
        r = int(max(c['x1']-c['x0'], c['y1']-c['y0'])/2)+26
        bx, by = cx+int(r*0.72), cy-int(r*0.72)
        cv2.circle(fac, (bx, by), 16, (245, 245, 248), -1, cv2.LINE_AA)
        cv2.circle(fac, (bx, by), 16, (120, 120, 128), 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(c['id']), 22, 0.5, 1)
        cv2.putText(fac, str(c['id']), (bx-tw//2, by+tht//2), FD, fs, (40, 40, 46), 1, cv2.LINE_AA)
    cv2.rectangle(fac, (0, 0), (W, 92), (30, 30, 34), -1)
    cv2.putText(fac, 'FACSIMILE  -  traced engravings in hull position', (26, 42), FD, 0.9, WHT, 2, cv2.LINE_AA)
    cv2.putText(fac, 'tracing only, same scale and position as panel 1', (26, 72), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, '03_engravings_facsimile.png'), fac)

    # ---------- 04 sign forms ----------
    SH = 150; items = []
    for c in cl:
        if tier(c['rel']) == 'weak': continue
        sub = (strokes[c['y0']-4:c['y1']+5, c['x0']-4:c['x1']+5]*255).astype(np.uint8)
        s = SH/sub.shape[0]
        im = cv2.resize(sub, (max(12, int(sub.shape[1]*s)), SH), interpolation=cv2.INTER_CUBIC)
        items.append((c, cv2.cvtColor(255-((im > 110)*255).astype(np.uint8), cv2.COLOR_GRAY2BGR)))
    GAP, MAXW = 26, 1900
    rws = [[]]; cur = 0
    for c, im in items:
        if cur+im.shape[1]+GAP > MAXW and rws[-1]:
            rws.append([]); cur = 0
        rws[-1].append((c, im)); cur += im.shape[1]+GAP
    strip = np.full((len(rws)*(SH+52)+96, MAXW+40, 3), 252, np.uint8)
    cv2.rectangle(strip, (0, 0), (strip.shape[1], 84), (30, 30, 34), -1)
    cv2.putText(strip, 'SIGN FORMS  -  normalised to equal height for comparison', (24, 38), FD, 0.82, WHT, 2, cv2.LINE_AA)
    cv2.putText(strip, 'strong + moderate groups only; proportions preserved, absolute size is not',
                (24, 66), FD, 0.44, (180, 200, 220), 1, cv2.LINE_AA)
    for r, row in enumerate(rws):
        x = 20; y = 96+r*(SH+52)
        for c, im in row:
            strip[y:y+SH, x:x+im.shape[1]] = im
            cv2.putText(strip, '#%d' % c['id'], (x, y+SH+26), FD, 0.55, (70, 70, 76), 1, cv2.LINE_AA)
            x += im.shape[1]+GAP
    cv2.imwrite(os.path.join(OUT, '04_sign_forms.png'), strip)

    # ---------- data ----------
    t = math.radians(eang)
    emin = np.array([math.cos(t), math.sin(t)]); emaj = np.array([-math.sin(t), math.cos(t)])
    sa, sb = ema/2., eMA/2.
    rowsD = []
    for c in cl:
        p = cv2.transform(np.array([[[c['cx'], c['cy']]]], np.float32), A).reshape(2)
        d = np.array([p[0]-ecx, p[1]-ecy]); u = float(d@emaj)/sa; v = float(d@emin)/sb
        rowsD.append(dict(id=c['id'], confidence=tier(c['rel']), rel_contrast=round(c['rel'], 2),
                          n_strokes=c['n_strokes'],
                          img1_x0=c['x0'], img1_y0=c['y0'],
                          img1_w=c['x1']-c['x0'], img1_h=c['y1']-c['y0'],
                          photo_x=round(float(p[0]), 1), photo_y=round(float(p[1]), 1),
                          photo_w=round((c['x1']-c['x0'])*SCALE, 1),
                          photo_h=round((c['y1']-c['y0'])*SCALE, 1),
                          r_norm=round(math.hypot(u, v), 4),
                          theta_deg=round(math.degrees(math.atan2(v, u)), 2)))
    with open(os.path.join(OUT, 'engravings.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rowsD[0].keys())); wr.writeheader(); wr.writerows(rowsD)
    print('wrote 4 panels + engravings.csv to', OUT)


if __name__ == '__main__':
    main()
