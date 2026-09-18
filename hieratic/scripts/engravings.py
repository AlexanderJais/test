#!/usr/bin/env python3
"""
Detect and trace the engravings on the hull, without assuming any script.

Method
  1. arrows erased and the hull enhanced (shared with pipeline.py)
  2. Hessian ridge filter (Frangi) run at BOTH polarities, sigma 3/5/7/9.
     These marks are relief: raking light gives each one a shadow on one side
     and a highlight on the other, ~4 px apart.  Taking only the dark ridge
     traces the shadow and throws away half the mark, so both are detected and
     the pair is closed into one solid body.
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
from channel import channel_between

OUT = sys.argv[1] if len(sys.argv) > 1 else ROOT
SRC = os.path.join(ROOT, 'source')
os.makedirs(OUT, exist_ok=True)
SCALE = 0.25889                      # this frame -> source-photo pixels
WHT = (252, 252, 252)
TIER = {'strong': (60, 190, 255), 'moderate': (110, 175, 110), 'weak': (150, 150, 150)}


def frangi_dark(I, scales=(3.0, 5.0, 7.0, 9.0), beta=0.5):
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
    # cuts on rel_contrast: the contrast of the figure's RIMS over the local
    # ridge-response floor.  Calibrated to the observed spread (1.9 - 5.5); the
    # floor is raised by including the highlight channel, so these sit lower
    # than a shadow-only run would put them.
    return 'strong' if r >= 4.3 else ('moderate' if r >= 3.1 else 'weak')


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
LIGHT_AZ = 180.0        # hull shading puts the light to the left (+-20 deg)
PAIR_AXIS = -5.0        # measured mean shadow->highlight direction, near +x
SCAN_SPAN, SCAN_STEP = 45, 15    # pair across a fan of directions, not just one
CHAN_W = 22             # widest rim-to-rim gap counted as one channel, px
RIM_PCT = 92.0          # percentile of the ridge response that counts as a rim
C_SHADOW = (70, 70, 255)     # BGR, shadow side
C_LIGHT = (255, 190, 60)     # BGR, lit rim
C_FIG = (36, 36, 42)         # the incised figure, on white
C_FIG_HI = (70, 210, 120)    # the incised figure, over the photo


def detect(clean, valid):
    """Trace the incised FIGURE, not its rims.

    Each mark is a sunken channel whose two raised rims take the light
    differently: one rim is lit, the other shadowed.  A ridge filter finds the
    rims - but the engraving is the channel BETWEEN them.  So rims are detected
    at both polarities and then paired: a pixel belongs to the figure when a
    shadowed rim lies behind it and a lit rim ahead of it along the illumination
    axis, closer together than CHAN_W.  Pairing runs over a fan of directions so
    channels of any orientation are recovered, not only those square to the light.

    Returns (shadow rims, lit rims, figure, groups).
    """
    g = cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY).astype(np.float32)
    pre = cv2.GaussianBlur(g, (0, 0), 1.8)
    dI = pre - cv2.GaussianBlur(pre, (0, 0), 26)
    Vd = frangi_dark(dI) * valid            # shadowed rim
    Vb = frangi_dark(-dI) * valid           # lit rim
    V = np.maximum(Vd, Vb)
    w = valid.astype(np.float32)
    vm = V.copy(); vm[valid == 0] = 0
    bg = cv2.GaussianBlur(vm, (0, 0), 55)/(cv2.GaussianBlur(w, (0, 0), 55)+1e-6)
    rel = V/(bg+1e-6); rel[valid == 0] = 0
    pc = lambda p: np.percentile(V[valid > 0], p)

    seed = ((V > pc(99.2)) & (valid > 0)).astype(np.uint8)
    grow = ((V > pc(92.0)) & (valid > 0)).astype(np.uint8)
    n, lab, st, _ = cv2.connectedComponentsWithStats(grow, 8)
    has = np.zeros(n, bool)
    for l in np.unique(lab[seed > 0]):
        if l: has[l] = True
    core = np.isin(lab, np.nonzero(has)[0]).astype(np.uint8); core[lab == 0] = 0
    env = cv2.morphologyEx(core, cv2.MORPH_CLOSE,
                           cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    env[valid == 0] = 0
    n2, l2, s2, _ = cv2.connectedComponentsWithStats(env, 8)
    keep = np.zeros_like(env)
    md, mb = Vd > pc(92.0), Vb > pc(92.0)
    for i in range(1, n2):
        m = (l2 == i)
        if s2[i, cv2.CC_STAT_AREA] < 150: continue
        if not (m & md).any() or not (m & mb).any():
            continue                        # relief shows BOTH a lit and a shadowed rim
        keep[m] = 1
    gate = cv2.dilate(keep, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15)))

    def rim(Vx):
        m = ((Vx > pc(RIM_PCT)) & (gate > 0)).astype(np.uint8)
        m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))
        nn, ll, ss, _ = cv2.connectedComponentsWithStats(m, 8)
        out = np.zeros_like(m)
        for k in range(1, nn):
            if ss[k, cv2.CC_STAT_AREA] >= 25: out[ll == k] = 1
        return out
    SD, HL = rim(Vd), rim(Vb)

    fig = np.zeros_like(SD)
    for a in range(-SCAN_SPAN, SCAN_SPAN+1, SCAN_STEP):
        fig |= channel_between(SD, HL, PAIR_AXIS+a, CHAN_W)
    fig = (fig * valid).astype(np.uint8)
    fig = cv2.morphologyEx(fig, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    nf, lf, sf, _ = cv2.connectedComponentsWithStats(fig, 8)
    out = np.zeros_like(fig)
    for k in range(1, nf):
        if sf[k, cv2.CC_STAT_AREA] >= 60: out[lf == k] = 1
    fig = out

    gl = cv2.dilate(fig, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35, 35)))
    cn, clab, _, _ = cv2.connectedComponentsWithStats(gl, 8)
    cl = []
    for k in range(1, cn):
        mm = (clab == k) & (fig > 0)
        if mm.sum() < 260: continue
        ys, xs = np.nonzero(mm)
        nn, _, _, _ = cv2.connectedComponentsWithStats(mm.astype(np.uint8), 8)
        # confidence = contrast of the RIMS that define this figure, against the
        # local response floor.  Measured on the rims, not the filled channel,
        # so it does not shift when the channel geometry changes.
        near = cv2.dilate(mm.astype(np.uint8), np.ones((11, 11), np.uint8)) > 0
        ev = near & ((SD > 0) | (HL > 0))
        if not ev.any(): ev = mm
        cl.append(dict(px=int(mm.sum()), n_strokes=int(nn-1),
                       x0=int(xs.min()), x1=int(xs.max()), y0=int(ys.min()), y1=int(ys.max()),
                       cx=float(xs.mean()), cy=float(ys.mean()), relief='sunken',
                       V=float(V[ev].mean()), rel=float(rel[ev].mean())))
    cl.sort(key=lambda c: -c['rel'])
    for k, c in enumerate(cl): c['id'] = k+1
    return SD, HL, fig, cl


def plural(n, w):
    return '%d %s%s' % (n, w, '' if n == 1 else 's')


def main():
    clean, enh, egg, valid = surfaces()
    SD, HL, strokes, cl = detect(clean, valid)
    print('%d mark groups, %d stroke px' % (len(cl), int(strokes.sum())))
    ell = hull_ellipse(cv2.imread(os.path.join(SRC, 'post_screenshot.webp'))[116:, 463:])
    (ecx, ecy), (eMA, ema), eang = ell
    x0, y0, x1, y1 = CHART_BOX2

    fig = strokes
    def rims_on(bg, sd, hl, alpha=0.75):
        out = bg.copy()
        h, d = hl > 0, sd > 0
        out[h] = ((1-alpha)*out[h] + alpha*np.array(C_LIGHT)).astype(np.uint8)
        out[d] = ((1-alpha)*out[d] + alpha*np.array(C_SHADOW)).astype(np.uint8)
        return out

    # ---------- 01 highlighted ----------
    img = enh.copy()
    img = rims_on(img, SD, HL, 0.34)                 # rims kept faint, for reference
    m = fig > 0
    img[m] = (0.20*img[m] + 0.80*np.array(C_FIG_HI)).astype(np.uint8)
    ed = cv2.morphologyEx(fig, cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)) > 0
    img[ed] = (250, 250, 250)
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
    img[0:118, :] = (img[0:118, :]*0.26).astype(np.uint8)
    cv2.putText(img, 'ENGRAVINGS ON THE HULL  -  the incised figures', (26, 44), FD, 0.95, WHT, 2, cv2.LINE_AA)
    cv2.putText(img, '%d figures.  Each mark is a sunken channel; its two raised rims take the light differently, '
                'and the figure is the channel between them.' % len(cl), (26, 74), FD, 0.45, (188, 204, 222), 1, cv2.LINE_AA)
    lx = 26
    for t, col in (('incised figure', C_FIG_HI), ('shadowed rim', C_SHADOW), ('lit rim', C_LIGHT)):
        cv2.rectangle(img, (lx, 88), (lx+22, 104), col, -1)
        cv2.putText(img, t, (lx+30, 101), FD, 0.44, (206, 216, 230), 1, cv2.LINE_AA)
        lx += int(cv2.getTextSize(t, FD, 0.44, 1)[0][0])+70
    cv2.imwrite(os.path.join(OUT, '01_engravings_highlighted.png'), img)

    # ---------- 02 plate ----------
    TW = TH = 250; PAD = 24; tiles = []
    for c in cl:
        bw, bh = c['x1']-c['x0'], c['y1']-c['y0']
        R = int(max(bw, bh)/2)+PAD; cx, cy = (c['x0']+c['x1'])//2, (c['y0']+c['y1'])//2
        def cut(im, fill=0):
            o = np.full((2*R, 2*R)+((3,) if im.ndim == 3 else ()), fill, im.dtype)
            sy0, sx0 = max(0, cy-R), max(0, cx-R)
            sy1, sx1 = min(im.shape[0], cy+R), min(im.shape[1], cx+R)
            o[sy0-(cy-R):sy1-(cy-R), sx0-(cx-R):sx1-(cx-R)] = im[sy0:sy1, sx0:sx1]
            return o
        ph = cv2.resize(cut(enh), (TW, TH), interpolation=cv2.INTER_CUBIC)
        sd = cv2.resize(cut(SD*255), (TW, TH), interpolation=cv2.INTER_CUBIC) > 110
        hl = cv2.resize(cut(HL*255), (TW, TH), interpolation=cv2.INTER_CUBIC) > 110
        fg = cv2.resize(cut(fig*255), (TW, TH), interpolation=cv2.INTER_CUBIC) > 110
        rim_pane = rims_on(ph, sd.astype(np.uint8), hl.astype(np.uint8), 0.78)
        fac = np.full((TH, TW, 3), 252, np.uint8); fac[fg] = C_FIG
        t = np.hstack([ph, rim_pane, fac])
        for xx in (TW, 2*TW): cv2.line(t, (xx, 0), (xx, TH), (150, 150, 156), 1)
        cv2.rectangle(t, (0, 0), (t.shape[1]-1, TH-1), (150, 150, 156), 1)
        hd = np.full((34, t.shape[1], 3), 248, np.uint8)
        cv2.circle(hd, (19, 17), 13, (34, 34, 38), -1, cv2.LINE_AA)
        cv2.circle(hd, (19, 17), 13, TIER[tier(c['rel'])], 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(c['id']), 20, 0.5, 1)
        cv2.putText(hd, str(c['id']), (19-tw//2, 17+tht//2), FD, fs, WHT, 1, cv2.LINE_AA)
        cv2.putText(hd, '%s   %s   %.0f x %.0f px on the photo' % (
            tier(c['rel']), plural(c['n_strokes'], 'part'), bw*SCALE, bh*SCALE),
            (40, 23), FD, 0.42, (60, 60, 66), 1, cv2.LINE_AA)
        tiles.append(np.vstack([hd, t]))
    cols = 2; TWt, THt = tiles[0].shape[1], tiles[0].shape[0]
    rows = math.ceil(len(tiles)/cols)
    plate = np.full((rows*(THt+14)+118, cols*(TWt+14)+14, 3), 244, np.uint8)
    cv2.rectangle(plate, (0, 0), (plate.shape[1], 104), (30, 30, 34), -1)
    cv2.putText(plate, 'ENGRAVING PLATE  -  %d incised figures' % len(cl), (24, 42), FD, 0.85, WHT, 2, cv2.LINE_AA)
    cv2.putText(plate, 'enhanced photo  |  the two rims of the channel  |  the incised figure between them.     '
                'Sizes are as the mark appears on the source photograph.', (24, 72), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    lx = 24
    for t, col in (('shadowed rim', C_SHADOW), ('lit rim', C_LIGHT)):
        cv2.rectangle(plate, (lx, 84), (lx+20, 98), col, -1)
        cv2.putText(plate, t, (lx+28, 96), FD, 0.42, (180, 200, 220), 1, cv2.LINE_AA)
        lx += int(cv2.getTextSize(t, FD, 0.42, 1)[0][0])+64
    for i, t in enumerate(tiles):
        r, cc = divmod(i, cols)
        plate[118+r*(THt+14):118+r*(THt+14)+THt, 14+cc*(TWt+14):14+cc*(TWt+14)+TWt] = t
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
    fac[fig > 0] = C_FIG
    for c in cl:
        cx, cy = int(c['cx']), int(c['cy'])
        r = int(max(c['x1']-c['x0'], c['y1']-c['y0'])/2)+26
        bx, by = cx+int(r*0.72), cy-int(r*0.72)
        cv2.circle(fac, (bx, by), 16, (245, 245, 248), -1, cv2.LINE_AA)
        cv2.circle(fac, (bx, by), 16, (120, 120, 128), 2, cv2.LINE_AA)
        fs, tw, tht = fit_text(str(c['id']), 22, 0.5, 1)
        cv2.putText(fac, str(c['id']), (bx-tw//2, by+tht//2), FD, fs, (40, 40, 46), 1, cv2.LINE_AA)
    cv2.rectangle(fac, (0, 0), (W, 92), (30, 30, 34), -1)
    cv2.putText(fac, 'FACSIMILE  -  incised figures in hull position', (26, 42), FD, 0.9, WHT, 2, cv2.LINE_AA)
    cv2.putText(fac, 'the channel between the rims, at the scale and position of panel 1',
                (26, 72), FD, 0.45, (180, 200, 220), 1, cv2.LINE_AA)
    cv2.imwrite(os.path.join(OUT, '03_engravings_facsimile.png'), fac)

    # ---------- 04 sign forms ----------
    SH = 170; items = []
    for c in cl:
        if tier(c['rel']) == 'weak': continue
        sl = (slice(max(0, c['y0']-5), c['y1']+6), slice(max(0, c['x0']-5), c['x1']+6))
        sub = (fig[sl]*255).astype(np.uint8)
        sc_ = SH/sub.shape[0]; wpx = max(14, int(sub.shape[1]*sc_))
        sub = cv2.resize(sub, (wpx, SH), interpolation=cv2.INTER_CUBIC) > 110
        im = np.full((SH, wpx, 3), 252, np.uint8); im[sub] = C_FIG
        items.append((c, im))
    GAP, MAXW = 28, 1900
    rws = [[]]; cur = 0
    for c, im in items:
        if cur+im.shape[1]+GAP > MAXW and rws[-1]:
            rws.append([]); cur = 0
        rws[-1].append((c, im)); cur += im.shape[1]+GAP
    strip = np.full((len(rws)*(SH+54)+96, MAXW+40, 3), 252, np.uint8)
    cv2.rectangle(strip, (0, 0), (strip.shape[1], 84), (30, 30, 34), -1)
    cv2.putText(strip, 'SIGN FORMS  -  normalised to equal height for comparison', (24, 38), FD, 0.82, WHT, 2, cv2.LINE_AA)
    cv2.putText(strip, 'strong + moderate figures.  Proportions preserved, absolute size is not.',
                (24, 66), FD, 0.44, (180, 200, 220), 1, cv2.LINE_AA)
    for r, row in enumerate(rws):
        xx = 20; yy = 96+r*(SH+54)
        for c, im in row:
            strip[yy:yy+SH, xx:xx+im.shape[1]] = im
            cv2.putText(strip, '#%d' % c['id'], (xx, yy+SH+26), FD, 0.55, (70, 70, 76), 1, cv2.LINE_AA)
            xx += im.shape[1]+GAP
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
                          relief=c.get('relief') or 'unclear',
                          r_norm=round(math.hypot(u, v), 4),
                          theta_deg=round(math.degrees(math.atan2(v, u)), 2)))
    with open(os.path.join(OUT, 'engravings.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rowsD[0].keys())); wr.writeheader(); wr.writerows(rowsD)
    print('wrote 4 panels + engravings.csv to', OUT)


if __name__ == '__main__':
    main()
