import cv2, numpy as np, math

def height_from_shading(dI, az_deg=180.0, lam=0.08, eps=1e-6):
    """Recover a height field from one directional light.

    Lambertian, low slope:  dI ~ d h / d u   with u the light azimuth.
    In Fourier that is  dI_hat = i (k.u) h_hat, so h_hat = dI_hat / (i k.u).
    The division blows up for frequencies perpendicular to u - those are the
    directions the light carries no information about - so it is damped by
    lam*|k|^2 (a smoothness prior).  Height ALONG u is measured; height ACROSS
    u is interpolated by the prior.  That is the honest limit of one light.
    """
    H, W = dI.shape
    ky = np.fft.fftfreq(H)[:, None]*2*np.pi
    kx = np.fft.fftfreq(W)[None, :]*2*np.pi
    u = np.array([math.cos(math.radians(az_deg)), math.sin(math.radians(az_deg))])
    ku = kx*u[0] + ky*u[1]
    k2 = kx**2 + ky**2
    F = np.fft.fft2(dI)
    h = np.real(np.fft.ifft2(F*(-1j*ku)/(ku**2 + lam*k2 + eps)))
    return h

def hillshade(h, az_deg, elev_deg=30.0, zf=1.0):
    gy, gx = np.gradient(h*zf)
    a = math.radians(az_deg); e = math.radians(elev_deg)
    L = np.array([math.cos(a)*math.cos(e), math.sin(a)*math.cos(e), math.sin(e)])
    n = np.dstack([-gx, -gy, np.ones_like(h)])
    n /= np.linalg.norm(n, axis=2, keepdims=True)
    return np.clip(n[:, :, 0]*L[0] + n[:, :, 1]*L[1] + n[:, :, 2]*L[2], 0, 1)

def ring_light(h, n=24, elev_deg=25.0, zf=1.0):
    """Relight the RECOVERED SURFACE from all round - the 'emboss from 360' idea,
    but applied to height, where it means something."""
    acc = np.zeros_like(h, np.float32)
    for i in range(n):
        acc += hillshade(h, 360.0*i/n, elev_deg, zf)
    return acc/n

def openness(h, radius=14, n=16):
    """Negative openness: how enclosed each point is.  Concave (incised) points
    score high.  Direction-free, the standard relief-visualisation operator."""
    H, W = h.shape
    acc = np.zeros((H, W), np.float32)
    for i in range(n):
        a = 2*math.pi*i/n
        dx, dy = math.cos(a), math.sin(a)
        best = np.full((H, W), -1e9, np.float32)
        for r in range(2, radius+1):
            M = np.float32([[1, 0, -dx*r], [0, 1, -dy*r]])
            sh = cv2.warpAffine(h, M, (W, H), flags=cv2.INTER_LINEAR,
                                borderMode=cv2.BORDER_REPLICATE)
            best = np.maximum(best, (sh-h)/r)
        acc += np.arctan(best)
    return acc/n

def local_relief(h, sigma=18.0):
    return h - cv2.GaussianBlur(h, (0, 0), sigma)
