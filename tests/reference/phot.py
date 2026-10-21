# SPDX-License-Identifier: GPL-3.0-or-later
"""Reference reduction, step 1: calibration + multi-aperture photometry.

Usage: python3 -I phot.py <lights_dir> <calib_dir> <out.npz>
"""
import sys, glob, os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage
from multiprocessing import Pool

LIGHTS, CAL, OUT = sys.argv[1:4]
RADII = np.array([8, 10, 12, 14, 16, 18, 20])
ANN_IN, ANN_OUT = 26, 40
SAT_ADU = 45000.0          # raw ADU; ~70% of the 16-bit range

DARK = fits.getdata(os.path.join(CAL, "master_dark_10s.fits")).astype(np.float32)
FLAT = fits.getdata(os.path.join(CAL, "master_flat_XOPBB.fits")).astype(np.float32)
HOT = np.load(os.path.join(CAL, "hotmask.npy"))
FLAT = np.where(FLAT > 0.5, FLAT, 1.0)
files = sorted(glob.glob(os.path.join(LIGHTS, "HD189733b_*.fit")))
NY, NX = DARK.shape
YY, XX = np.mgrid[0:NY, 0:NX]


def calibrate(path):
    raw = fits.getdata(path).astype(np.float32)
    img = (raw - DARK) / FLAT
    if HOT.any():  # replace hot pixels with local median
        med = ndimage.median_filter(img, size=3)
        img = np.where(HOT, med, img)
    return raw, img


def coarse_shift(ref_small, img):
    small = ndimage.zoom(img - np.median(img), 0.25, order=1)
    F = np.fft.fft2(ref_small) * np.conj(np.fft.fft2(small))
    cc = np.fft.ifft2(F).real
    y, x = np.unravel_index(cc.argmax(), cc.shape)
    if y > cc.shape[0] / 2: y -= cc.shape[0]
    if x > cc.shape[1] / 2: x -= cc.shape[1]
    return -4.0 * x, -4.0 * y


def centroid(img, x, y, r=16, it=4):
    for _ in range(it):
        xi, yi = int(round(x)), int(round(y))
        x0, x1, y0, y1 = max(xi - 45, 0), min(xi + 46, NX), max(yi - 45, 0), min(yi + 46, NY)
        cut = img[y0:y1, x0:x1]; yy, xx = YY[y0:y1, x0:x1], XX[y0:y1, x0:x1]
        rr = np.hypot(xx - x, yy - y)
        bg = np.median(cut[(rr > ANN_IN) & (rr < ANN_OUT)])
        w = np.clip(cut - bg, 0, None) * (rr < r)
        if w.sum() <= 0: break
        x, y = (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    return x, y


def measure(img, raw, x, y):
    xi, yi = int(round(x)), int(round(y))
    x0, x1, y0, y1 = max(xi - 45, 0), min(xi + 46, NX), max(yi - 45, 0), min(yi + 46, NY)
    cut = img[y0:y1, x0:x1]; rc = raw[y0:y1, x0:x1]
    rr = np.hypot(XX[y0:y1, x0:x1] - x, YY[y0:y1, x0:x1] - y)
    _, bg, bgsd = sigma_clipped_stats(cut[(rr > ANN_IN) & (rr < ANN_OUT)], sigma=3)
    sub = cut - bg
    fl = np.array([sub[rr < R].sum() for R in RADII])
    npx = np.array([(rr < R).sum() for R in RADII])
    tot = sub[rr < 22].sum()
    rs = np.arange(1, 22)
    ee = np.array([sub[rr < k].sum() for k in rs]) / tot if tot > 0 else np.ones(len(rs))
    r90 = float(rs[min(np.searchsorted(ee, 0.9), len(rs) - 1)])
    return fl, npx, bg, bgsd, float(rc[rr < 12].max()), r90


# --- reference frame & star list ---
REF_IDX = 99  # frame 100: cooler settled, early in the run
_, ref = calibrate(files[REF_IDX])
ref_small = ndimage.zoom(ref - np.median(ref), 0.25, order=1)


def find_stars(img, nmax=40):
    sm = ndimage.gaussian_filter(img, 4)
    mx = ndimage.maximum_filter(sm, size=31)
    _, med, sd = sigma_clipped_stats(sm)
    pk = (sm == mx) & (sm > med + 25 * sd)
    ys, xs = np.nonzero(pk)
    order = np.argsort(sm[ys, xs])[::-1]
    return [(float(xs[i]), float(ys[i])) for i in order[:nmax]]


def init():
    pass


def process(i):
    path = files[i]
    raw, img = calibrate(path)
    h = fits.getheader(path)
    dx, dy = coarse_shift(ref_small, img)
    res = []
    for (sx, sy) in STARS:
        x, y = sx + dx, sy + dy
        if not (60 < x < NX - 60 and 60 < y < NY - 60):
            res.append(None); continue
        x, y = centroid(img, x, y)
        res.append((x, y) + measure(img, raw, x, y))
    return dict(i=i, date=h["DATE-OBS"], exp=float(h["EXPTIME"]), temp=float(h.get("CCD-TEMP", np.nan)),
                dx=dx, dy=dy, res=res)


cands = find_stars(ref)
# keep stars that are not near the edge in the reference frame
STARS = [(x, y) for x, y in cands if 70 < x < NX - 70 and 70 < y < NY - 70]
print("reference stars:", len(STARS))

if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        out = p.map(process, range(len(files)), chunksize=8)
    ns, nr = len(STARS), len(RADII)
    nf = len(out)
    flux = np.full((nf, ns, nr), np.nan); npx = np.full((nf, ns, nr), np.nan)
    xy = np.full((nf, ns, 2), np.nan); bg = np.full((nf, ns), np.nan); bgsd = np.full((nf, ns), np.nan)
    peak = np.full((nf, ns), np.nan); r90 = np.full((nf, ns), np.nan)
    for o in out:
        i = o["i"]
        for s, r in enumerate(o["res"]):
            if r is None: continue
            xy[i, s] = r[0], r[1]; flux[i, s] = r[2]; npx[i, s] = r[3]
            bg[i, s], bgsd[i, s], peak[i, s], r90[i, s] = r[4], r[5], r[6], r[7]
    np.savez(OUT, flux=flux, npx=npx, xy=xy, bg=bg, bgsd=bgsd, peak=peak, r90=r90, radii=RADII,
             stars=np.array(STARS), date=np.array([o["date"] for o in out]),
             exp=np.array([o["exp"] for o in out]), temp=np.array([o["temp"] for o in out]),
             shift=np.array([(o["dx"], o["dy"]) for o in out]),
             files=np.array([os.path.basename(f) for f in files]))
    print("saved", OUT, flux.shape)
