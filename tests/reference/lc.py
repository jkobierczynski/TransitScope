# SPDX-License-Identifier: GPL-3.0-or-later
"""Reference reduction, step 2: light curve, BJD_TDB timing, ensemble selection.

Usage: python3 -I lc.py phot.npz lc.npz
"""
import sys, json
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

PHOT, OUT = sys.argv[1:3]
d = np.load(PHOT)
flux, radii, stars, peak, bg = d["flux"], d["radii"], d["stars"], d["peak"], d["bg"]
nf, ns, nr = flux.shape

# ---- timing: DATE-OBS = exposure start (UTC, PC clock) -> mid-exposure BJD_TDB
SITE = EarthLocation(lat=50.93 * u.deg, lon=5.40 * u.deg, height=50 * u.m)  # Diepenbeek
TARGET = SkyCoord("20h00m43.713s +22d42m39.07s")
t_start = Time(list(d["date"]), scale="utc", format="isot", location=SITE)
t_mid = t_start + d["exp"] / 2 * u.s
bjd = (t_mid.tdb + t_mid.light_travel_time(TARGET, location=SITE)).jd
jd_utc_start = t_start.jd
alt = TARGET.transform_to(AltAz(obstime=t_mid, location=SITE)).alt.deg
airmass = 1 / np.sin(np.radians(alt))

# ---- ephemeris for masking (Agol et al. 2010)
T0, P = 2454279.436714, 2.21857567
n = np.round((bjd.mean() - T0) / P)
tc_pred = T0 + n * P
T14 = 0.07527  # days (~108 min)
oot = np.abs(bjd - tc_pred) > T14 / 2 + 0.003

# ---- frame quality flags
TGT = 0
temp = d["temp"]
good = np.isfinite(flux[:, TGT, :]).all(1)
good &= np.abs(temp + 5.0) < 0.35          # cooler regulated
good &= peak[:, TGT] < 45000               # linear regime
flags = dict(cooler=int((np.abs(temp + 5.0) >= 0.35).sum()), saturated=int((peak[:, TGT] >= 45000).sum()))

# ---- comparison candidates: valid in every frame, reasonably bright, away from M27 (x > 1150 in ref frame)
valid_all = np.isfinite(flux[:, :, 2]).all(0)
rel = np.nanmedian(flux[:, :, 2], 0) / np.nanmedian(flux[:, TGT, 2])
cand = [s for s in range(ns) if s != TGT and valid_all[s] and rel[s] > 0.03 and stars[s, 0] < 1150]


def detrend_rms(y, t, m, bin_n=10):
    """RMS of y[m] about a quadratic in time, unbinned and binned (~2.2 min)."""
    c = np.polyfit(t[m] - t[m].mean(), y[m], 2)
    r = y[m] - np.polyval(c, t[m] - t[m].mean())
    k = len(r) // bin_n
    rb = r[: k * bin_n].reshape(k, bin_n).mean(1)
    return r.std(), rb.std()


def ensemble_lc(ri, comps, m):
    F = flux[:, comps, ri]
    Fn = F / np.nanmedian(F[m], 0)
    # weights from each comp's own OOT scatter against the others
    w = []
    for j in range(len(comps)):
        others = np.delete(np.arange(len(comps)), j)
        ref = Fn[:, others].mean(1) if len(others) else np.ones(nf)
        w.append(1 / detrend_rms(Fn[:, j] / ref, bjd, m)[0] ** 2)
    w = np.array(w) / np.sum(w)
    refl = (Fn * w).sum(1)
    lc = flux[:, TGT, ri] / refl
    return lc / np.median(lc[m]), w


import batman
from scipy.optimize import least_squares
_bp = batman.TransitParams()
_bp.t0, _bp.per, _bp.rp, _bp.a, _bp.inc, _bp.ecc, _bp.w = tc_pred, P, 0.155, 8.84, 85.7, 0., 90.
_bp.u, _bp.limb_dark = [0.50, 0.19], "quadratic"
_bm = batman.TransitModel(_bp, bjd)


def model_rms(lc, m):
    """Residual RMS after fitting a fixed-geometry transit (free Tc, Rp/Rs) times a quadratic baseline.
    Uses in-transit points through a model, so comparison choice cannot 'shape' the transit for free."""
    tt = bjd - bjd.mean()
    def res(p):
        _bp.t0, _bp.rp = p[0], p[1]
        return (lc - _bm.light_curve(_bp) * (p[2] + p[3] * tt + p[4] * tt ** 2))[m]
    s = least_squares(res, [tc_pred, 0.155, 1, 0, 0], x_scale=[1e-3, 1e-2, 1e-3, 1e-2, 1e-1])
    return res(s.x).std()


best = None
m = good
log = []
for ri, R in enumerate(radii):
    comps = list(cand)
    lc, _ = ensemble_lc(ri, comps, good & oot)
    cur = model_rms(lc, m)
    while len(comps) > 3:
        trial = []
        for c in comps:
            cc = [x for x in comps if x != c]
            trial.append((model_rms(ensemble_lc(ri, cc, good & oot)[0], m), c))
        bmin, cdrop = min(trial)
        if bmin < cur * 0.998: comps.remove(cdrop); cur = bmin
        else: break
    lc, _ = ensemble_lc(ri, comps, good & oot)
    rms_oot, _ = detrend_rms(lc, bjd, good & oot)
    log.append(dict(R=int(R), n=len(comps), model_rms=round(cur, 5), oot_rms=round(rms_oot, 5)))
    if best is None or cur < best[0]:
        best = (cur, ri, comps)
m = good & oot

_, ri, comps = best
lc, w = ensemble_lc(ri, comps, m)
# per-point error estimate: photon + read + sky + Young scintillation (D=20 cm)
gain, rn_e = 0.26, 4.0
Ft = flux[:, TGT, ri] * gain
npx = d["npx"][:, TGT, ri]
sig_ph = np.sqrt(Ft + npx * (bg[:, TGT] * gain + rn_e ** 2)) / Ft
sig_sc = 0.09 * 20 ** (-2 / 3) * airmass ** 1.75 * np.exp(-50 / 8000) / np.sqrt(2 * d["exp"])
sig = np.sqrt(sig_ph ** 2 + sig_sc ** 2)

np.savez(OUT, bjd=bjd, jd_utc_start=jd_utc_start, flux=lc, sig=sig, good=good, oot=oot,
         airmass=airmass, r90=d["r90"][:, TGT], xy=d["xy"][:, TGT], temp=temp, tc_pred=tc_pred,
         epoch=n, radius=radii[ri], comps=np.array(comps), weights=w)
summ = dict(radius=int(radii[ri]), comps=[int(c) for c in comps], weights=[round(float(x), 3) for x in w],
            comp_xy=[[int(v) for v in stars[c]] for c in comps], n_good=int(good.sum()), flags=flags,
            tc_pred=float(tc_pred), epoch=int(n), aperture_scan=log,
            sig_photon_med=float(np.median(sig_ph)), sig_scint_range=[float(sig_sc.min()), float(sig_sc.max())])
print(json.dumps(summ, indent=1, default=float))
