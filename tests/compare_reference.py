# SPDX-License-Identifier: GPL-3.0-or-later
"""Regression check: C++ transitscope output vs the Python reference reduction.

python3 -I compare_reference.py phot_cpp.csv phot.npz lc.npz
Exit code 1 if any tolerance is exceeded.

Note on comparison stars: the Python reference started each centroid from an FFT shift quantised
to 4 px, the C++ code from the tracked sub-pixel offset. For faint stars the 4-iteration centroid
does not fully converge, so their positions legitimately differ by up to a few pixels. The checks
therefore compare the target star exactly, and the comparison stars at the level that matters:
the final ensemble light curve (same stars, weights and aperture as lc.npz).
"""
import sys, csv
import numpy as np

cpp_csv, ref_npz, lc_npz = sys.argv[1:4]
R = np.load(ref_npz); L = np.load(lc_npz)
nf, ns, nr = R["flux"].shape
radii = list(R["radii"])
flux = np.full((nf, ns, nr), np.nan); xy = np.full((nf, ns, 2), np.nan)
bjd = np.full(nf, np.nan); am = np.full(nf, np.nan); r90 = np.full((nf, ns), np.nan); bg = np.full((nf, ns), np.nan)
with open(cpp_csv) as fh:
    for row in csv.DictReader(fh):
        i, s = int(row["frame"]), int(row["star"])
        bjd[i], am[i] = float(row["bjd_tdb"]), float(row["airmass"])
        if row["valid"] != "1": continue
        xy[i, s] = float(row["x"]), float(row["y"]); bg[i, s] = float(row["bg"]); r90[i, s] = float(row["r90"])
        flux[i, s] = [float(row[f"flux_r{r:g}"]) for r in radii]

ok = True
def report(name, val, tol, unit=""):
    global ok
    flag = "OK  " if val <= tol else "FAIL"
    ok &= bool(val <= tol)
    print(f"{flag} {name:<52s} {val:.3g} {unit} (tol {tol:g})")

report("BJD_TDB max |C++ (ERFA) - astropy|", np.nanmax(np.abs(bjd - L["bjd"])) * 86400, 0.05, "s")
report("airmass max |diff|", np.nanmax(np.abs(am - L["airmass"])), 1e-3)
d0 = np.hypot(*(xy[:, 0] - R["xy"][:, 0]).T)
report("target centroid max offset", np.nanmax(d0), 0.01, "px")
rel0 = np.abs(flux[:, 0] - R["flux"][:, 0]) / R["flux"][:, 0]
# hard-edged apertures: a 1e-5 px centroid difference can move one boundary pixel in/out (~3e-4)
report("target flux max relative diff (all radii)", np.nanmax(rel0), 1e-3)
report("target flux median relative diff", np.nanmedian(rel0), 1e-6)
report("target background max |diff|", np.nanmax(np.abs(bg[:, 0] - R["bg"][:, 0])), 0.05, "ADU")
report("target r90 mismatching frames", float(np.sum(r90[:, 0] != R["r90"][:, 0])), 0)

# ensemble light curve with the reference's own aperture, comparison stars and weights
ri = radii.index(L["radius"]); comps = list(L["comps"]); w = L["weights"]; m = L["good"] & L["oot"]
F = flux[:, comps, ri]; Fn = F / np.nanmedian(F[m], 0)
lc = flux[:, 0, ri] / (Fn * w).sum(1); lc /= np.median(lc[m])
g = L["good"]
dlc = np.abs(lc[g] - L["flux"][g])
report("ensemble light curve max |diff|", np.nanmax(dlc) * 1e3, 0.5, "mmag")
report("ensemble light curve rms diff", np.sqrt(np.nanmean(dlc ** 2)) * 1e3, 0.1, "mmag")
print(f"     (light-curve scatter itself is ~6.4 mmag per frame)")
sys.exit(0 if ok else 1)
