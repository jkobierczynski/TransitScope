#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Independent numpy re-implementation of `transitscope ld` (table interpolation, wavelength
interpolation, photon-weighted passband average) and comparison with its JSON output.

usage: check_ld.py ld.json [DATA_DIR]      exit 1 if |delta| > 2e-4 in any coefficient
"""
import json
import sys

import numpy as np

res = json.load(open(sys.argv[1]))
data = sys.argv[2] if len(sys.argv) > 2 else "data"
star = res["star"]
LAM = {"B": 438, "V": 545, "R": 641, "I": 798, "J": 1220, "Sloang": 483, "Sloanr": 626, "Sloani": 767, "Sloanz": 910}

tab = {}
for line in open(f"{data}/ld/claret2011_atlas_quad.csv"):
    if line.startswith("#") or line.startswith("method"):
        continue
    m, b, feh, t, g, a, bb = line.strip().split(",")
    tab[(m, b, float(feh), float(t), float(g))] = (float(a), float(bb))
fehs, teffs, loggs = (sorted({k[i] for k in tab}) for i in (2, 3, 4))


def trilinear(m, band):
    out = np.zeros(2)
    idx = []
    for v, g in ((star["feh"], fehs), (star["teff"], teffs), (star["logg"], loggs)):
        i = min(max(np.searchsorted(g, v, side="right") - 1, 0), len(g) - 2)
        idx.append((i, (v - g[i]) / (g[i + 1] - g[i]), g))
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = np.prod([f if d else 1 - f for d, (_, f, _) in zip((di, dj, dk), idx)])
                if w:
                    key = (m, band) + tuple(g[i + d] for d, (i, _, g) in zip((di, dj, dk), idx))
                    out += w * np.array(tab[key])
    return out


def curve(path):
    d = np.array([[float(x) for x in l.split(",")[:2]] for l in open(path) if l[:1].isdigit()])
    return d[:, 0], d[:, 1]


fl, ft = curve(f"{data}/passbands/astrodon_xop_bb.csv")
ql, qr = curve(f"{data}/passbands/sony_icx285al_qe.csv")
lo, hi = max(fl[0], ql[0]), min(fl[-1], ql[-1])
grid = np.arange(lo, hi + 1e-9, 1.0)
lm = grid * 1e-9
w = np.interp(grid, fl, ft, 0, 0) * np.interp(grid, ql, qr, 0, 0) * lm**-4 / np.expm1(1.438776877e-2 / (lm * star["teff"]))

bad = 0
for m in "LF":
    bands = sorted(LAM, key=LAM.get)
    L = np.array([LAM[b] for b in bands])
    c = np.array([trilinear(m, b) for b in bands])
    u1 = (w * np.interp(grid, L, c[:, 0])).sum() / w.sum()
    u2 = (w * np.interp(grid, L, c[:, 1])).sum() / w.sum()
    d1, d2 = u1 - res[f"u1_{m}"], u2 - res[f"u2_{m}"]
    print(f"method {m}: numpy u1 {u1:.5f} u2 {u2:.5f}   C++ - numpy: {-d1:+.1e} {-d2:+.1e}")
    bad += abs(d1) > 2e-4 or abs(d2) > 2e-4
lam_eff = (w * grid).sum() / w.sum()
print(f"lambda_eff numpy {lam_eff:.2f} nm, C++ {res['lambda_eff_nm']:.2f} nm")
sys.exit(1 if bad or abs(lam_eff - res["lambda_eff_nm"]) > 0.5 else 0)
