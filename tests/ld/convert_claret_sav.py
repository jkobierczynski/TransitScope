#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Convert the Claret & Bloemen (2011) quadratic limb-darkening tables, as distributed in
EXOFASTv2's quadld/ directory (IDL .sav, ATLAS models, microturbulence 2 km/s), into the CSV
read by `transitscope ld`.  Data-preparation only; not needed at runtime.

usage: convert_claret_sav.py EXOFAST_QUADLD_DIR OUT.csv
"""
import sys
import numpy as np
from scipy.io import readsav

BANDS = ["B", "V", "R", "I", "J", "Sloang", "Sloanr", "Sloani", "Sloanz"]
TEFF_MAX = 10000.0

src, out = sys.argv[1], sys.argv[2]
rows = 0
with open(out, "w") as f:
    f.write("# Quadratic limb-darkening coefficients I(mu)/I(1) = 1 - a(1-mu) - b(1-mu)^2\n"
            "# Claret A., Bloemen S., 2011, A&A 529, A75 (VizieR J/A+A/529/A75), ATLAS model atmospheres,\n"
            "# microturbulent velocity 2 km/s; method L = least-squares fit, F = flux-conservation fit.\n"
            "# Converted from EXOFASTv2 quadld/ATLAS.<method>.2.<band>.sav (Eastman et al. 2013, 2019);\n"
            f"# Teff <= {TEFF_MAX:.0f} K only; grid cells missing in the original tables are omitted.\n"
            "method,band,feh,teff,logg,a,b\n")
    for method in "LF":
        for band in BANDS:
            s = readsav(f"{src}/ATLAS.{method}.2.{band}.sav")
            fehs, teffs, loggs = s["fehs"], s["teffs"], s["loggs"]
            A, B = s["quadlda"], s["quadldb"]
            assert A.shape == (len(fehs), len(teffs), len(loggs)), A.shape
            for i, feh in enumerate(fehs):
                for j, teff in enumerate(teffs):
                    if teff > TEFF_MAX:
                        continue
                    for k, logg in enumerate(loggs):
                        a, b = A[i, j, k], B[i, j, k]
                        if not (np.isfinite(a) and np.isfinite(b)):
                            continue
                        f.write(f"{method},{band},{feh:g},{teff:g},{logg:g},{a:.4f},{b:.4f}\n")
                        rows += 1
print(f"{rows} rows -> {out}")
