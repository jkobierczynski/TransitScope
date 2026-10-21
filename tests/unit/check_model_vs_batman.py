# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare transit_model.cpp with batman over a grid (python3 -I check_model_vs_batman.py build/model_dump)."""
import sys, subprocess, numpy as np, batman
exe = sys.argv[1]; rng = np.random.default_rng(3)
rows, ref = [], []
cases = [(rp, a, inc, u1, u2) for rp in (0.01, 0.05, 0.1, 0.155, 0.3, 0.5, 0.7, 1.2)
         for a, inc in ((8.86, 85.71), (8.86, 90.0), (5.0, 80.0), (15.0, 88.0)) for u1, u2 in ((0.5, 0.19), (0.0, 0.0), (0.9, -0.1), (0.2, 0.6))]
for rp, a, inc, u1, u2 in cases:
    t = np.linspace(-0.12, 0.12, 401)
    bp = batman.TransitParams(); bp.t0, bp.per, bp.rp, bp.a, bp.inc, bp.ecc, bp.w = 0., 2.21857567, rp, a, inc, 0., 90.
    bp.u, bp.limb_dark = [u1, u2], "quadratic"
    ref.append(batman.TransitModel(bp, t).light_curve(bp))
    rows += [f"0 2.21857567 {rp} {a} {inc} {u1} {u2} {x:.17g}" for x in t]
out = subprocess.run([exe], input="\n".join(rows), capture_output=True, text=True).stdout.split()
mine = np.array(out, float); ref = np.concatenate(ref)
d = np.abs(mine - ref); i = np.argmax(d)
print(f"{len(cases)} parameter sets x 401 times: max |C++ - batman| = {d.max():.2e} (case {i // 401}: {cases[i // 401]}), median {np.median(d):.1e}")
sys.exit(0 if d.max() < 1e-6 else 1)
