# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare a transitscope fit JSON with a reference fit JSON (fit.py output).
python3 -I compare_fit.py cpp_fit.json fitC.json [max_sigma=0.2]"""
import sys, json
A, B = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
lim = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
worst = 0.0
for k in ("Tc", "rp", "a", "inc", "u1", "u2", "lnf"):
    sig = 0.5 * (B[k]["lo"] + B[k]["hi"])
    z = (A[k]["med"] - B[k]["med"]) / sig
    err = 0.5 * (A[k]["lo"] + A[k]["hi"]) / sig - 1
    worst = max(worst, abs(z))
    print(f"{k:4s} {A[k]['med']:.6f} vs {B[k]['med']:.6f}  diff {z:+.3f} sigma, error bar ratio {err:+.1%}")
print(f"worst |diff| = {worst:.3f} sigma (limit {lim})")
sys.exit(0 if worst <= lim else 1)
