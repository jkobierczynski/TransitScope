# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare a transitscope light-curve CSV with the Python reference lc.npz."""
import sys, numpy as np
L = np.load(sys.argv[2])
lines = [l for l in open(sys.argv[1]) if not l.startswith("#")]
d = np.genfromtxt(lines, delimiter=",", names=True)
g = L["good"]
ok = (d["good"].astype(bool) == g).all() and (d["oot"].astype(bool) == L["oot"]).all()
print("good/oot masks identical:", ok)
df = np.abs(d["flux"][g] - L["flux"][g]); ds = np.abs(d["sigma"][g] / L["sig"][g] - 1)
print(f"flux max |diff| {df.max()*1e3:.4f} mmag, rms {np.sqrt(np.mean(df**2))*1e3:.4f} mmag; sigma max rel diff {ds.max():.2e}")
sys.exit(0 if ok and df.max() < float(sys.argv[3]) else 1)
