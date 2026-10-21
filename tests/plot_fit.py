# SPDX-License-Identifier: GPL-3.0-or-later
"""Plot a transitscope fit: python3 -I plot_fit.py model.csv fit.json out.png [title]"""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
m = np.genfromtxt(sys.argv[1], delimiter=",", names=True); r = json.load(open(sys.argv[2]))
title = sys.argv[4] if len(sys.argv) > 4 else "transitscope fit"
tc = r["Tc"]["med"]; x = (m["bjd_tdb"] - tc) * 24; fd, mod = m["flux_detrended"], m["model_transit"]
res = fd - mod; nb = 22; k = len(x) // nb
xb = x[:k*nb].reshape(k, nb).mean(1); fb = fd[:k*nb].reshape(k, nb).mean(1)
eb = fd[:k*nb].reshape(k, nb).std(1) / np.sqrt(nb); rb = res[:k*nb].reshape(k, nb).mean(1)
INK, INK2, GRID, SURF, S1, S2 = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.08), facecolor=SURF)
for a in (a1, a2): a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.8); a.set_axisbelow(True)
a1.plot(x, fd, ".", ms=2.5, color="#b9b8b2", label=f"10 s frames ({len(x)})")
a1.errorbar(xb, fb, eb, fmt="o", ms=5, color=S1, ecolor=S1, elinewidth=1, mec=SURF, mew=1, label="5-min bins")
o = np.argsort(x); a1.plot(x[o], mod[o], color=S2, lw=2, label="Transit model")
a1.set_ylabel("Relative flux (detrended)"); a1.set_ylim(0.955, 1.02)
a1.legend(frameon=False, loc="lower left", fontsize=9, labelcolor=INK2)
a1.set_title(title, loc="left", color=INK, fontsize=11)
e = 0.5 * (r["Tc"]["lo"] + r["Tc"]["hi"])
a1.text(0.99, 0.04, f"Tc = {tc:.5f} ± {e:.5f} BJD_TDB\nRp/Rs = {r['rp']['med']:.4f} ± {0.5*(r['rp']['lo']+r['rp']['hi']):.4f}\n"
        f"O–C = {r['O_C_min']:+.1f} ± {r['Tc_err_inflated_min']:.1f} min   rms = {r['rms_resid']*1000:.1f} mmag/frame",
        transform=a1.transAxes, ha="right", va="bottom", color=INK, fontsize=9)
a2.plot(x, res * 1000, ".", ms=2.5, color="#b9b8b2"); a2.plot(xb, rb * 1000, "o", ms=5, color=S1, mec=SURF, mew=1)
a2.axhline(0, color=INK2, lw=1); a2.set_ylabel("Residual (ppt)"); a2.set_xlabel("Hours from mid-transit"); a2.set_ylim(-25, 25)
fig.savefig(sys.argv[3], dpi=150, bbox_inches="tight", facecolor=SURF)
