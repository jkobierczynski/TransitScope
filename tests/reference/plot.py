# SPDX-License-Identifier: GPL-3.0-or-later
import numpy as np, json, batman, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L = np.load("lc.npz"); g = L["good"]
t, f = L["bjd"][g], L["flux"][g]
best = np.load("fitC_best.npy"); res = json.load(open("fitC.json"))
tt = t - t.mean()
bp = batman.TransitParams()
bp.t0, bp.per, bp.rp, bp.a, bp.inc = best[0], 2.21857567, best[1], best[2], best[3]
bp.ecc, bp.w, bp.u, bp.limb_dark = 0, 90, [best[4], best[5]], "quadratic"
base = best[6] + best[7] * tt + best[8] * tt ** 2
fd = f / base                                   # detrended
tf = np.linspace(t.min(), t.max(), 2000)
mod = batman.TransitModel(bp, tf).light_curve(bp)
resid = fd - batman.TransitModel(bp, t).light_curve(bp)
x = (t - best[0]) * 24; xf = (tf - best[0]) * 24
nb = 22                                          # ~5 min bins
k = len(t) // nb
xb = x[:k * nb].reshape(k, nb).mean(1); fb = fd[:k * nb].reshape(k, nb).mean(1)
eb = fd[:k * nb].reshape(k, nb).std(1) / np.sqrt(nb); rb = resid[:k * nb].reshape(k, nb).mean(1)

INK, INK2, GRID, SURF, S1, S2 = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.08),
                             facecolor=SURF)
for a in (a1, a2):
    a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.8); a.set_axisbelow(True)
a1.plot(x, fd, ".", ms=2.5, color="#b9b8b2", label="10 s frames (706)")
a1.errorbar(xb, fb, eb, fmt="o", ms=5, color=S1, ecolor=S1, elinewidth=1, mec=SURF, mew=1, label="5-min bins")
a1.plot(xf, mod, color=S2, lw=2, label="Transit model")
a1.set_ylabel("Relative flux (detrended)")
a1.set_ylim(0.955, 1.02)
a1.legend(frameon=False, loc="lower left", fontsize=9, labelcolor=INK2)
tc = res["Tc"]
a1.set_title("HD 189733 b — 27 Sep 2014, Diepenbeek · 200 mm f/4 + Atik 314L + Astrodon XOP-BB", loc="left", color=INK, fontsize=11)
a1.text(0.99, 0.04, f"Tc = {tc['med']:.5f} ± {0.5*(tc['lo']+tc['hi']):.5f} BJD_TDB\nRp/Rs = {res['rp']['med']:.4f} ± {0.5*(res['rp']['lo']+res['rp']['hi']):.4f}\n"
         f"O–C = {res['O_C_min']:+.1f} ± {res['Tc_err_inflated_min']:.1f} min   rms = {res['rms_resid']*1000:.1f} mmag/frame",
         transform=a1.transAxes, ha="right", va="bottom", color=INK, fontsize=9)
a2.plot(x, resid * 1000, ".", ms=2.5, color="#b9b8b2")
a2.plot(xb, rb * 1000, "o", ms=5, color=S1, mec=SURF, mew=1)
a2.axhline(0, color=INK2, lw=1)
a2.set_ylabel("Residual (ppt)"); a2.set_xlabel("Hours from mid-transit")
a2.set_ylim(-25, 25)
fig.savefig("hd189733b_2014_reanalysis.png", dpi=150, bbox_inches="tight", facecolor=SURF)
print("binned resid rms ppt", rb.std() * 1000)
