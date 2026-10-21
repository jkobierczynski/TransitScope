# SPDX-License-Identifier: GPL-3.0-or-later
"""Reference reduction, step 3: MCMC transit fit.

Usage: python3 -I fit.py lc.npz <free|constrained> out_prefix
  free        : Tc, Rp/Rs, a/Rs, inc free (uniform priors), LD with Gaussian priors
  constrained : a/Rs and inc with literature Gaussian priors (Agol et al. 2010) -> robust Tc and Rp/Rs
"""
import sys, json
import numpy as np, emcee, batman

L = np.load(sys.argv[1]); MODE = sys.argv[2]; PFX = sys.argv[3]
g = L["good"]
t, f, s0, X = L["bjd"][g], L["flux"][g], L["sig"][g], L["airmass"][g]
P, tc0 = 2.21857567, float(L["tc_pred"])
tt = t - t.mean()
bp = batman.TransitParams()
bp.t0, bp.per, bp.rp, bp.a, bp.inc, bp.ecc, bp.w, bp.u, bp.limb_dark = tc0, P, .155, 8.86, 85.7, 0, 90, [.5, .19], "quadratic"
bm = batman.TransitModel(bp, t)
LD_PRIOR = (0.50, 0.19, 0.10)        # approx. quadratic LD for K1V (Teff~5050 K) in a ~R-band passband (XOP-BB proxy)
A_PRIOR, I_PRIOR = (8.863, 0.05), (85.71, 0.05)
names = ["Tc", "rp", "a", "inc", "u1", "u2", "c0", "c1", "c2", "lnf"]


def model(p):
    bp.t0, bp.rp, bp.a, bp.inc, bp.u = p[0], p[1], p[2], p[3], [p[4], p[5]]
    return bm.light_curve(bp) * (p[6] + p[7] * tt + p[8] * tt ** 2)


def lnprob(p):
    Tc, rp, a, inc, u1, u2, c0, c1, c2, lnf = p
    if not (abs(Tc - tc0) < 0.03 and 0.05 < rp < 0.3 and 3 < a < 20 and 75 < inc <= 90 and -1 < lnf < 2): return -np.inf
    if not (0 < u1 + u2 < 1 and u1 > 0): return -np.inf
    if a * np.cos(np.radians(inc)) > 1 + rp: return -np.inf
    lp = -0.5 * ((u1 - LD_PRIOR[0]) ** 2 + (u2 - LD_PRIOR[1]) ** 2) / LD_PRIOR[2] ** 2
    if MODE == "constrained":
        lp += -0.5 * ((a - A_PRIOR[0]) / A_PRIOR[1]) ** 2 - 0.5 * ((inc - I_PRIOR[0]) / I_PRIOR[1]) ** 2
    s = s0 * np.exp(lnf)
    r = f - model(p)
    return lp - 0.5 * np.sum(r ** 2 / s ** 2 + 2 * np.log(s))


p0 = np.array([tc0, .155, 8.86, 85.7, .5, .19, 1, 0, 0, 0.])
nw, nd = 40, len(p0)
rng = np.random.default_rng(1)
pos = p0 + rng.normal(0, 1, (nw, nd)) * np.array([5e-4, 2e-3, .05, .05, .02, .02, 1e-4, 1e-3, 1e-2, .01])
sam = emcee.EnsembleSampler(nw, nd, lnprob)
sam.run_mcmc(pos, int(__import__("os").environ.get("NSTEP", 6000)), progress=False)
ch = sam.get_chain(discard=2000, thin=10, flat=True)
best = ch[np.argmax(sam.get_log_prob(discard=2000, thin=10, flat=True))]

# red-noise beta (Winn et al. 2008 time-averaging) on residuals
res = f - model(best)
betas = []
for n in (5, 10, 15, 20, 30):
    k = len(res) // n
    sb = res[: k * n].reshape(k, n).mean(1).std()
    betas.append(sb / (res.std() / np.sqrt(n) * np.sqrt(k / (k - 1))))
beta = max(1.0, float(np.median(betas)))

out = {}
for i, nm in enumerate(names):
    q = np.percentile(ch[:, i], [16, 50, 84])
    out[nm] = dict(med=float(q[1]), lo=float(q[1] - q[0]), hi=float(q[2] - q[1]))
out["beta"] = beta
out["Tc_err_inflated_min"] = float(np.mean([out["Tc"]["lo"], out["Tc"]["hi"]]) * beta * 1440)
out["rms_resid"] = float(res.std())
out["depth_pct"] = float(np.median(ch[:, 1] ** 2) * 100)
out["b"] = float(np.median(ch[:, 2] * np.cos(np.radians(ch[:, 3]))))
out["acceptance"] = float(np.mean(sam.acceptance_fraction))
try: out["autocorr_max"] = float(np.max(sam.get_autocorr_time(discard=2000, quiet=True)))
except Exception: pass
out["O_C_min"] = float((out["Tc"]["med"] - tc0) * 1440)
np.save(PFX + "_chain.npy", ch); np.save(PFX + "_best.npy", best)
json.dump(out, open(PFX + ".json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ["Tc", "rp", "a", "inc", "u1", "beta", "Tc_err_inflated_min", "rms_resid", "depth_pct", "b", "O_C_min", "acceptance"]}, indent=1))
