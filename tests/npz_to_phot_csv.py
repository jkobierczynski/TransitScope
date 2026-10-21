# SPDX-License-Identifier: GPL-3.0-or-later
"""Write the Python reference photometry (phot.npz + lc.npz timing) in the transitscope reduce CSV format,
so the C++ light-curve and fit stages can be tested on exactly the same input as the Python ones."""
import sys, numpy as np
R = np.load(sys.argv[1]); L = np.load(sys.argv[2]); out = sys.argv[3]
nf, ns, nr = R["flux"].shape; radii = R["radii"]
with open(out, "w") as o:
    o.write("frame,file,date_obs,exptime,ccd_temp,bjd_tdb,airmass,star,valid,x,y,bg,bg_sd,peak_raw,r90,"
            + ",".join(f"flux_r{r:g}" for r in radii) + "," + ",".join(f"npix_r{r:g}" for r in radii) + "\n")
    for i in range(nf):
        for s in range(ns):
            v = int(np.isfinite(R["flux"][i, s, 0]))
            vals = [i, R["files"][i], R["date"][i], R["exp"][i], R["temp"][i], repr(float(L["bjd"][i])), repr(float(L["airmass"][i])),
                    s, v, R["xy"][i, s, 0], R["xy"][i, s, 1], R["bg"][i, s], R["bgsd"][i, s], R["peak"][i, s], R["r90"][i, s]]
            o.write(",".join(str(x) for x in vals) + "," + ",".join(repr(float(x)) for x in R["flux"][i, s]) + ","
                    + ",".join(str(int(x)) if np.isfinite(x) else "0" for x in R["npx"][i, s]) + "\n")
