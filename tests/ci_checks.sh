#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Checks that need no image data (CI on Linux, macOS and Windows/MSYS2). They start from the Python
# reference photometry in tests/reference (phot.npz), so the 2.4 GB of frames are not needed.
#
# usage: tests/ci_checks.sh BUILD_DIR [WORK_DIR]
#   needs python3 with numpy; with the batman-package module the transit model is checked too
set -euo pipefail
B=$(cd "$1" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
W=${2:-$B/ci-work}
mkdir -p "$W"
PY=${PYTHON:-python3}
step() { echo; echo "=== $*"; }

step "project file"
(cd "$W" && "$B/project_test")

step "Python reference photometry -> CSV"
"$PY" -I "$HERE/npz_to_phot_csv.py" "$HERE/reference/phot.npz" "$HERE/reference/lc.npz" "$W/phot_ref.csv"

step "ensemble light curve vs the Python reference (its 0.2 % selection rule)"
"$B/transitscope" lightcurve --phot "$W/phot_ref.csv" --T0 2454279.436714 --P 2.21857567 --T14 0.07527 \
  --target 0 --exclude-stars 10,32,36 --gain 0.26 --read-noise 4 --aperture-cm 20 --height 50 \
  --min-improvement 0.002 --out "$W/lc_ref.csv"
"$PY" -I "$HERE/compare_lightcurve.py" "$W/lc_ref.csv" "$HERE/reference/lc.npz" 1e-9

step "MCMC fit vs emcee (constrained, 16000 steps)"
"$B/transitscope" fit --lc "$W/lc_ref.csv" --mode constrained --out "$W/fitC.json"
"$PY" -I "$HERE/compare_fit.py" "$W/fitC.json" "$HERE/reference/fitC.json"

step "limb darkening vs an independent numpy calculation"
"$B/transitscope" ld --teff 5050 --logg 4.59 --feh -0.03 --teff-err 50 --logg-err 0.02 --feh-err 0.08 --out "$W/ld.json"
"$PY" -I "$HERE/ld/check_ld.py" "$W/ld.json" "$ROOT/data"

step "fit API: progress, cancelling, in-memory data, thread count"
(cd "$W" && "$B/fit_api_test" "$W/phot_ref.csv")

step "comparison-star selection stability"
"$B/comp_stability_test" "$W/phot_ref.csv"

if "$PY" -I -c "import batman" 2>/dev/null; then   # -I as below: user site-packages do not count
  step "transit model vs batman"
  "$PY" -I "$HERE/unit/check_model_vs_batman.py" "$B/model_dump"
else
  step "transit model vs batman: skipped (pip install batman-package)"
fi

echo
echo "all checks passed"
