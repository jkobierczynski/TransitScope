#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Full HD 189733 b reduction with transitscope: masters -> reduce -> lightcurve -> ld -> fit.
#
# usage: tests/run_end_to_end.sh BUILD_DIR LIGHTS_DIR FLATS_DIR BIAS_DIR DARKS_DIR OUT_DIR [backend]
#
# Frames are found by name, case-insensitively, with extension .fit, .fits or .fts:
#   lights LIGHTS_DIR/<LIGHT_PREFIX>*   LIGHT_PREFIX (environment) defaults to HD189733b_
#   flats  FLATS_DIR/*flat*   (every FITS file in FLATS_DIR if none has "flat" in its name and
#                              FLATS_DIR is not the lights folder)
#   bias   BIAS_DIR/*bias*    (same fallback)
#   darks  DARKS_DIR/*dark*   (same fallback)
#
# THREADS (environment) sets the worker threads for reduce and fit; default 0 = all cores.
set -euo pipefail
shopt -s nullglob nocaseglob

if (($# < 6)); then
    echo "usage: $0 BUILD_DIR LIGHTS_DIR FLATS_DIR BIAS_DIR DARKS_DIR OUT_DIR [cfitsio|xisfconv]" >&2
    exit 2
fi
B=$1; L=${2%/}; FLAT=${3%/}; BIAS=${4%/}; DARK=${5%/}; O=${6%/}; BE=${7:-cfitsio}
for d in "$L" "$FLAT" "$BIAS" "$DARK"; do [[ -d $d ]] || { echo "run_end_to_end: not a directory: $d" >&2; exit 2; }; done
LIGHT_PREFIX=${LIGHT_PREFIX:-HD189733b_}
THREADS=${THREADS:-0}
HERE=$(cd "$(dirname "$0")" && pwd)

# fits_in DIR PATTERN  -> FITS files in DIR whose name matches *PATTERN*, sorted
fits_in() {
    local d=$1 p=$2 f out=()
    for f in "$d"/*"$p"*.fit "$d"/*"$p"*.fits "$d"/*"$p"*.fts; do [[ -f $f ]] && out+=("$f"); done
    ((${#out[@]})) && printf '%s\n' "${out[@]}" | sort
    return 0
}
need() {  # need WHAT FILES...
    local what=$1; shift
    if (($# == 0)); then
        echo "run_end_to_end: no $what frames found. FITS files present:" >&2
        for d in "${SEARCHED[@]}"; do echo "  $d:" >&2; fits_in "$d" "" | head -5 | sed 's/^/    /' >&2; done
        echo "Set the directory or name explicitly (see the header of $0)." >&2
        exit 1
    fi
    echo "$what: $# frames ($(basename "$1") ... $(basename "${@: -1}"))"
}

mapfile -t BIAS_F < <(fits_in "$BIAS" bias); ((${#BIAS_F[@]})) || mapfile -t BIAS_F < <(fits_in "$BIAS" "")
SEARCHED=("$BIAS"); need bias "${BIAS_F[@]}"
mapfile -t DARK_F < <(fits_in "$DARK" dark); ((${#DARK_F[@]})) || mapfile -t DARK_F < <(fits_in "$DARK" "")
SEARCHED=("$DARK"); need dark "${DARK_F[@]}"
SEARCHED=("$FLAT")
mapfile -t FLAT_F < <(fits_in "$FLAT" flat)
if ((${#FLAT_F[@]} == 0)) && [[ $(cd "$FLAT" && pwd) != $(cd "$L" && pwd) ]]; then
    mapfile -t FLAT_F < <(fits_in "$FLAT" "")
fi
need flat "${FLAT_F[@]}"
mapfile -t LIGHT_F < <(fits_in "$L" "")
LIGHT_N=0; for f in "${LIGHT_F[@]}"; do [[ ${f##*/} == "$LIGHT_PREFIX"* ]] && LIGHT_N=$((LIGHT_N + 1)); done
if ((LIGHT_N == 0)); then
    SEARCHED=("$L"); echo "run_end_to_end: no light frames named $LIGHT_PREFIX* in $L (set LIGHT_PREFIX)." >&2; need light
fi
echo "light: $LIGHT_N frames ($LIGHT_PREFIX*)"

mkdir -p "$O"
"$B/transitscope" master --backend "$BE" --type bias --out "$O/master_bias.fits" "${BIAS_F[@]}"
"$B/transitscope" master --backend "$BE" --type dark --out "$O/master_dark_10s.fits" "${DARK_F[@]}"
"$B/transitscope" master --backend "$BE" --type flat --bias "$O/master_bias.fits" --out "$O/master_flat.fits" "${FLAT_F[@]}"
"$B/transitscope" reduce --threads "$THREADS" --backend "$BE" --lights "$L" --prefix "$LIGHT_PREFIX" \
  --dark "$O/master_dark_10s.fits" --flat "$O/master_flat.fits" --bias "$O/master_bias.fits" --hot-adu 50 \
  --ra 20h00m43.713s --dec +22d42m39.07s --lat 50.93 --lon 5.40 --height 50 \
  --stars "$HERE/ref_stars.csv" --ref-index 99 --target 0 --out "$O/phot.csv"
"$B/transitscope" lightcurve --phot "$O/phot.csv" --T0 2454279.436714 --P 2.21857567 --T14 0.07527 \
  --target 0 --exclude-stars 10,32,36 --gain 0.26 --read-noise 4 --aperture-cm 20 --height 50 --out "$O/lc.csv"
"$B/transitscope" ld --teff 5050 --logg 4.59 --feh -0.03 --teff-err 50 --logg-err 0.02 --feh-err 0.08 --out "$O/ld.json"
"$B/transitscope" fit --threads "$THREADS" --lc "$O/lc.csv" --mode constrained --ld-json "$O/ld.json" --a-prior 8.863,0.05 --inc-prior 85.71,0.05 --out "$O/fit_constrained.json" --model-csv "$O/model_constrained.csv"
"$B/transitscope" fit --threads "$THREADS" --lc "$O/lc.csv" --mode free --ld-json "$O/ld.json" --a-prior 8.863,0.05 --inc-prior 85.71,0.05 --out "$O/fit_free.json" --model-csv "$O/model_free.csv"
