# SPDX-License-Identifier: GPL-3.0-or-later
"""Write small FITS files that exercise reader edge cases (run: python3 -I make_fits_edge.py fits_edge)."""
import sys, os, numpy as np, warnings
from astropy.io import fits
warnings.simplefilter("ignore")
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(0)
def hdr(h):
    h["OBJECT"] = ("O'Brien's field", "embedded quote")
    h["LONGSTR"] = "x" * 150                       # CONTINUE cards
    h["HIERARCH ESO DET CHIP TEMP"] = -120.25      # HIERARCH
    h["FLAG"] = True
    h["NEG"] = -42
    h["SMALL"] = 1.234567890123e-12
    h["CRVAL1"] = 300.18213541666667
    h["CD1_1"] = -4.0168e-4
    h["DATE-OBS"] = "2014-09-27T20:38:03.125"
    h["EMPTY"] = ""
    h.add_comment("a comment card"); h.add_history("a history card")
    return h
u16 = rng.integers(0, 65535, (64, 80)).astype(np.uint16)
h = fits.PrimaryHDU(u16); hdr(h.header); h.writeto(f"{out}/uint16_bzero.fits", overwrite=True)
i16 = rng.integers(-32768, 32767, (64, 80)).astype(np.int16)
h = fits.PrimaryHDU(i16); hdr(h.header); h.writeto(f"{out}/int16_negative.fits", overwrite=True)
f32 = rng.normal(100, 10, (64, 80)).astype(np.float32); f32[3, 5] = np.nan; f32[10, :] = np.nan
h = fits.PrimaryHDU(f32); hdr(h.header); h.writeto(f"{out}/float32_nan.fits", overwrite=True)
f64 = rng.normal(0, 1e6, (64, 80))
h = fits.PrimaryHDU(f64); hdr(h.header); h.writeto(f"{out}/float64.fits", overwrite=True)
i32 = rng.integers(-2**31, 2**31 - 1, (64, 80)).astype(np.int32)
h = fits.PrimaryHDU(i32); hdr(h.header); h.writeto(f"{out}/int32.fits", overwrite=True)
sc = fits.PrimaryHDU(rng.integers(0, 1000, (64, 80)).astype(np.int16)); sc.header["BSCALE"] = 0.5; sc.header["BZERO"] = 10.0
hdr(sc.header); sc.writeto(f"{out}/int16_bscale.fits", overwrite=True)
blank = rng.integers(0, 1000, (64, 80)).astype(np.int16); blank[0, 0] = -999
b = fits.PrimaryHDU(blank); b.header["BLANK"] = -999; hdr(b.header); b.writeto(f"{out}/int16_blank.fits", overwrite=True)
fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(u16, name="SCI"), fits.BinTableHDU.from_columns([fits.Column("a", "E", array=np.arange(3.))]),
              fits.ImageHDU(f32, name="ERR")]).writeto(f"{out}/multi_ext.fits", overwrite=True)
fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(u16.astype(np.int32), compression_type="RICE_1")]).writeto(f"{out}/rice.fits.fz", overwrite=True)
c3 = rng.integers(0, 65535, (3, 64, 80)).astype(np.uint16)
fits.PrimaryHDU(c3).writeto(f"{out}/rgb_cube.fits", overwrite=True)
print(sorted(os.listdir(out)))
