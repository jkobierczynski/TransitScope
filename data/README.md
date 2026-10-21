# Data files

| file | content | source | accuracy / status |
|---|---|---|---|
| `ld/claret2011_atlas_quad.csv` | quadratic limb-darkening coefficients, ATLAS atmospheres, v_turb = 2 km/s, methods L and F, bands B V R I J g' r' i' z', Teff 3500–10000 K, log g 0–5, [Fe/H] −5 … +1 | Claret A., Bloemen S., 2011, A&A 529, A75 (VizieR J/A+A/529/A75), converted from the EXOFASTv2 `quadld/` IDL save files (BSD-3-Clause, J. Eastman) with `tests/ld/convert_claret_sav.py` | as published; grid cells missing in the original are omitted |
| `passbands/sony_icx285al_qe.csv` | ICX285AL relative spectral response, 400–1000 nm | Sony ICX285AL datasheet, "Spectral Sensitivity Characteristics" (p. 9), https://www.1stvision.com/cameras/sensor_specs/ICX285.pdf | digitised from the datasheet graph: ±0.03 in response, ±3 nm |
| `passbands/astrodon_xop_bb.csv` | Astrodon Exo-Planet BB filter transmission | **model**: vendor states transmission starts near 500 nm and continues into the near-IR; no published table was found. Logistic edge at 500 nm, 0.95 plateau | `transitscope ld` adds a ±10 nm edge shift to the error budget |

Replace either passband file with a measured curve (two columns: `wavelength_nm,value`, `#` comments
allowed) and pass it with `--filter FILE` / `--qe FILE`.
