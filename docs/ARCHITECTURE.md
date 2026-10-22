# Architecture

Three layers: **libtransit** (all science, no UI) → **transitscope** CLI (scriptable, headless, regression-testable) → **transitscope-gui** (Qt 6 Widgets, `gui/`; edits the same project file the CLI reads and runs `run_project` in a worker thread).

## Module map — Python reference → C++

| Python reference | libtransit module | Status |
|---|---|---|
| `astropy.io.fits` | `image_reader.h` + `io/cfitsio_reader.cpp`, `io/xisfconv_reader.cpp` | both done; bit-identical on all test files except BLANK handling |
| `np.median` stacks, master frames | `image.h` (`median_combine`), `calib.h` (`make_master_flat`) | done |
| `(raw - dark) / flat`, hot-pixel 3×3 median | `calib.h` (`Calibration::apply`) | done |
| `sigma_clipped_stats` | `stats.h` | done (astropy semantics) |
| `find_stars` (Gaussian + max filter) | `photometry.h` (`find_stars`) | done |
| FFT coarse shift | — replaced by sequential tracking (`reduce.h`) | done |
| reduce loop (read, calibrate, track, measure, CSV) | `reduce.h`: `reduce_sequence`, `write_phot_csv`, `to_phot_table`, `find_light_frames`, `load_calibration` | done in 0.5.0; CLI output unchanged |
| `centroid`, `measure` (apertures, annulus, r90, peak) | `photometry.h` | done |
| `astropy.time` → BJD_TDB, AltAz → airmass | `timing.h` (ERFA: `eraUtctai`, `eraTaitt`, `eraDtdb`, `eraEpv00`, `eraAtco13`) | done |
| reduce CSV | `phot_table.h` | done |
| `lc.py`: OOT mask, weighted ensemble, model-residual comp selection, noise model | `ensemble.h` (+ `lsq.h`: polyfit, Levenberg–Marquardt) | done; identical to < 1e-9 |
| batman quadratic model | `transit_model.h` (Mandel–Agol, AGM + Bulirsch elliptic integrals) | done; 2.4e-8 vs batman |
| `fit.py`: emcee, priors, β red-noise factor | `mcmc.h` (stretch move, red/blue split, emcee autocorrelation estimator), `fit.h` (`fit_transit` with progress/cancel, `FitData`, posterior samples, `fit_curves`) | done; within 0.06 σ of emcee |
| EXOFAST `quadld` + passband weighting | `limb_darkening.h` (Claret 2011 tables, filter × QE × blackbody average, error budget) | done; numpy check < 1e-6 |
| `plot.py` | `gui/PlotWidget` (QPainter), `gui/ResultsView` | done in 0.8.0 |

## Design decisions

- **Backend-neutral I/O.** The pipeline only sees `ImageReader`; `--backend cfitsio|xisfconv` picks one per run, and `fitscompare` diffs both on the same file. Both backends number images the same way (n-th HDU holding a ≥ 2-D image) and deliver physical values in FITS row order. The project is GPL-3.0-or-later; libxisfconv (LGPL-3.0-or-later) links cleanly.
- **Sequential tracking, parallel measuring.** Start at the reference frame and walk outward, re-centroiding the target and carrying the offset (follows the ~150 px drift in the 2014 data without an FFT dependency). Only this step is serial; reading/calibrating and measuring run on a FIFO thread pool (`thread_pool.h`) with per-thread readers. Calibration is `const` and thread-safe.
- **The CLI is a thin front end.** `reduce_sequence(config, calibration, progress)` does the whole reduction and returns the results in memory (`ReduceResult`, or a `PhotTable` via `to_phot_table`); the CLI only parses flags and writes the CSV. The progress callback runs on the calling thread once per frame and can cancel the run (it then waits for queued work and throws `Cancelled`), which is what a GUI worker thread needs. `fit_transit` works the same way, per MCMC step (`progress.h`), takes its points in memory (`fit_data_from`) or from the CSV (`load_fit_data`), and returns the posterior samples and the best-fit curves for plotting.
- **Bit-exact optimisation.** Distance tests use squared distances on sub-boxes that contain exactly the pixels the full-box loops accepted, visited in the same order, so all sums are bit-identical to the reference implementation; the encircled-energy profile is binned once per pixel.
- **Timing.** DATE-OBS is taken as UTC exposure start; mid-exposure → TAI → TT → TDB (+ periodic TDB−TT) → Rømer delay from `eraEpv00` (BCRS Earth position). Omitted: observer geocentric offset (< 21 ms) and Shapiro delay (< 0.1 ms). Agrees with astropy to 0.02 s.
- **Fit.** Quadratic-LD transit × quadratic baseline in time, noise scale exp(lnf) on the per-frame noise model, Gaussian LD priors; `constrained` adds a/Rs and i priors. Walkers are evaluated in parallel threads, with all random numbers drawn serially, so a seed gives the same chain on any thread count.
- **Project file.** `Project` holds the configs of every stage (`ReduceConfig`, `EnsembleConfig`, `LdSettings`, `FitConfig`) plus frame sources and paths. Values shared by stages (period, site height, target star) are written once and copied into each config on load. The reader rejects unknown keys and wrong types and reports every problem with its line number; `project_to_toml` writes a file that loads back to the same project (save is atomic). Comments are not preserved on save; the GUI writes the canonical layout.
- **GUI.** `gui/SettingsEditor` takes its list of keys, their types and tooltips from `project_template()`, so the file format has one definition; the tree is turned back into TOML and parsed with `parse_project` on every edit, and `check_project_inputs` fills the Problems list. `gui/ImageView` shows the reference frame in FITS orientation and edits the star list; `gui/MainWindow` keeps the star indices in `target_star`, `exclude_stars` and `comparison_stars` consistent when stars are removed. Runs use `run_project` on a copy of the project in a `QThread`; progress and log lines are posted to the GUI thread (throttled), cancelling sets an atomic flag the progress callback returns.
- **Long-format CSV between stages.** `reduce` writes one row per frame × star × all radii, so aperture and comparison choices are made later without re-measuring images (the GUI's instant comparison-star toggling relies on this cache).

## Known limitations / roadmap

1. Done in 0.2.0: ensemble light curve and MCMC fit (see README for the regression results).
2. Comparison selection is greedy and keeps any star that lowers the residual RMS by > 0.2 %; a faint, noisy star can tip in or out on tiny photometry differences. Options: a minimum weight, or selecting on a held-out half of the frames.
3. **Fractional-pixel apertures** (exact circle/pixel overlap). Hard-edged apertures make flux jump by ~3e-4 when one boundary pixel moves in/out — visible in the regression test at r = 8 px.
4. Defocus mode: donut-template centroid and encircled-energy sizing (`r90` is already logged).
5. Done in 0.4.0: limb-darkening priors from filter × camera QE (`transitscope ld`). Open: a measured XOP-BB transmission curve, PHOENIX tables for M dwarfs, and power-2 / non-linear laws.
6. Done in 0.7.0: TOML project file (`project.h`: `load_project`, `save_project`, `run_project`; `transitscope run`, `project new|check`). Open: plate solving via local `solve-field`; Gaia comparison-star candidates; Exoplanet Watch (AAVSO) export. ETD and ExoClock export done in 0.9.0 (`export.h`, `transitscope export`, Results tab).
7. Done in 0.3.0: pipelined multithreaded `reduce` and faster photometry kernels (see README, Performance).
8. The MCMC fit now dominates the run time (~40 s per mode on 2 cores). 16 000 steps is > 100 autocorrelation times; ~6 000 would do, or an automatic stop at 50 τ.
