// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <string>
#include <vector>

#include "transit/phot_table.h"
#include "transit/transit_model.h"

namespace transit {

struct Ephemeris { double T0 = 0, P = 1, T14 = 0.1; };  // BJD_TDB, days

struct EnsembleConfig {
    size_t target = 0;
    Ephemeris eph;
    double oot_margin = 0.003;            // days beyond T14/2 that still count as in transit
    double cooler_tol = 0.35;             // |CCD-TEMP - median| above this flags the frame
    double sat_adu = 45000;               // raw peak of the target above this flags the frame
    double min_rel_flux = 0.03;           // comparison candidates: at least this fraction of the target
    std::vector<size_t> exclude_stars;    // never use as comparison (e.g. on nebulosity)
    std::vector<double> radii;            // empty = try all radii in the table
    std::vector<size_t> comps;            // empty = automatic selection
    size_t min_comps = 3;
    // a star is dropped from the ensemble only if that lowers the light curve's residual RMS by more
    // than this fraction (a margin well above the RMS changes small photometry differences cause)
    double min_improvement = 0.01;
    // fixed geometry used to score comparison sets (free Tc, Rp/Rs, quadratic baseline)
    TransitParams score_model{0, 2.21857567, 0.155, 8.84, 85.7, 0.50, 0.19};
    // noise model
    double gain = 0.26, read_noise_e = 4.0, aperture_cm = 20.0, site_height_m = 0.0;
};

struct ApertureScore { double radius; size_t ncomps; double model_rms, oot_rms; };

struct EnsembleResult {
    double radius = 0;
    std::vector<size_t> comps;
    std::vector<double> weights;
    std::vector<double> flux, sigma;        // per frame (NaN where the target is not measured)
    std::vector<bool> good, oot;
    double tc_pred = 0;
    long epoch = 0;
    int n_flag_cooler = 0, n_flag_saturated = 0;
    std::vector<ApertureScore> scan;
    double sig_photon_median = 0, sig_scint_min = 0, sig_scint_max = 0;
};

EnsembleResult build_light_curve(const PhotTable& t, const EnsembleConfig& cfg);

/// Writes BJD_TDB, flux, sigma, airmass, good, oot (+ a comment header with the choices made).
void write_light_curve_csv(const std::string& path, const PhotTable& t, const EnsembleResult& r);

}  // namespace transit
