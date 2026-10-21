// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "transit/ensemble.h"
#include "transit/phot_table.h"
#include "transit/progress.h"

namespace transit {

struct Gaussian { double mu, sigma; };

struct FitConfig {
    enum class Mode { Free, Constrained } mode = Mode::Constrained;
    double period = 2.21857567;
    Gaussian u1{0.50, 0.10}, u2{0.19, 0.10};     // quadratic limb darkening priors
    Gaussian a{8.863, 0.05}, inc{85.71, 0.05};   // used in Constrained mode (Agol et al. 2010)
    double rp0 = 0.155, a0 = 8.86, inc0 = 85.7;  // walker start
    size_t walkers = 40, steps = 16000, discard = 2000, thin = 10;
    uint64_t seed = 1;
    unsigned threads = 0;
};

/// Points to fit: the good frames of a light curve.
struct FitData {
    std::vector<size_t> frame;                    // frame index in the light curve
    std::vector<double> t, f, sigma;              // BJD_TDB, relative flux, per-point noise model
    double tc_pred = NAN;                         // predicted mid-transit time (from the ephemeris)
};

/// From the light-curve CSV of `transitscope lightcurve` (good frames, tc_pred from its header).
FitData load_fit_data(const std::string& lc_csv);
/// The same, in memory, from build_light_curve's result (no CSV rounding).
FitData fit_data_from(const PhotTable& phot, const EnsembleResult& lc);

/// Limb-darkening priors from the JSON of `transitscope ld` (u1, u2, u1_sigma, u2_sigma).
void read_ld_json(const std::string& path, FitConfig& cfg);

struct Quantiles { double med, lo, hi; };

struct FitResult {
    std::vector<std::string> names;               // Tc rp a inc u1 u2 c0 c1 c2 lnf
    std::vector<Quantiles> q;
    std::vector<double> best;                     // highest-posterior sample
    double beta = 1, tc_err_inflated_min = 0, rms_resid = 0, depth_pct = 0, b = 0, o_c_min = 0;
    double acceptance = 0, autocorr_max = 0, tc_pred = 0;
    size_t n_points = 0;
    std::vector<std::vector<double>> samples;     // posterior samples after discard and thinning: [sample][param]
    std::vector<double> sample_lnprob;            // log posterior of each sample
};

/// Per-point curves of the best-fit model, for plotting.
struct FitCurves {
    std::vector<double> t, flux, flux_detrended, model_transit, residual;
};

/// t, f, sigma: good frames only (BJD_TDB, relative flux, per-point noise model).
/// progress(step, steps) after every MCMC step; returning false throws Cancelled.
FitResult fit_transit(const std::vector<double>& t, const std::vector<double>& f, const std::vector<double>& sigma,
                      double tc_pred, const FitConfig& cfg, const ProgressFn& progress = {});
/// The same on FitData (uses its tc_pred).
FitResult fit_transit(const FitData& d, const FitConfig& cfg, const ProgressFn& progress = {});

FitCurves fit_curves(const std::vector<double>& t, const std::vector<double>& f, const FitResult& r, const FitConfig& cfg);

void write_fit_json(const std::string& path, const FitResult& r, const FitConfig& cfg);
/// time, flux, detrended flux, model, residual (for plotting)
void write_fit_model_csv(const std::string& path, const std::vector<double>& t, const std::vector<double>& f,
                         const FitResult& r, const FitConfig& cfg);

}  // namespace transit
