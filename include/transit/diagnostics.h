// SPDX-License-Identifier: GPL-3.0-or-later
// Red-noise and systematics diagnostics of a fitted light curve's residuals:
//   - RMS of the residuals binned over n consecutive points against the white-noise expectation
//     (time-averaging method, Pont et al. 2006, Winn et al. 2008): beta(n) = rms(n) / white(n)
//   - autocorrelation of the residuals at lags 1..L
//   - correlation of the residuals with per-frame quantities that commonly drive systematics:
//     airmass, sky background, image size (r90), drift in x and y, CCD temperature, and the summed
//     comparison-star flux (transparency). Significance uses the effective number of points
//     N / beta^2, so red noise does not make chance correlations look significant.
#pragma once
#include <string>
#include <vector>

#include "transit/fit.h"
#include "transit/phot_table.h"

namespace transit {

struct BinRms {
    size_t n = 1;             // points per bin
    size_t bins = 0;          // number of bins
    double rms = 0;           // RMS of the bin means
    double white = 0;         // expected for white noise: rms(1)/sqrt(n) * sqrt(bins/(bins-1))
    double beta = 1;          // rms / white
};

struct Correlation {
    std::string name, unit;
    size_t n = 0;
    double r = 0;             // Pearson correlation of residual and quantity (NaN if the quantity is constant)
    double threshold = 0;     // |r| above this is significant (3 sigma, effective N)
    double slope_ppt = 0;     // residual change per unit of the quantity, ppt
    double span_ppt = 0;      // residual change across the 5-95 % range of the quantity, ppt
    bool flagged = false;
};

struct Diagnostics {
    size_t n = 0;
    double rms_ppt = 0;                   // unbinned residual RMS
    std::vector<BinRms> bins;
    double beta = 1;                      // median beta(n) for n = 5, 10, 15, 20, 30, at least 1 (as in the fit)
    double beta_max = 1;
    double n_eff = 0;                     // n / beta^2
    std::vector<double> acf;              // lags 1..acf.size()
    double acf_threshold = 0;             // 3 / sqrt(n)
    std::vector<Correlation> correlations;
    std::vector<std::string> warnings;    // plain-language findings
};

/// A per-point quantity to correlate with (same length and order as the residuals).
struct Quantity { std::string name, unit; std::vector<double> values; };

Diagnostics compute_diagnostics(const std::vector<double>& residual, const std::vector<Quantity>& quantities);

/// The per-frame quantities of the fitted points: frames = FitData::frame, star = target, the
/// comparison stars and aperture radius of the light curve (lc.csv header or EnsembleResult).
std::vector<Quantity> frame_quantities(const PhotTable& phot, const std::vector<size_t>& frames, size_t target,
                                       const std::vector<size_t>& comps, double radius);

/// Comparison stars and radius from the header of lc.csv.
void read_lc_header(const std::string& lc_csv, std::vector<size_t>& comps, double& radius);

/// Everything for one fitted light curve: residuals of the best fit, and the frame quantities from
/// phot.csv (phot may be null: then only the time-correlation part is computed).
Diagnostics diagnose_fit(const PhotTable* phot, const std::string& lc_csv, const FitData& data,
                         const std::vector<double>& residual, size_t target);

/// The residual column of model_<mode>.csv (write_fit_model_csv).
std::vector<double> load_residuals(const std::string& model_csv);

void write_diagnostics_json(const std::string& path, const Diagnostics& d);
std::string diagnostics_summary(const Diagnostics& d);   // a few lines for the log / terminal

}  // namespace transit
