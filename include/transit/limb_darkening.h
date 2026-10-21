// SPDX-License-Identifier: GPL-3.0-or-later
// Quadratic limb-darkening coefficients for an arbitrary passband (filter x detector QE),
// from the Claret & Bloemen (2011) ATLAS tables for standard bands.
//
//   1. trilinear interpolation of the tables in [Fe/H], Teff, log g, per standard band
//   2. coefficients vs wavelength: piecewise linear through the bands' effective wavelengths
//   3. photon-weighted average over the passband: w(l) = T(l) QE(l) B_l(Teff) l
//      (exact for the disc-integrated I(mu) if the coefficients are weighted by I_l(1), which
//      the blackbody approximates)
//   4. error budget: L vs F fitting method, wavelength interpolation (leave-one-band-out),
//      passband uncertainty (filter edge shift, QE tilt), stellar parameters, model floor.
#pragma once
#include <map>
#include <string>
#include <vector>

namespace transit {

/// Two-column curve (wavelength in nm, response). '#' lines and a non-numeric header are skipped.
struct Curve {
    std::vector<double> nm, y;
    double at(double lambda_nm) const;   // linear interpolation, 0 outside the tabulated range
};
Curve load_curve_csv(const std::string& path);

/// A data file: `name` itself if it exists, else <data dir>/<sub>/<name>[.csv], where the data dirs
/// are $TRANSIT_DATA_DIR, the source tree's data/ (build time) and ./data. Throws if not found.
std::string find_data_file(const std::string& name, const std::string& sub);

struct StarParams { double teff, logg, feh; };

class ClaretTable {
public:
    void load(const std::string& path);   // data/ld/claret2011_atlas_quad.csv
    /// (a, b) for one band and method ('L' or 'F'); throws if a needed grid cell is missing
    /// or the star is outside the grid.
    std::pair<double, double> at(char method, const std::string& band, const StarParams& s) const;
    std::vector<std::string> bands() const;
private:
    std::vector<double> feh_, teff_, logg_;
    // key: method + band; value: a/b on the full grid, NaN where the table has no entry
    std::map<std::string, std::vector<double>> a_, b_;
    size_t idx(size_t i, size_t j, size_t k) const { return (i * teff_.size() + j) * logg_.size() + k; }
};

/// Effective wavelength (nm) of the standard bands used for the wavelength interpolation.
double band_wavelength_nm(const std::string& band);

struct LdConfig {
    std::vector<std::string> bands = {"B", "V", "R", "I", "J", "Sloang", "Sloanr", "Sloani", "Sloanz"};
    double filter_shift_nm = 10;   // +- uncertainty of the filter cut-on
    double qe_tilt = 0.10;         // +- fractional QE tilt across 400-1000 nm
    StarParams star_err{0, 0, 0};  // 1-sigma stellar parameter errors
    double model_floor = 0.05;     // ATLAS vs observed limb darkening, added in quadrature
    double grid_step_nm = 1;
};

struct LdBudget { double method = 0, interp = 0, passband = 0, stellar = 0, floor = 0, total = 0; };

struct LdResult {
    double u1 = 0, u2 = 0;               // mean of the L and F results
    double u1_L = 0, u2_L = 0, u1_F = 0, u2_F = 0;
    double lambda_eff_nm = 0;            // photon-weighted mean wavelength of the passband
    double lambda_lo_nm = 0, lambda_hi_nm = 0;   // 5 % / 95 % points of the photon distribution
    LdBudget s1, s2;
    struct Band { std::string name; double lambda_nm, aL, bL, aF, bF, loo_a, loo_b; bool in_passband; };
    std::vector<Band> per_band;
};

LdResult passband_limb_darkening(const ClaretTable& tab, const StarParams& star, const Curve& filter,
                                 const Curve& qe, const LdConfig& cfg = {});

void write_ld_json(const std::string& path, const LdResult& r, const StarParams& star, const LdConfig& cfg);

}  // namespace transit
