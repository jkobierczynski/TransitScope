// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <string>
#include <vector>

namespace transit {

/// The `transitscope reduce` CSV, loaded as dense frame x star (x radius) arrays.
/// Invalid measurements are NaN.
struct PhotTable {
    size_t nf = 0, ns = 0, nr = 0;
    std::vector<double> radii;
    std::vector<std::string> file, date_obs;
    std::vector<double> exptime, temp, bjd, airmass;          // per frame
    std::vector<double> x, y, bg, bg_sd, peak, r90;           // per frame x star
    std::vector<double> flux, npix;                           // per frame x star x radius

    size_t fs(size_t f, size_t s) const { return f * ns + s; }
    size_t fsr(size_t f, size_t s, size_t r) const { return (f * ns + s) * nr + r; }
    double F(size_t f, size_t s, size_t r) const { return flux[fsr(f, s, r)]; }
    int radius_index(double R) const;
};

PhotTable load_phot_csv(const std::string& path);

}  // namespace transit
