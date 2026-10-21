// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <vector>

#include "transit/image.h"

namespace transit {

struct ApertureConfig {
    std::vector<double> radii{8, 10, 12, 14, 16, 18, 20};
    double ann_in = 26, ann_out = 40;    // background annulus, px
    double centroid_radius = 16;         // flux-weighted centroid window, px
    int centroid_iter = 4;
    int half_box = 45;                   // cut-out half size around a star
    double edge_margin = 60;             // stars closer than this to an edge are skipped
};

struct StarMeasure {
    bool valid = false;
    double x = 0, y = 0;
    std::vector<double> flux;            // background-subtracted sum per radius (hard-edged pixels)
    std::vector<int> npix;               // pixels in each aperture
    double bg = 0, bg_sd = 0;            // sigma-clipped median / std in the annulus
    double peak_raw = 0;                 // max raw ADU within r < 12 (linearity check)
    double r90 = 0;                      // 90% encircled-energy radius, px (donut-safe size metric)
};

/// Flux-weighted centroid with a median annulus background, iterated.
void centroid(const Image& img, double& x, double& y, const ApertureConfig& cfg);

/// Measure one star at (x, y). `raw` is the uncalibrated frame for the peak check.
StarMeasure measure(const Image& img, const Image& raw, double x, double y, const ApertureConfig& cfg);

/// Simple peak finder for the reference frame (Gaussian smooth + local maxima).
struct Detection { double x, y, peak; };
std::vector<Detection> find_stars(const Image& img, double sigma_px = 4, int max_box = 31,
                                  double nsigma = 25, size_t nmax = 40);

}  // namespace transit
