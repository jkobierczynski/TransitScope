// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <vector>

namespace transit {

struct ClippedStats { double mean = 0, median = 0, stddev = 0; size_t n = 0; };

double median(std::vector<float> v);  // by value: partially sorts a copy
/// Same semantics as astropy.stats.sigma_clipped_stats(sigma=3, maxiters=5):
/// median centre, population std, iterate until stable or maxiters.
ClippedStats sigma_clipped_stats(std::vector<float> v, double sigma = 3.0, int maxiters = 5);

}  // namespace transit
