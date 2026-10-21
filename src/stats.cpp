// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/stats.h"

#include <algorithm>
#include <cmath>

namespace transit {

double median(std::vector<float> v) {
    if (v.empty()) return NAN;
    const size_t n = v.size();
    auto mid = v.begin() + n / 2;
    std::nth_element(v.begin(), mid, v.end());
    double m = *mid;
    if (n % 2 == 0) m = 0.5 * (m + *std::max_element(v.begin(), mid));
    return m;
}

ClippedStats sigma_clipped_stats(std::vector<float> v, double sigma, int maxiters) {
    ClippedStats s;
    for (int it = 0; it <= maxiters; ++it) {
        if (v.empty()) return s;
        s.median = median(v);
        double sum = 0, sum2 = 0;
        for (float x : v) sum += x;
        s.mean = sum / v.size();
        for (float x : v) sum2 += (x - s.mean) * (x - s.mean);
        s.stddev = std::sqrt(sum2 / v.size());
        s.n = v.size();
        if (it == maxiters) break;
        const double lo = s.median - sigma * s.stddev, hi = s.median + sigma * s.stddev;
        const size_t before = v.size();
        v.erase(std::remove_if(v.begin(), v.end(), [&](float x) { return x < lo || x > hi; }), v.end());
        if (v.size() == before) break;
    }
    return s;
}

}  // namespace transit
