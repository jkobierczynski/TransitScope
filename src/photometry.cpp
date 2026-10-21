// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/photometry.h"

#include <algorithm>
#include <cmath>

#include "transit/stats.h"

namespace transit {
namespace {

struct Box { int x0, x1, y0, y1; };  // [x0, x1) x [y0, y1)

Box box_around(const Image& img, double x, double y, int half) {
    // numpy reference uses int(round(x)): round-half-to-even == nearbyint in the default FP mode
    const int xi = static_cast<int>(std::nearbyint(x)), yi = static_cast<int>(std::nearbyint(y));
    return {std::max(xi - half, 0), std::min(xi + half + 1, img.nx), std::max(yi - half, 0),
            std::min(yi + half + 1, img.ny)};
}

// Box `b` clipped to the pixels that can lie within `r` of (x, y). Iterating the sub-box in
// row-major order visits exactly the same pixels, in the same order, as the full box with a
// distance test, so sums are bit-identical to the unrestricted loops.
Box shrink(const Box& b, double x, double y, double r) {
    return {std::max(b.x0, static_cast<int>(std::floor(x - r))), std::min(b.x1, static_cast<int>(std::ceil(x + r)) + 1),
            std::max(b.y0, static_cast<int>(std::floor(y - r))), std::min(b.y1, static_cast<int>(std::ceil(y + r)) + 1)};
}

std::vector<float> annulus(const Image& img, const Box& b0, double x, double y, double rin, double rout) {
    const Box b = shrink(b0, x, y, rout);
    const double r2in = rin * rin, r2out = rout * rout;
    std::vector<float> v;
    v.reserve(static_cast<size_t>(3.2 * (r2out - r2in)) + 16);
    for (int yy = b.y0; yy < b.y1; ++yy) {
        const double dy2 = (yy - y) * (yy - y);
        const float* row = &img.px[static_cast<size_t>(yy) * img.nx];
        for (int xx = b.x0; xx < b.x1; ++xx) {
            const double d2 = (xx - x) * (xx - x) + dy2;
            if (d2 > r2in && d2 < r2out) v.push_back(row[xx]);
        }
    }
    return v;
}

// reflect-mode index (scipy 'reflect': d c b a | a b c d | d c b a)
inline int reflect(int i, int n) {
    while (i < 0 || i >= n) i = (i < 0) ? -i - 1 : 2 * n - i - 1;
    return i;
}

Image gaussian_blur(const Image& in, double sigma) {
    const int r = static_cast<int>(4.0 * sigma + 0.5);
    std::vector<double> k(2 * r + 1);
    double s = 0;
    for (int i = -r; i <= r; ++i) s += k[i + r] = std::exp(-0.5 * i * i / (sigma * sigma));
    for (auto& v : k) v /= s;
    Image tmp(in.nx, in.ny), out(in.nx, in.ny);
    for (int y = 0; y < in.ny; ++y)
        for (int x = 0; x < in.nx; ++x) {
            double a = 0;
            for (int i = -r; i <= r; ++i) a += k[i + r] * in.at(reflect(x + i, in.nx), y);
            tmp.at(x, y) = static_cast<float>(a);
        }
    for (int y = 0; y < in.ny; ++y)
        for (int x = 0; x < in.nx; ++x) {
            double a = 0;
            for (int i = -r; i <= r; ++i) a += k[i + r] * tmp.at(x, reflect(y + i, in.ny));
            out.at(x, y) = static_cast<float>(a);
        }
    return out;
}

Image max_filter(const Image& in, int size) {
    const int h = size / 2;
    Image tmp(in.nx, in.ny), out(in.nx, in.ny);
    for (int y = 0; y < in.ny; ++y)
        for (int x = 0; x < in.nx; ++x) {
            float m = -INFINITY;
            for (int i = -h; i <= h; ++i) m = std::max(m, in.at(reflect(x + i, in.nx), y));
            tmp.at(x, y) = m;
        }
    for (int y = 0; y < in.ny; ++y)
        for (int x = 0; x < in.nx; ++x) {
            float m = -INFINITY;
            for (int i = -h; i <= h; ++i) m = std::max(m, tmp.at(x, reflect(y + i, in.ny)));
            out.at(x, y) = m;
        }
    return out;
}

}  // namespace

void centroid(const Image& img, double& x, double& y, const ApertureConfig& cfg) {
    const double r2c = cfg.centroid_radius * cfg.centroid_radius;
    for (int it = 0; it < cfg.centroid_iter; ++it) {
        const Box b = box_around(img, x, y, cfg.half_box);
        const double bg = median(annulus(img, b, x, y, cfg.ann_in, cfg.ann_out));
        const Box c = shrink(b, x, y, cfg.centroid_radius);
        double sw = 0, sx = 0, sy = 0;
        for (int yy = c.y0; yy < c.y1; ++yy) {
            const double dy2 = (yy - y) * (yy - y);
            const float* row = &img.px[static_cast<size_t>(yy) * img.nx];
            for (int xx = c.x0; xx < c.x1; ++xx) {
                if ((xx - x) * (xx - x) + dy2 >= r2c) continue;
                const double w = std::max(0.0, row[xx] - bg);
                sw += w; sx += w * xx; sy += w * yy;
            }
        }
        if (sw <= 0) break;
        x = sx / sw;
        y = sy / sw;
    }
}

StarMeasure measure(const Image& img, const Image& raw, double x, double y, const ApertureConfig& cfg) {
    StarMeasure m;
    m.x = x; m.y = y;
    const Box b = box_around(img, x, y, cfg.half_box);
    const auto st = sigma_clipped_stats(annulus(img, b, x, y, cfg.ann_in, cfg.ann_out));
    m.bg = st.median; m.bg_sd = st.stddev;

    const size_t nr = cfg.radii.size();
    constexpr int kEE = 22;                       // encircled-energy profile out to r = 22 px
    std::vector<double> r2(nr);
    double rmax = kEE;
    for (size_t i = 0; i < nr; ++i) { r2[i] = cfg.radii[i] * cfg.radii[i]; rmax = std::max(rmax, cfg.radii[i]); }
    m.flux.assign(nr, 0.0); m.npix.assign(nr, 0);
    std::vector<double> ring(kEE + 1, 0.0);       // ring[k]: pixels with k-1 <= r < k
    double tot = 0, peak = -INFINITY;
    const Box c = shrink(b, x, y, rmax);
    for (int yy = c.y0; yy < c.y1; ++yy) {
        const double dy2 = (yy - y) * (yy - y);
        const float* row = &img.px[static_cast<size_t>(yy) * img.nx];
        const float* rrow = &raw.px[static_cast<size_t>(yy) * raw.nx];
        for (int xx = c.x0; xx < c.x1; ++xx) {
            const double d2 = (xx - x) * (xx - x) + dy2;
            const double v = row[xx] - m.bg;
            for (size_t i = 0; i < nr; ++i)
                if (d2 < r2[i]) { m.flux[i] += v; ++m.npix[i]; }
            if (d2 < kEE * kEE) {
                tot += v;
                const int k = static_cast<int>(std::sqrt(d2)) + 1;   // smallest integer k with r < k
                ring[std::min(k, kEE)] += v;
                if (d2 < 144.0) peak = std::max(peak, static_cast<double>(rrow[xx]));
            }
        }
    }
    m.peak_raw = peak;
    m.r90 = kEE - 1;
    if (tot > 0) {
        double cum = 0;
        for (int k = 1; k < kEE; ++k) {
            cum += ring[k];
            if (cum / tot >= 0.9) { m.r90 = k; break; }
        }
    }
    m.valid = true;
    return m;
}

std::vector<Detection> find_stars(const Image& img, double sigma_px, int max_box, double nsigma, size_t nmax) {
    const Image sm = gaussian_blur(img, sigma_px);
    const Image mx = max_filter(sm, max_box);
    const auto st = sigma_clipped_stats(sm.px);
    const double thr = st.median + nsigma * st.stddev;
    std::vector<Detection> det;
    for (int y = 0; y < sm.ny; ++y)
        for (int x = 0; x < sm.nx; ++x)
            if (sm.at(x, y) == mx.at(x, y) && sm.at(x, y) > thr) det.push_back({double(x), double(y), sm.at(x, y)});
    std::sort(det.begin(), det.end(), [](auto& a, auto& b) { return a.peak > b.peak; });
    if (det.size() > nmax) det.resize(nmax);
    return det;
}

}  // namespace transit
