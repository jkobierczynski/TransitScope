// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace transit {

/// Single-plane float image, row-major, (0,0) = first pixel of the first row.
/// Pixel *centres* are at integer coordinates, matching the numpy reference.
struct Image {
    int nx = 0, ny = 0;
    std::vector<float> px;

    Image() = default;
    Image(int w, int h, float v = 0.f) : nx(w), ny(h), px(static_cast<size_t>(w) * h, v) {}
    float& at(int x, int y) { return px[static_cast<size_t>(y) * nx + x]; }
    float at(int x, int y) const { return px[static_cast<size_t>(y) * nx + x]; }
    bool empty() const { return px.empty(); }
    bool same_shape(const Image& o) const { return nx == o.nx && ny == o.ny; }
};

/// Pixel-wise median of a stack (all frames must share a shape).
Image median_combine(const std::vector<Image>& stack);

/// Divide by the median so the image is normalised to ~1.
void normalise_by_median(Image& img);

}  // namespace transit
