// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/image.h"

#include <algorithm>
#include <stdexcept>

#include "transit/stats.h"

namespace transit {

Image median_combine(const std::vector<Image>& stack) {
    if (stack.empty()) throw std::invalid_argument("median_combine: empty stack");
    Image out(stack[0].nx, stack[0].ny);
    for (const auto& im : stack)
        if (!im.same_shape(out)) throw std::invalid_argument("median_combine: shape mismatch");
    const size_t n = stack.size();
    std::vector<float> col(n);
    for (size_t i = 0; i < out.px.size(); ++i) {
        for (size_t k = 0; k < n; ++k) col[k] = stack[k].px[i];
        auto mid = col.begin() + n / 2;
        std::nth_element(col.begin(), mid, col.end());
        float m = *mid;
        if (n % 2 == 0) m = 0.5f * (m + *std::max_element(col.begin(), mid));
        out.px[i] = m;
    }
    return out;
}

void normalise_by_median(Image& img) {
    const float m = static_cast<float>(median(img.px));
    if (m == 0) throw std::runtime_error("normalise_by_median: zero median");
    for (auto& v : img.px) v /= m;
}

}  // namespace transit
