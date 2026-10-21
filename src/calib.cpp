// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/calib.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <stdexcept>

#include "transit/image_reader.h"
#include "transit/progress.h"

namespace transit {

void Calibration::build_hot_mask(const Image& bias, float threshold_adu) {
    if (!dark.same_shape(bias)) throw std::invalid_argument("hot mask: dark/bias shape mismatch");
    hot.assign(dark.px.size(), 0);
    hot_index_.clear();
    for (size_t i = 0; i < hot.size(); ++i) {
        hot[i] = (dark.px[i] - bias.px[i]) > threshold_adu;
        if (hot[i]) hot_index_.push_back(i);
    }
}

Image Calibration::apply(const Image& raw) const {
    if (!dark.empty() && !raw.same_shape(dark)) throw std::invalid_argument("calibrate: light/dark shape mismatch");
    if (!flat.empty() && !raw.same_shape(flat)) throw std::invalid_argument("calibrate: light/flat shape mismatch");
    Image out(raw.nx, raw.ny);
    const size_t n = out.px.size();
    float* o = out.px.data();
    const float* r = raw.px.data();
    if (!dark.empty() && !flat.empty()) {      // one pass: (raw - dark) / flat, same operation order
        const float *d = dark.px.data(), *f = flat.px.data();
        for (size_t i = 0; i < n; ++i) o[i] = (r[i] - d[i]) / (f[i] > 0.5f ? f[i] : 1.f);
    } else {
        for (size_t i = 0; i < n; ++i) o[i] = r[i];
        if (!dark.empty()) { const float* d = dark.px.data(); for (size_t i = 0; i < n; ++i) o[i] -= d[i]; }
        if (!flat.empty()) { const float* f = flat.px.data(); for (size_t i = 0; i < n; ++i) o[i] /= (f[i] > 0.5f ? f[i] : 1.f); }
    }
    if (!hot.empty()) {
        // 3x3 median of the calibrated image (edges: nearest pixel, = scipy median_filter mode 'reflect').
        // All medians are taken from the unrepaired image first, then written, as before.
        std::vector<size_t> local;                       // `hot` set by hand without build_hot_mask
        const std::vector<size_t>* idx = &hot_index_;
        if (hot_index_.empty()) {
            for (size_t i = 0; i < hot.size(); ++i) if (hot[i]) local.push_back(i);
            idx = &local;
        }
        const auto& hot_index = *idx;
        std::vector<float> rep(hot_index.size());
        std::array<float, 9> w{};
        for (size_t h = 0; h < hot_index.size(); ++h) {
            const int x = static_cast<int>(hot_index[h] % out.nx), y = static_cast<int>(hot_index[h] / out.nx);
            int k = 0;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                    w[k++] = out.at(std::clamp(x + dx, 0, out.nx - 1), std::clamp(y + dy, 0, out.ny - 1));
            std::nth_element(w.begin(), w.begin() + 4, w.end());
            rep[h] = w[4];
        }
        for (size_t h = 0; h < hot_index.size(); ++h) out.px[hot_index[h]] = rep[h];
    }
    return out;
}

Image make_master_flat(const std::vector<Image>& flats, const Image& bias) {
    std::vector<Image> sub;
    sub.reserve(flats.size());
    for (const auto& f : flats) {
        if (!f.same_shape(bias)) throw std::invalid_argument("master flat: shape mismatch");
        Image s = f;
        for (size_t i = 0; i < s.px.size(); ++i) s.px[i] -= bias.px[i];
        normalise_by_median(s);
        sub.push_back(std::move(s));
    }
    Image m = median_combine(sub);
    normalise_by_median(m);
    return m;
}

Image build_master(const std::string& backend, const std::string& type, const std::vector<std::string>& files,
                   const std::string& bias_path, const std::function<bool(size_t, size_t)>& progress) {
    if (type != "bias" && type != "dark" && type != "flat") throw std::invalid_argument("master type must be bias, dark or flat");
    if (files.empty()) throw std::runtime_error("no " + type + " frames");
    auto r = make_reader(backend);
    auto load = [&](const std::string& f) { r->open(f); return r->pixels(0); };
    std::vector<Image> stack;
    stack.reserve(files.size());
    for (size_t i = 0; i < files.size(); ++i) {
        stack.push_back(load(files[i]));
        if (progress && !progress(i + 1, files.size())) throw Cancelled();
    }
    if (type == "flat") {
        if (bias_path.empty()) throw std::invalid_argument("a master flat needs the master bias");
        return make_master_flat(stack, load(bias_path));
    }
    return median_combine(stack);
}

void write_master(const std::string& path, const Image& m, const std::string& type, size_t ncombine) {
    write_fits_float(path, m, {{"NCOMBINE", std::to_string(ncombine), "frames combined", HeaderCard::Type::Int},
                               {"IMAGETYP", "master " + type, "", HeaderCard::Type::String}});
}

std::vector<std::string> find_calibration_frames(const std::string& dir, const std::string& word,
                                                 const std::string& lights_dir) {
    namespace fs = std::filesystem;
    auto lower = [](std::string s) { for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); return s; };
    std::vector<std::string> named, all;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (!e.is_regular_file()) continue;
        const std::string ext = lower(e.path().extension().string());
        if (ext != ".fit" && ext != ".fits" && ext != ".fts") continue;
        all.push_back(e.path().string());
        if (lower(e.path().filename().string()).find(lower(word)) != std::string::npos) named.push_back(e.path().string());
    }
    std::error_code ec;
    const bool is_lights = !lights_dir.empty() && fs::equivalent(dir, lights_dir, ec);
    auto& out = !named.empty() || is_lights ? named : all;
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace transit
