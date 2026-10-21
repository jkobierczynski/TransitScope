// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "transit/image.h"

namespace transit {

/// Master frames for one camera setup. `dark` must match the light exposure
/// time and temperature (it then also carries the bias).
struct Calibration {
    Image dark;                  // bias + thermal, same exposure as lights
    Image flat;                  // bias-subtracted, normalised to 1
    std::vector<uint8_t> hot;    // 1 = hot pixel, replaced by 3x3 median after calibration

    /// hot = (dark - bias) > threshold_adu
    void build_hot_mask(const Image& bias, float threshold_adu);
    /// (raw - dark) / flat, flat values <= 0.5 treated as 1, hot pixels repaired.
    Image apply(const Image& raw) const;

private:
    std::vector<size_t> hot_index_;   // positions of `hot`, filled by build_hot_mask (apply() is const and thread-safe)
};

/// Master flat from raw flats: median of (flat - bias), normalised.
Image make_master_flat(const std::vector<Image>& flats, const Image& bias);

/// Master frame from files: "bias" and "dark" are median stacks, "flat" uses make_master_flat with
/// the master bias at `bias_path`. Written by write_master with NCOMBINE and IMAGETYP cards.
/// progress(done, total) per frame read; returning false throws Cancelled.
Image build_master(const std::string& backend, const std::string& type, const std::vector<std::string>& files,
                   const std::string& bias_path = {}, const std::function<bool(size_t, size_t)>& progress = {});
void write_master(const std::string& path, const Image& m, const std::string& type, size_t ncombine);

/// Calibration frames in `dir` (.fit/.fits/.fts, sorted): those whose name contains `word`
/// (case-insensitive); if none does and `dir` is not `lights_dir`, every FITS file in `dir`.
std::vector<std::string> find_calibration_frames(const std::string& dir, const std::string& word,
                                                 const std::string& lights_dir = {});

}  // namespace transit
