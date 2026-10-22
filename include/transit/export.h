// SPDX-License-Identifier: GPL-3.0-or-later
// Light curves for the amateur transit databases, in the plain-text layout their upload pages read
// (the layout HOPS writes for them):
//   ExoClock (www.exoclock.space)  JD_UTC at exposure start, relative flux, flux uncertainty
//   ETD (var.astro.cz/ETD)          JD_UTC (geocentric) at mid-exposure, relative flux, uncertainty
// Both get the comparison-ensemble light curve of the good frames, not detrended: the databases
// fit their own transit model and baseline. JD_UTC is recomputed from DATE-OBS of each frame.
#pragma once
#include <string>

#include "transit/project.h"

namespace transit {

enum class ExportTarget { ExoClock, ETD };
const char* export_target_name(ExportTarget t);   // "ExoClock", "ETD"
ExportTarget parse_export_target(const std::string& s);   // "exoclock" or "etd" (any case)

struct ExportResult {
    std::string path;
    size_t points = 0;
    double exposure_s = 0;
    std::string planet, filter;
    std::string instructions;   // what to choose in the upload form, and where to upload
};

/// Default file name in the project's output folder, e.g. "exoclock_HD189733b.txt".
std::string default_export_path(const Project& p, ExportTarget t);

/// Writes the file from the project's phot.csv and lc.csv (run the lightcurve stage first).
ExportResult export_light_curve(const Project& p, ExportTarget t, const std::string& path);

}  // namespace transit
