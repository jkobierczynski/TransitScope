// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/export.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "transit/fit.h"
#include "transit/phot_table.h"
#include "transit/timing.h"

namespace transit {

const char* export_target_name(ExportTarget t) { return t == ExportTarget::ExoClock ? "ExoClock" : "ETD"; }

ExportTarget parse_export_target(const std::string& s) {
    std::string l = s;
    for (auto& c : l) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (l == "exoclock") return ExportTarget::ExoClock;
    if (l == "etd") return ExportTarget::ETD;
    throw std::invalid_argument("export target must be exoclock or etd, not '" + s + "'");
}

namespace {
std::string planet_name(const Project& p) {
    if (!p.planet.empty()) return p.planet;
    std::string t = p.target_name;
    t.erase(std::remove(t.begin(), t.end(), ' '), t.end());
    return t.empty() ? std::string("planet") : t + "b";
}
std::string safe_file_part(std::string s) {
    for (auto& c : s)
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '-' && c != '_' && c != '.') c = '_';
    return s;
}
}  // namespace

std::string default_export_path(const Project& p, ExportTarget t) {
    return p.output(std::string(t == ExportTarget::ExoClock ? "exoclock_" : "etd_") + safe_file_part(planet_name(p)) + ".txt");
}

ExportResult export_light_curve(const Project& p, ExportTarget t, const std::string& path) {
    const PhotTable phot = load_phot_csv(p.output("phot.csv"));
    const FitData lc = load_fit_data(p.output("lc.csv"));   // good frames: frame index, flux, sigma
    if (lc.t.empty()) throw std::runtime_error("no good frames in " + p.output("lc.csv"));

    ExportResult r;
    r.path = path;
    r.planet = planet_name(p);
    r.filter = p.filter_name;
    std::vector<double> exp;
    for (size_t f : lc.frame)
        if (f < phot.nf && std::isfinite(phot.exptime[f])) exp.push_back(phot.exptime[f]);
    if (!exp.empty()) {
        std::nth_element(exp.begin(), exp.begin() + exp.size() / 2, exp.end());
        r.exposure_s = exp[exp.size() / 2];
    }
    char expbuf[32];
    std::snprintf(expbuf, sizeof expbuf, "%g", r.exposure_s);

    const bool exoclock = t == ExportTarget::ExoClock;
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o << "#Planet: " << r.planet << "\n"
      << (exoclock ? "#Time format: JD_UTC\n#Time stamp: Exposure start\n#Flux format: Flux\n"
                   : "#Time format: JD_UTC (geocentric)\n#Time stamp: Mid-exposure\n#Flux format: Flux (in flux)\n")
      << "#Filter: " << r.filter << "\n"
      << "#Exposure time in seconds: " << expbuf << "\n"
      << "#Time                   Flux            Flux uncertainty\n";
    for (size_t i = 0; i < lc.frame.size(); ++i) {
        const size_t f = lc.frame[i];
        if (f >= phot.nf) throw std::runtime_error("lc.csv and phot.csv do not belong together (frame " + std::to_string(f) + ")");
        const FrameTime ft = frame_time(phot.date_obs[f], phot.exptime[f], p.reduce.target, p.reduce.site);
        const double jd = exoclock ? ft.jd_utc_start : ft.jd_utc_start + phot.exptime[f] / 2 / 86400.0;
        char line[128];
        std::snprintf(line, sizeof line, "%.10f\t%.10f\t%.10f\n", jd, lc.f[i], lc.sigma[i]);
        o << line;
        ++r.points;
    }
    if (!o) throw std::runtime_error("error writing " + path);

    std::ostringstream ins;
    if (exoclock) {
        ins << "Upload at https://www.exoclock.space/upload/ (log in; your telescope must be registered).\n"
            << "Choose in the form:\n"
            << "  Planet: " << r.planet << "\n  Time format: JD_UTC\n  Time stamp: Exposure start\n"
            << "  Flux format: Flux\n  Filter: " << (r.filter.empty() ? "(the filter you used)" : r.filter) << "\n"
            << "  Exposure time in seconds: " << expbuf << "\n";
    } else {
        ins << "Upload at https://var.astro.cz/ETD (log in; Var-Astro light-curve upload).\n"
            << "Choose in the form:\n"
            << "  Planet: " << r.planet << "\n  Time: JD, geocentric (UTC), mid-exposure\n"
            << "  Brightness column: in flux\n  Filter: " << (r.filter.empty() ? "(the filter you used)" : r.filter) << "\n";
    }
    ins << "The file holds " << r.points << " good frames, comparison-ensemble flux, not detrended.";
    r.instructions = ins.str();
    return r;
}

}  // namespace transit
