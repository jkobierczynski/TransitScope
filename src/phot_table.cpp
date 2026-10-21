// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/phot_table.h"

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace transit {

int PhotTable::radius_index(double R) const {
    for (size_t i = 0; i < radii.size(); ++i)
        if (std::fabs(radii[i] - R) < 1e-9) return static_cast<int>(i);
    return -1;
}

static std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> v;
    std::string cur;
    for (char c : s) {
        if (c == ',') { v.push_back(cur); cur.clear(); }
        else if (c != '\r') cur += c;
    }
    v.push_back(cur);
    return v;
}

static double num(const std::string& s) {
    if (s.empty() || s == "nan" || s == "-nan") return NAN;
    return std::stod(s);
}

PhotTable load_phot_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(in, line);
    const auto head = split(line);
    std::map<std::string, size_t> col;
    std::vector<size_t> fcol, ncol;
    PhotTable t;
    for (size_t i = 0; i < head.size(); ++i) {
        col[head[i]] = i;
        if (head[i].rfind("flux_r", 0) == 0) { fcol.push_back(i); t.radii.push_back(std::stod(head[i].substr(6))); }
        if (head[i].rfind("npix_r", 0) == 0) ncol.push_back(i);
    }
    for (const char* k : {"frame", "file", "star", "valid", "bjd_tdb", "airmass", "x", "y", "bg"})
        if (!col.count(k)) throw std::runtime_error(path + ": missing column " + k);
    t.nr = t.radii.size();

    std::vector<std::vector<std::string>> rows;
    size_t nf = 0, ns = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        rows.push_back(split(line));
        nf = std::max(nf, static_cast<size_t>(std::stoul(rows.back()[col["frame"]])) + 1);
        ns = std::max(ns, static_cast<size_t>(std::stoul(rows.back()[col["star"]])) + 1);
    }
    t.nf = nf; t.ns = ns;
    t.file.assign(nf, ""); t.date_obs.assign(nf, "");
    for (auto* v : {&t.exptime, &t.temp, &t.bjd, &t.airmass}) v->assign(nf, NAN);
    for (auto* v : {&t.x, &t.y, &t.bg, &t.bg_sd, &t.peak, &t.r90}) v->assign(nf * ns, NAN);
    t.flux.assign(nf * ns * t.nr, NAN); t.npix.assign(nf * ns * t.nr, NAN);
    auto get = [&](const std::vector<std::string>& r, const char* k) {
        auto it = col.find(k);
        return it == col.end() ? NAN : num(r[it->second]);
    };
    for (const auto& r : rows) {
        const size_t f = std::stoul(r[col["frame"]]), s = std::stoul(r[col["star"]]);
        t.file[f] = r[col["file"]];
        if (col.count("date_obs")) t.date_obs[f] = r[col["date_obs"]];
        t.exptime[f] = get(r, "exptime"); t.temp[f] = get(r, "ccd_temp");
        t.bjd[f] = get(r, "bjd_tdb"); t.airmass[f] = get(r, "airmass");
        if (r[col["valid"]] != "1") continue;
        const size_t k = t.fs(f, s);
        t.x[k] = get(r, "x"); t.y[k] = get(r, "y"); t.bg[k] = get(r, "bg"); t.bg_sd[k] = get(r, "bg_sd");
        t.peak[k] = get(r, "peak_raw"); t.r90[k] = get(r, "r90");
        for (size_t i = 0; i < t.nr; ++i) {
            t.flux[t.fsr(f, s, i)] = num(r[fcol[i]]);
            if (i < ncol.size()) t.npix[t.fsr(f, s, i)] = num(r[ncol[i]]);
        }
    }
    return t;
}

}  // namespace transit
