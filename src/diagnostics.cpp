// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/diagnostics.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace transit {

namespace {

double mean(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0) / double(v.size()); }

double pop_std(const std::vector<double>& v) {
    const double m = mean(v);
    double s = 0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / double(v.size()));
}

BinRms bin_rms(const std::vector<double>& res, size_t n, double rms1) {
    BinRms b;
    b.n = n;
    b.bins = res.size() / n;
    std::vector<double> bm(b.bins);
    for (size_t j = 0; j < b.bins; ++j)
        bm[j] = std::accumulate(res.begin() + long(j * n), res.begin() + long((j + 1) * n), 0.0) / double(n);
    b.rms = pop_std(bm);
    b.white = rms1 / std::sqrt(double(n)) * std::sqrt(double(b.bins) / double(b.bins - 1));
    b.beta = b.rms / b.white;
    return b;
}

double quantile(std::vector<double> v, double q) {
    std::sort(v.begin(), v.end());
    const double pos = q * double(v.size() - 1);
    const size_t i = size_t(pos);
    return i + 1 < v.size() ? v[i] + (pos - double(i)) * (v[i + 1] - v[i]) : v.back();
}

std::string fmt(const char* f, double a) {
    char buf[64];
    std::snprintf(buf, sizeof buf, f, a);
    return buf;
}

}  // namespace

Diagnostics compute_diagnostics(const std::vector<double>& res, const std::vector<Quantity>& quantities) {
    Diagnostics d;
    d.n = res.size();
    if (d.n < 40) throw std::invalid_argument("diagnostics need at least 40 points");
    const double rms1 = pop_std(res);
    d.rms_ppt = rms1 * 1000;

    // binned RMS, from single points up to bins of n/8 points (at least 8 bins)
    for (size_t n : {1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 60, 80, 100, 150, 200, 300})
        if (res.size() / n >= 8) d.bins.push_back(bin_rms(res, n, rms1));
    // beta exactly as the fit computes it (median over n = 5, 10, 15, 20, 30, at least 1)
    std::vector<double> betas;
    for (size_t n : {5, 10, 15, 20, 30}) betas.push_back(bin_rms(res, n, rms1).beta);
    std::sort(betas.begin(), betas.end());
    d.beta = std::max(1.0, betas[betas.size() / 2]);
    d.beta_max = 1;
    for (const auto& b : d.bins) d.beta_max = std::max(d.beta_max, b.beta);
    d.n_eff = double(d.n) / (d.beta * d.beta);

    // autocorrelation
    const double m = mean(res);
    double c0 = 0;
    for (double x : res) c0 += (x - m) * (x - m);
    const size_t lags = std::min<size_t>(30, d.n / 10);
    for (size_t k = 1; k <= lags; ++k) {
        double c = 0;
        for (size_t i = 0; i + k < d.n; ++i) c += (res[i] - m) * (res[i + k] - m);
        d.acf.push_back(c / c0);
    }
    d.acf_threshold = 3.0 / std::sqrt(double(d.n));

    // correlations with per-frame quantities
    for (const auto& q : quantities) {
        Correlation c;
        c.name = q.name;
        c.unit = q.unit;
        std::vector<double> x, y;
        for (size_t i = 0; i < res.size() && i < q.values.size(); ++i)
            if (std::isfinite(q.values[i]) && std::isfinite(res[i])) { x.push_back(q.values[i]); y.push_back(res[i]); }
        c.n = x.size();
        c.threshold = 3.0 / std::sqrt(std::max(1.0, double(c.n) / (d.beta * d.beta)));
        if (c.n < 20 || pop_std(x) == 0) {
            c.r = NAN;
            d.correlations.push_back(c);
            continue;
        }
        const double mx = mean(x), my = mean(y);
        double sxy = 0, sxx = 0, syy = 0;
        for (size_t i = 0; i < c.n; ++i) {
            sxy += (x[i] - mx) * (y[i] - my);
            sxx += (x[i] - mx) * (x[i] - mx);
            syy += (y[i] - my) * (y[i] - my);
        }
        c.r = sxy / std::sqrt(sxx * syy);
        c.slope_ppt = sxy / sxx * 1000;
        c.span_ppt = c.slope_ppt * (quantile(x, 0.95) - quantile(x, 0.05));
        c.flagged = std::fabs(c.r) > c.threshold;
        d.correlations.push_back(c);
    }

    // findings in words
    if (d.beta > 1.5)
        d.warnings.push_back("time-correlated noise: binned residuals scatter " + fmt("%.2f", d.beta) +
                             "x more than white noise would; the timing error is inflated by this factor");
    else if (d.beta > 1.2)
        d.warnings.push_back("mild time-correlated noise (beta " + fmt("%.2f", d.beta) + ")");
    size_t run = 0;
    while (run < d.acf.size() && d.acf[run] > d.acf_threshold) ++run;
    if (run > 0)
        d.warnings.push_back("neighbouring residuals are correlated (lag-1 autocorrelation " + fmt("%.2f", d.acf[0]) +
                             ", significant up to lag " + std::to_string(run) + ")");
    // flagged quantities, strongest first; quantities that vary together with a stronger one
    // (|r| > 0.7 between them, directly or through another) are named with it, since the residuals
    // cannot tell them apart
    std::vector<size_t> flagged;
    for (size_t i = 0; i < d.correlations.size(); ++i) if (d.correlations[i].flagged) flagged.push_back(i);
    std::sort(flagged.begin(), flagged.end(), [&](size_t a, size_t b) {
        return std::fabs(d.correlations[a].r) / d.correlations[a].threshold > std::fabs(d.correlations[b].r) / d.correlations[b].threshold;
    });
    auto qcorr = [&](size_t a, size_t b) {
        const auto &x = quantities[a].values, &y = quantities[b].values;
        std::vector<double> u, v;
        for (size_t i = 0; i < x.size() && i < y.size() && i < res.size(); ++i)
            if (std::isfinite(x[i]) && std::isfinite(y[i])) { u.push_back(x[i]); v.push_back(y[i]); }
        if (u.size() < 3) return 0.0;
        const double mu = mean(u), mv = mean(v);
        double suv = 0, suu = 0, svv = 0;
        for (size_t i = 0; i < u.size(); ++i) { suv += (u[i] - mu) * (v[i] - mv); suu += (u[i] - mu) * (u[i] - mu); svv += (v[i] - mv) * (v[i] - mv); }
        return suu > 0 && svv > 0 ? suv / std::sqrt(suu * svv) : 0.0;
    };
    std::vector<bool> named(d.correlations.size(), false);
    for (size_t i : flagged) {
        if (named[i]) continue;
        named[i] = true;
        const auto& c = d.correlations[i];
        std::string w = "residuals correlate with " + c.name + " (r = " + fmt("%+.2f", c.r) + ", limit " +
                        fmt("%.2f", c.threshold) + "): " + fmt("%+.2f", c.span_ppt) + " ppt across its observed range";
        // everything connected to i through |r| > 0.7 (in a night most drifts follow time)
        std::string together;
        std::vector<size_t> group{i};
        for (size_t g = 0; g < group.size(); ++g)
            for (size_t j : flagged)
                if (!named[j] && std::fabs(qcorr(group[g], j)) > 0.7) {
                    named[j] = true;
                    group.push_back(j);
                    together += (together.empty() ? "" : ", ") + d.correlations[j].name;
                }
        if (!together.empty())
            w += "; " + together + (group.size() > 2 ? " vary" : " varies") +
                 " together with it through the night, so the residuals cannot tell which one is the cause";
        d.warnings.push_back(w + (group.size() > 1 ? "; consider detrending against one of them" : "; consider detrending against it"));
    }
    return d;
}

std::vector<Quantity> frame_quantities(const PhotTable& t, const std::vector<size_t>& frames, size_t target,
                                       const std::vector<size_t>& comps, double radius) {
    const int ri = t.radius_index(radius);
    std::vector<Quantity> q = {{"airmass", "", {}},       {"sky background", "ADU", {}}, {"image size r90", "px", {}},
                               {"x drift", "px", {}},     {"y drift", "px", {}},         {"CCD temperature", "C", {}},
                               {"comparison flux", "relative", {}}};
    std::vector<double> comp_sum;
    for (size_t f : frames) {
        const size_t k = t.fs(f, target);
        q[0].values.push_back(t.airmass[f]);
        q[1].values.push_back(t.bg[k]);
        q[2].values.push_back(t.r90[k]);
        q[3].values.push_back(t.x[k]);
        q[4].values.push_back(t.y[k]);
        q[5].values.push_back(t.temp[f]);
        double s = 0;
        bool ok = ri >= 0 && !comps.empty();
        for (size_t c : comps) {
            const double v = ri >= 0 ? t.F(f, c, size_t(ri)) : NAN;
            ok &= std::isfinite(v);
            s += v;
        }
        comp_sum.push_back(ok ? s : NAN);
    }
    // drift and transparency relative to their medians
    for (int i : {3, 4}) {
        std::vector<double> v;
        for (double x : q[size_t(i)].values) if (std::isfinite(x)) v.push_back(x);
        if (v.empty()) continue;
        const double med = quantile(v, 0.5);
        for (double& x : q[size_t(i)].values) x -= med;
    }
    std::vector<double> v;
    for (double x : comp_sum) if (std::isfinite(x)) v.push_back(x);
    const double med = v.empty() ? NAN : quantile(v, 0.5);
    for (double x : comp_sum) q[6].values.push_back(x / med);
    return q;
}

Diagnostics diagnose_fit(const PhotTable* phot, const std::string& lc_csv, const FitData& data,
                         const std::vector<double>& residual, size_t target) {
    std::vector<Quantity> q;
    if (phot) {
        std::vector<size_t> comps;
        double radius = NAN;
        read_lc_header(lc_csv, comps, radius);
        for (size_t f : data.frame)
            if (f >= phot->nf) throw std::runtime_error("lc.csv does not belong to this phot.csv (frame " + std::to_string(f) + ")");
        q = frame_quantities(*phot, data.frame, target, comps, radius);
    }
    return compute_diagnostics(residual, q);
}

std::vector<double> load_residuals(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(in, line);
    std::stringstream hs(line);
    int col = -1, i = 0;
    for (std::string h; std::getline(hs, h, ','); ++i) if (h == "residual") col = i;
    if (col < 0) throw std::runtime_error(path + ": no residual column");
    std::vector<double> r;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string v;
        for (int k = 0; k <= col; ++k) std::getline(ss, v, ',');
        r.push_back(std::stod(v));
    }
    return r;
}

void read_lc_header(const std::string& path, std::vector<size_t>& comps, double& radius) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    comps.clear();
    radius = NAN;
    for (std::string line; std::getline(in, line) && !line.empty() && line[0] == '#';) {
        if (line.rfind("# aperture_radius_px=", 0) == 0) radius = std::stod(line.substr(21));
        if (line.rfind("# comparison_stars=", 0) == 0) {
            std::stringstream ss(line.substr(19));
            for (std::string s; std::getline(ss, s, ';');) if (!s.empty()) comps.push_back(std::stoul(s));
        }
    }
}

void write_diagnostics_json(const std::string& path, const Diagnostics& d) {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(6);
    auto num = [](double v) { return std::isfinite(v) ? fmt("%.6g", v) : std::string("null"); };
    auto esc = [](const std::string& s) {
        std::string r;
        for (char c : s) { if (c == '"' || c == '\\') r += '\\'; r += c; }
        return r;
    };
    o << "{\n  \"n\": " << d.n << ",\n  \"rms_ppt\": " << num(d.rms_ppt) << ",\n  \"beta\": " << num(d.beta)
      << ",\n  \"beta_max\": " << num(d.beta_max) << ",\n  \"n_eff\": " << num(d.n_eff) << ",\n  \"bins\": [\n";
    for (size_t i = 0; i < d.bins.size(); ++i) {
        const auto& b = d.bins[i];
        o << "    {\"n\": " << b.n << ", \"bins\": " << b.bins << ", \"rms_ppt\": " << num(b.rms * 1000) << ", \"white_ppt\": "
          << num(b.white * 1000) << ", \"beta\": " << num(b.beta) << "}" << (i + 1 < d.bins.size() ? "," : "") << "\n";
    }
    o << "  ],\n  \"acf_threshold\": " << num(d.acf_threshold) << ",\n  \"acf\": [";
    for (size_t i = 0; i < d.acf.size(); ++i) o << (i ? ", " : "") << num(d.acf[i]);
    o << "],\n  \"correlations\": [\n";
    for (size_t i = 0; i < d.correlations.size(); ++i) {
        const auto& c = d.correlations[i];
        o << "    {\"name\": \"" << esc(c.name) << "\", \"unit\": \"" << esc(c.unit) << "\", \"n\": " << c.n << ", \"r\": "
          << num(c.r) << ", \"threshold\": " << num(c.threshold) << ", \"slope_ppt\": " << num(c.slope_ppt)
          << ", \"span_ppt\": " << num(c.span_ppt) << ", \"flagged\": " << (c.flagged ? "true" : "false") << "}"
          << (i + 1 < d.correlations.size() ? "," : "") << "\n";
    }
    o << "  ],\n  \"warnings\": [";
    for (size_t i = 0; i < d.warnings.size(); ++i) o << (i ? ",\n    " : "\n    ") << "\"" << esc(d.warnings[i]) << "\"";
    o << (d.warnings.empty() ? "" : "\n  ") << "]\n}\n";
}

std::string diagnostics_summary(const Diagnostics& d) {
    std::ostringstream o;
    o << "residuals: rms " << fmt("%.2f", d.rms_ppt) << " ppt, beta " << fmt("%.2f", d.beta) << " (max "
      << fmt("%.2f", d.beta_max) << "), lag-1 autocorrelation " << fmt("%+.2f", d.acf.empty() ? 0 : d.acf[0]) << "\n";
    o << "correlation with:";
    for (const auto& c : d.correlations)
        o << "  " << c.name << " " << (std::isfinite(c.r) ? fmt("%+.2f", c.r) : std::string("n/a")) << (c.flagged ? "(!)" : "");
    o << "\n";
    if (d.warnings.empty()) o << "no significant red noise or systematics found\n";
    for (const auto& w : d.warnings) o << "WARNING: " << w << "\n";
    return o.str();
}

}  // namespace transit
