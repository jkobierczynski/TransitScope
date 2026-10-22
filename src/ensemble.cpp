// SPDX-License-Identifier: GPL-3.0-or-later
// Port of the reference lc.py: quality flags, out-of-transit mask, inverse-variance weighted
// comparison ensemble, comparison selection scored by fixed-geometry transit-model residuals.
#include "transit/ensemble.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <stdexcept>

#include "transit/lsq.h"
#include "transit/stats.h"

namespace transit {
namespace {

double nanmedian(std::vector<double> v) {
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double stddev(const std::vector<double>& v) {
    if (v.empty()) return NAN;
    const double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double s = 0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / v.size());
}

// RMS of y about a quadratic in time over the masked points (time centred on the masked mean).
double detrend_rms(const std::vector<double>& y, const std::vector<double>& t, const std::vector<bool>& m) {
    double tm = 0;
    size_t n = 0;
    for (size_t i = 0; i < t.size(); ++i) if (m[i]) { tm += t[i]; ++n; }
    tm /= n;
    std::vector<double> tc(t.size());
    for (size_t i = 0; i < t.size(); ++i) tc[i] = t[i] - tm;
    const auto c = polyfit(tc, y, m, 2);
    std::vector<double> r;
    for (size_t i = 0; i < t.size(); ++i) if (m[i]) r.push_back(y[i] - polyval(c, tc[i]));
    return stddev(r);
}

struct Ctx {
    const PhotTable& t;
    const EnsembleConfig& cfg;
    std::vector<bool> good, oot, good_oot;
};

std::vector<double> ensemble_lc(const Ctx& c, size_t ri, const std::vector<size_t>& comps, std::vector<double>* wout) {
    const auto& t = c.t;
    const size_t nf = t.nf, nc = comps.size();
    std::vector<std::vector<double>> Fn(nc, std::vector<double>(nf));
    for (size_t j = 0; j < nc; ++j) {
        std::vector<double> sel;
        for (size_t f = 0; f < nf; ++f) if (c.good_oot[f]) sel.push_back(t.F(f, comps[j], ri));
        const double med = nanmedian(sel);
        for (size_t f = 0; f < nf; ++f) Fn[j][f] = t.F(f, comps[j], ri) / med;
    }
    std::vector<double> w(nc);
    for (size_t j = 0; j < nc; ++j) {
        std::vector<double> y(nf);
        for (size_t f = 0; f < nf; ++f) {
            double ref = 0;
            for (size_t k = 0; k < nc; ++k) if (k != j) ref += Fn[k][f];
            ref = nc > 1 ? ref / (nc - 1) : 1.0;
            y[f] = Fn[j][f] / ref;
        }
        const double r = detrend_rms(y, t.bjd, c.good_oot);
        w[j] = 1.0 / (r * r);
    }
    const double ws = std::accumulate(w.begin(), w.end(), 0.0);
    for (auto& x : w) x /= ws;
    std::vector<double> lc(nf);
    for (size_t f = 0; f < nf; ++f) {
        double ref = 0;
        for (size_t j = 0; j < nc; ++j) ref += Fn[j][f] * w[j];
        lc[f] = t.F(f, c.cfg.target, ri) / ref;
    }
    std::vector<double> sel;
    for (size_t f = 0; f < nf; ++f) if (c.good_oot[f]) sel.push_back(lc[f]);
    const double med = nanmedian(sel);
    for (auto& x : lc) x /= med;
    if (wout) *wout = w;
    return lc;
}

// Residual RMS (good frames) after fitting a fixed-geometry transit x quadratic baseline
// with free Tc and Rp/Rs. Time for the baseline is centred on the mean over all frames.
double model_rms(const Ctx& c, const std::vector<double>& lc, double tc_pred) {
    const auto& t = c.t;
    double tmean = 0;
    for (double x : t.bjd) tmean += x;
    tmean /= t.nf;
    std::vector<double> tt(t.nf), model;
    for (size_t f = 0; f < t.nf; ++f) tt[f] = t.bjd[f] - tmean;
    TransitParams tp = c.cfg.score_model;
    auto res = [&](const std::vector<double>& p, std::vector<double>& r) {
        tp.t0 = p[0]; tp.rp = p[1];
        light_curve(tp, t.bjd, model);
        r.clear();
        for (size_t f = 0; f < t.nf; ++f)
            if (c.good[f]) r.push_back(lc[f] - model[f] * (p[2] + p[3] * tt[f] + p[4] * tt[f] * tt[f]));
    };
    const std::vector<double> p0{tc_pred, c.cfg.score_model.rp, 1, 0, 0}, sc{1e-3, 1e-2, 1e-3, 1e-2, 1e-1};
    const auto p = levenberg_marquardt(res, p0, sc);
    std::vector<double> r;
    res(p, r);
    return stddev(r);
}

}  // namespace

EnsembleResult build_light_curve(const PhotTable& t, const EnsembleConfig& cfg) {
    if (cfg.target >= t.ns) throw std::invalid_argument("target index out of range");
    EnsembleResult R;
    // ephemeris -> predicted Tc for this night
    const double bmean = std::accumulate(t.bjd.begin(), t.bjd.end(), 0.0) / t.nf;
    R.epoch = std::lround((bmean - cfg.eph.T0) / cfg.eph.P);
    R.tc_pred = cfg.eph.T0 + R.epoch * cfg.eph.P;

    Ctx c{t, cfg, {}, {}, {}};
    const double tmed = nanmedian(t.temp);
    c.good.assign(t.nf, true); c.oot.assign(t.nf, false); c.good_oot.assign(t.nf, false);
    for (size_t f = 0; f < t.nf; ++f) {
        bool ok = true;
        for (size_t r = 0; r < t.nr; ++r) ok &= std::isfinite(t.F(f, cfg.target, r));
        if (std::isfinite(tmed) && std::isfinite(t.temp[f]) && std::fabs(t.temp[f] - tmed) >= cfg.cooler_tol) { ok = false; ++R.n_flag_cooler; }
        if (!(t.peak[t.fs(f, cfg.target)] < cfg.sat_adu)) { ok = false; ++R.n_flag_saturated; }
        c.good[f] = ok;
        c.oot[f] = std::fabs(t.bjd[f] - R.tc_pred) > cfg.eph.T14 / 2 + cfg.oot_margin;
        c.good_oot[f] = c.good[f] && c.oot[f];
    }

    // candidate comparison stars: measured in every frame, bright enough, not excluded
    const int r_ref = std::max(0, t.radius_index(12) >= 0 ? t.radius_index(12) : static_cast<int>(t.nr / 2));
    std::vector<double> tflux;
    for (size_t f = 0; f < t.nf; ++f) tflux.push_back(t.F(f, cfg.target, r_ref));
    const double tmedflux = nanmedian(tflux);
    std::vector<size_t> cand;
    for (size_t s = 0; s < t.ns; ++s) {
        if (s == cfg.target) continue;
        if (std::find(cfg.exclude_stars.begin(), cfg.exclude_stars.end(), s) != cfg.exclude_stars.end()) continue;
        bool all = true;
        std::vector<double> fl;
        for (size_t f = 0; f < t.nf; ++f) { const double v = t.F(f, s, r_ref); all &= std::isfinite(v); fl.push_back(t.F(f, s, r_ref)); }
        if (all && nanmedian(fl) / tmedflux > cfg.min_rel_flux) cand.push_back(s);
    }
    if (!cfg.comps.empty()) cand = cfg.comps;
    if (cand.size() < cfg.min_comps) throw std::runtime_error("too few comparison stars");

    std::vector<size_t> radii_idx;
    if (cfg.radii.empty()) for (size_t i = 0; i < t.nr; ++i) radii_idx.push_back(i);
    else for (double Rr : cfg.radii) {
        const int i = t.radius_index(Rr);
        if (i < 0) throw std::invalid_argument("radius not in the photometry table");
        radii_idx.push_back(static_cast<size_t>(i));
    }

    double best = INFINITY;
    size_t best_ri = radii_idx[0];
    std::vector<size_t> best_comps;
    for (size_t ri : radii_idx) {
        auto comps = cand;
        double cur = model_rms(c, ensemble_lc(c, ri, comps, nullptr), R.tc_pred);
        while (cfg.comps.empty() && comps.size() > cfg.min_comps) {
            double bmin = INFINITY;
            size_t drop = 0;
            for (size_t k = 0; k < comps.size(); ++k) {
                auto cc = comps;
                cc.erase(cc.begin() + k);
                const double v = model_rms(c, ensemble_lc(c, ri, cc, nullptr), R.tc_pred);
                if (v < bmin) { bmin = v; drop = k; }
            }
            if (bmin < cur * (1.0 - cfg.min_improvement)) { comps.erase(comps.begin() + drop); cur = bmin; }
            else break;
        }
        const auto lc = ensemble_lc(c, ri, comps, nullptr);
        R.scan.push_back({t.radii[ri], comps.size(), cur, detrend_rms(lc, t.bjd, c.good_oot)});
        if (cur < best) { best = cur; best_ri = ri; best_comps = comps; }
    }

    R.radius = t.radii[best_ri];
    R.comps = best_comps;
    R.flux = ensemble_lc(c, best_ri, best_comps, &R.weights);
    R.good = c.good; R.oot = c.oot;

    // noise model: photon + sky + read noise, plus Young (1967) scintillation
    R.sigma.assign(t.nf, NAN);
    std::vector<double> sph, ssc;
    for (size_t f = 0; f < t.nf; ++f) {
        const double Ft = t.F(f, cfg.target, best_ri) * cfg.gain;
        const double npx = t.npix[t.fsr(f, cfg.target, best_ri)];
        const double bg = t.bg[t.fs(f, cfg.target)];
        const double sp = std::sqrt(Ft + npx * (bg * cfg.gain + cfg.read_noise_e * cfg.read_noise_e)) / Ft;
        const double sc = 0.09 * std::pow(cfg.aperture_cm, -2.0 / 3.0) * std::pow(t.airmass[f], 1.75) *
                          std::exp(-cfg.site_height_m / 8000.0) / std::sqrt(2.0 * t.exptime[f]);
        R.sigma[f] = std::hypot(sp, sc);
        if (std::isfinite(sp)) sph.push_back(sp);
        ssc.push_back(sc);
    }
    R.sig_photon_median = nanmedian(sph);
    R.sig_scint_min = *std::min_element(ssc.begin(), ssc.end());
    R.sig_scint_max = *std::max_element(ssc.begin(), ssc.end());
    return R;
}

void write_light_curve_csv(const std::string& path, const PhotTable& t, const EnsembleResult& r) {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(16);
    o << "# TransitScope light curve\n# aperture_radius_px=" << r.radius << "\n# comparison_stars=";
    for (size_t i = 0; i < r.comps.size(); ++i) o << (i ? ";" : "") << r.comps[i];
    o << "\n# weights=";
    for (size_t i = 0; i < r.weights.size(); ++i) o << (i ? ";" : "") << r.weights[i];
    o << "\n# tc_pred=" << r.tc_pred << " epoch=" << r.epoch << "\n";
    o << "frame,bjd_tdb,flux,sigma,airmass,good,oot\n";
    for (size_t f = 0; f < t.nf; ++f)
        o << f << ',' << t.bjd[f] << ',' << r.flux[f] << ',' << r.sigma[f] << ',' << t.airmass[f] << ','
          << (r.good[f] ? 1 : 0) << ',' << (r.oot[f] ? 1 : 0) << '\n';
}

}  // namespace transit
