// SPDX-License-Identifier: GPL-3.0-or-later
// Port of the reference fit.py: quadratic-LD transit x quadratic baseline, noise scale exp(lnf),
// Gaussian LD priors, optional literature priors on a/Rs and inclination, affine-invariant MCMC,
// red-noise beta factor by time-averaging (Winn et al. 2008).
#include "transit/fit.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>
#include <cstdio>
#include <stdexcept>

#include "transit/mcmc.h"
#include "transit/transit_model.h"

namespace transit {
namespace {

enum P { TC, RP, A, INC, U1, U2, C0, C1, C2, LNF, NP };

double stddev(const std::vector<double>& v) {
    const double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    double s = 0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / v.size());
}

struct Model {
    const std::vector<double>& t;
    std::vector<double> tt;
    double period;
    Model(const std::vector<double>& time, double per) : t(time), tt(time.size()), period(per) {
        const double m = std::accumulate(t.begin(), t.end(), 0.0) / t.size();
        for (size_t i = 0; i < t.size(); ++i) tt[i] = t[i] - m;
    }
    void eval(const std::vector<double>& p, std::vector<double>& out) const {
        TransitParams tp{p[TC], period, p[RP], p[A], p[INC], p[U1], p[U2]};
        light_curve(tp, t, out);
        for (size_t i = 0; i < t.size(); ++i) out[i] *= p[C0] + p[C1] * tt[i] + p[C2] * tt[i] * tt[i];
    }
};

}  // namespace

FitResult fit_transit(const std::vector<double>& t, const std::vector<double>& f, const std::vector<double>& s0,
                      double tc0, const FitConfig& cfg, const ProgressFn& progress) {
    if (!std::isfinite(tc0)) throw std::invalid_argument("fit: no predicted Tc");
    if (t.size() != f.size() || t.size() != s0.size() || t.empty()) throw std::invalid_argument("fit: bad input sizes");
    const Model model(t, cfg.period);
    const bool constrained = cfg.mode == FitConfig::Mode::Constrained;

    auto lnprob = [&](const std::vector<double>& p) -> double {
        if (!(std::fabs(p[TC] - tc0) < 0.03 && p[RP] > 0.05 && p[RP] < 0.3 && p[A] > 3 && p[A] < 20 && p[INC] > 75 &&
              p[INC] <= 90 && p[LNF] > -1 && p[LNF] < 2))
            return -INFINITY;
        if (!(p[U1] + p[U2] > 0 && p[U1] + p[U2] < 1 && p[U1] > 0)) return -INFINITY;
        if (p[A] * std::cos(p[INC] * M_PI / 180) > 1 + p[RP]) return -INFINITY;
        double lp = -0.5 * (std::pow((p[U1] - cfg.u1.mu) / cfg.u1.sigma, 2) + std::pow((p[U2] - cfg.u2.mu) / cfg.u2.sigma, 2));
        if (constrained)
            lp += -0.5 * std::pow((p[A] - cfg.a.mu) / cfg.a.sigma, 2) - 0.5 * std::pow((p[INC] - cfg.inc.mu) / cfg.inc.sigma, 2);
        thread_local std::vector<double> m;
        model.eval(p, m);
        const double sc = std::exp(p[LNF]);
        double ll = 0;
        for (size_t i = 0; i < t.size(); ++i) {
            const double s = s0[i] * sc, r = f[i] - m[i];
            ll += r * r / (s * s) + 2 * std::log(s);
        }
        return lp - 0.5 * ll;
    };

    const std::vector<double> p0{tc0, cfg.rp0, cfg.a0, cfg.inc0, cfg.u1.mu, cfg.u2.mu, 1, 0, 0, 0};
    const std::vector<double> sc{5e-4, 2e-3, .05, .05, .02, .02, 1e-4, 1e-3, 1e-2, .01};
    std::mt19937_64 rng(cfg.seed * 7919 + 17);
    std::normal_distribution<double> N(0, 1);
    std::vector<std::vector<double>> init(cfg.walkers, std::vector<double>(NP));
    for (auto& w : init)
        do { for (int d = 0; d < NP; ++d) w[d] = p0[d] + N(rng) * sc[d]; } while (!std::isfinite(lnprob(w)));

    EnsembleSampler sam(cfg.walkers, NP, lnprob, cfg.seed, 2.0, cfg.threads);
    sam.run(init, cfg.steps, progress);
    const auto ch = sam.flat_chain(cfg.discard, cfg.thin);
    const auto lp = sam.flat_lnprob(cfg.discard, cfg.thin);

    FitResult R;
    R.names = {"Tc", "rp", "a", "inc", "u1", "u2", "c0", "c1", "c2", "lnf"};
    R.best = ch[std::max_element(lp.begin(), lp.end()) - lp.begin()];
    R.tc_pred = tc0;
    R.n_points = t.size();
    for (int d = 0; d < NP; ++d) {
        std::vector<double> v;
        for (const auto& s : ch) v.push_back(s[d]);
        std::sort(v.begin(), v.end());
        auto pct = [&](double q) {  // numpy 'linear' percentile
            const double pos = q / 100.0 * (v.size() - 1);
            const size_t i = static_cast<size_t>(pos);
            return i + 1 < v.size() ? v[i] + (pos - i) * (v[i + 1] - v[i]) : v.back();
        };
        const double a = pct(16), m = pct(50), b = pct(84);
        R.q.push_back({m, m - a, b - m});
    }
    std::vector<double> m, res(t.size());
    model.eval(R.best, m);
    for (size_t i = 0; i < t.size(); ++i) res[i] = f[i] - m[i];
    R.rms_resid = stddev(res);
    std::vector<double> betas;
    for (size_t n : {5, 10, 15, 20, 30}) {
        const size_t k = res.size() / n;
        std::vector<double> bm(k);
        for (size_t j = 0; j < k; ++j) bm[j] = std::accumulate(res.begin() + j * n, res.begin() + (j + 1) * n, 0.0) / n;
        betas.push_back(stddev(bm) / (R.rms_resid / std::sqrt(double(n)) * std::sqrt(double(k) / (k - 1))));
    }
    std::sort(betas.begin(), betas.end());
    R.beta = std::max(1.0, betas[betas.size() / 2]);
    R.tc_err_inflated_min = 0.5 * (R.q[TC].lo + R.q[TC].hi) * R.beta * 1440.0;
    std::vector<double> depth, bb;
    for (const auto& s : ch) { depth.push_back(s[RP] * s[RP]); bb.push_back(s[A] * std::cos(s[INC] * M_PI / 180)); }
    auto med = [](std::vector<double> v) { std::sort(v.begin(), v.end()); const size_t n = v.size(); return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]); };
    R.depth_pct = med(depth) * 100;
    R.b = med(bb);
    R.o_c_min = (R.q[TC].med - tc0) * 1440.0;
    R.acceptance = sam.acceptance_fraction();
    const auto tau = sam.autocorr_time(cfg.discard);
    R.autocorr_max = *std::max_element(tau.begin(), tau.end());
    R.samples = ch;
    R.sample_lnprob = lp;
    return R;
}

FitResult fit_transit(const FitData& d, const FitConfig& cfg, const ProgressFn& progress) {
    return fit_transit(d.t, d.f, d.sigma, d.tc_pred, cfg, progress);
}

FitData load_fit_data(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    FitData d;
    for (std::string line; std::getline(in, line);) {
        if (line.rfind("# tc_pred=", 0) == 0) d.tc_pred = std::stod(line.substr(10));
        if (line.empty() || line[0] == '#' || line.rfind("frame", 0) == 0) continue;
        long fr; double bjd, fl, sg, am; int good, oot;
        if (std::sscanf(line.c_str(), "%ld,%lf,%lf,%lf,%lf,%d,%d", &fr, &bjd, &fl, &sg, &am, &good, &oot) == 7 && good) {
            d.frame.push_back(static_cast<size_t>(fr));
            d.t.push_back(bjd);
            d.f.push_back(fl);
            d.sigma.push_back(sg);
        }
    }
    return d;
}

FitData fit_data_from(const PhotTable& phot, const EnsembleResult& lc) {
    FitData d;
    d.tc_pred = lc.tc_pred;
    for (size_t i = 0; i < phot.nf && i < lc.good.size(); ++i) {
        if (!lc.good[i]) continue;
        d.frame.push_back(i);
        d.t.push_back(phot.bjd[i]);
        d.f.push_back(lc.flux[i]);
        d.sigma.push_back(lc.sigma[i]);
    }
    return d;
}

void read_ld_json(const std::string& path, FitConfig& cfg) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string s = ss.str();
    auto num = [&](const std::string& key) {
        const std::string k = "\"" + key + "\":";
        const auto p = s.find(k);
        if (p == std::string::npos) throw std::runtime_error(path + ": no key " + key);
        return std::stod(s.substr(p + k.size()));
    };
    cfg.u1 = {num("u1"), num("u1_sigma")};
    cfg.u2 = {num("u2"), num("u2_sigma")};
}

FitCurves fit_curves(const std::vector<double>& t, const std::vector<double>& f, const FitResult& r, const FitConfig& cfg) {
    const Model model(t, cfg.period);
    std::vector<double> m;
    model.eval(r.best, m);
    FitCurves c;
    c.t = t;
    c.flux = f;
    for (size_t i = 0; i < t.size(); ++i) {
        const double base = r.best[C0] + r.best[C1] * model.tt[i] + r.best[C2] * model.tt[i] * model.tt[i];
        c.flux_detrended.push_back(f[i] / base);
        c.model_transit.push_back(m[i] / base);
        c.residual.push_back(f[i] - m[i]);
    }
    return c;
}

void write_fit_json(const std::string& path, const FitResult& r, const FitConfig& cfg) {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(15);
    o << "{\n  \"mode\": \"" << (cfg.mode == FitConfig::Mode::Free ? "free" : "constrained") << "\",\n";
    for (size_t i = 0; i < r.names.size(); ++i)
        o << "  \"" << r.names[i] << "\": {\"med\": " << r.q[i].med << ", \"lo\": " << r.q[i].lo << ", \"hi\": " << r.q[i].hi << "},\n";
    o << "  \"best\": [";
    for (size_t i = 0; i < r.best.size(); ++i) o << (i ? ", " : "") << r.best[i];
    o << "],\n  \"beta\": " << r.beta << ",\n  \"Tc_err_inflated_min\": " << r.tc_err_inflated_min
      << ",\n  \"rms_resid\": " << r.rms_resid << ",\n  \"depth_pct\": " << r.depth_pct << ",\n  \"b\": " << r.b
      << ",\n  \"O_C_min\": " << r.o_c_min << ",\n  \"tc_pred\": " << r.tc_pred << ",\n  \"acceptance\": " << r.acceptance
      << ",\n  \"autocorr_max\": " << r.autocorr_max << ",\n  \"n_points\": " << r.n_points << ",\n  \"walkers\": "
      << cfg.walkers << ",\n  \"steps\": " << cfg.steps << ",\n  \"discard\": " << cfg.discard << ",\n  \"thin\": " << cfg.thin
      << ",\n  \"seed\": " << cfg.seed << "\n}\n";
}

void write_fit_model_csv(const std::string& path, const std::vector<double>& t, const std::vector<double>& f,
                         const FitResult& r, const FitConfig& cfg) {
    const FitCurves c = fit_curves(t, f, r, cfg);
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(12);
    o << "bjd_tdb,flux,flux_detrended,model_transit,residual\n";
    for (size_t i = 0; i < t.size(); ++i)
        o << c.t[i] << ',' << c.flux[i] << ',' << c.flux_detrended[i] << ',' << c.model_transit[i] << ',' << c.residual[i] << '\n';
}

}  // namespace transit
