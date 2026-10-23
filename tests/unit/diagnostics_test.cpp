// SPDX-License-Identifier: GPL-3.0-or-later
// Red-noise and systematics diagnostics on synthetic residuals:
//   white noise          -> beta ~ 1, no autocorrelation, no correlation flagged
//   AR(1) noise (phi .7) -> beta > 1.5, lag-1 autocorrelation ~ .7 and flagged
//   white + airmass term -> airmass flagged with the injected slope, the other quantity not
// and, with a phot.csv/lc.csv argument, beta equals the fit's beta on the same residuals.
#include <cmath>
#include <cstdio>
#include <random>
#include <string>

#include "transit/diagnostics.h"
#include "transit/fit.h"

using namespace transit;

static int failures = 0;
#define CHECK(cond, ...)                                                       \
    do {                                                                       \
        if (!(cond)) { std::printf("FAIL %s: ", #cond); std::printf(__VA_ARGS__); std::printf("\n"); ++failures; } \
    } while (0)

int main(int argc, char** argv) {
    const size_t N = 600;
    std::mt19937_64 rng(42);
    std::normal_distribution<double> g(0, 1);

    std::vector<double> airmass(N), other(N);
    for (size_t i = 0; i < N; ++i) {
        airmass[i] = 1.05 + 0.6 * double(i) / double(N) + 0.002 * g(rng);
        other[i] = g(rng);
    }
    auto quantities = [&] { return std::vector<Quantity>{{"airmass", "", airmass}, {"other", "", other}, {"constant", "", std::vector<double>(N, 3.0)}}; };

    // white noise
    {
        std::vector<double> r(N);
        for (auto& x : r) x = 0.002 * g(rng);
        const Diagnostics d = compute_diagnostics(r, quantities());
        std::printf("white:   beta %.3f  acf1 %+.3f  r(airmass) %+.3f\n", d.beta, d.acf[0], d.correlations[0].r);
        CHECK(d.beta < 1.25, "beta %.3f", d.beta);
        CHECK(std::fabs(d.acf[0]) < d.acf_threshold, "acf1 %.3f", d.acf[0]);
        CHECK(!d.correlations[0].flagged && !d.correlations[1].flagged, "a correlation flagged");
        CHECK(std::isnan(d.correlations[2].r) && !d.correlations[2].flagged, "constant quantity");
        CHECK(std::fabs(d.rms_ppt - 2.0) < 0.2, "rms %.3f", d.rms_ppt);
        CHECK(d.bins.front().n == 1 && std::fabs(d.bins.front().beta - std::sqrt((N - 1.0) / N)) < 1e-12, "bin 1");
        CHECK(N / d.bins.back().n >= 8, "last bin");
    }
    // AR(1)
    {
        std::vector<double> r(N);
        double prev = 0;
        for (auto& x : r) x = prev = 0.7 * prev + 0.002 * g(rng);
        const Diagnostics d = compute_diagnostics(r, quantities());
        std::printf("AR(1):   beta %.3f  acf1 %+.3f  warnings %zu\n", d.beta, d.acf[0], d.warnings.size());
        CHECK(d.beta > 1.5, "beta %.3f", d.beta);
        CHECK(std::fabs(d.acf[0] - 0.7) < 0.1, "acf1 %.3f", d.acf[0]);
        CHECK(d.warnings.size() >= 2, "warnings %zu", d.warnings.size());
        CHECK(std::fabs(d.n_eff - double(N) / (d.beta * d.beta)) < 1e-9, "n_eff");
    }
    // airmass systematic: 3 ppt per unit airmass
    {
        std::vector<double> r(N);
        for (size_t i = 0; i < N; ++i) r[i] = 0.003 * (airmass[i] - 1.3) + 0.0005 * g(rng);
        const Diagnostics d = compute_diagnostics(r, quantities());
        const auto& c = d.correlations[0];
        std::printf("airmass: r %+.3f (limit %.3f) slope %.2f ppt span %.2f ppt; other r %+.3f\n", c.r, c.threshold,
                    c.slope_ppt, c.span_ppt, d.correlations[1].r);
        CHECK(c.flagged && c.r > 0.5, "airmass r %.3f", c.r);
        CHECK(std::fabs(c.slope_ppt - 3.0) < 0.3, "slope %.3f", c.slope_ppt);
        CHECK(!d.correlations[1].flagged, "other flagged");
        // a second quantity that follows airmass is flagged too, but reported in the same warning
        std::vector<Quantity> q2 = quantities();
        std::vector<double> drift(N);
        for (size_t i = 0; i < N; ++i) drift[i] = -2.0 * airmass[i] + 0.3 * g(rng);
        q2.push_back({"drift", "px", drift});
        const Diagnostics d2 = compute_diagnostics(r, q2);
        size_t mentions = 0;
        for (const auto& w : d2.warnings) mentions += w.find("airmass") != std::string::npos;
        std::printf("grouped: drift flagged %d, %zu warnings, airmass in %zu\n", int(d2.correlations[3].flagged), d2.warnings.size(), mentions);
        CHECK(d2.correlations[3].flagged, "drift not flagged");
        CHECK(mentions == 1 && d2.warnings.back().find("drift") != std::string::npos, "airmass and drift not in one warning");
    }
    // too short
    {
        bool threw = false;
        try { compute_diagnostics(std::vector<double>(10, 0.0), {}); } catch (const std::invalid_argument&) { threw = true; }
        CHECK(threw, "short input accepted");
    }
    // beta as the fit computes it, on a real light curve
    if (argc > 1) {
        const FitData data = load_fit_data(argv[1]);
        FitConfig cfg;
        cfg.steps = 1500;
        cfg.discard = 500;
        const FitResult fr = fit_transit(data, cfg);
        const FitCurves fc = fit_curves(data.t, data.f, fr, cfg);
        const Diagnostics d = compute_diagnostics(fc.residual, {});
        std::printf("real:    fit beta %.6f  diagnostics beta %.6f  rms %.4f / %.4f ppt\n", fr.beta, d.beta, fr.rms_resid * 1000, d.rms_ppt);
        CHECK(std::fabs(fr.beta - d.beta) < 1e-6, "beta %.6f vs %.6f", fr.beta, d.beta);
    }
    std::printf(failures ? "%d failures\n" : "diagnostics: all checks passed\n", failures);
    return failures ? 1 : 0;
}
