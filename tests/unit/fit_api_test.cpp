// SPDX-License-Identifier: GPL-3.0-or-later
// Checks of the libtransit light-curve and fit API (the calls a GUI makes):
//   - fit_data_from (in memory) gives the points load_fit_data reads from the light-curve CSV
//   - MCMC progress is reported once per step, in order, ending at steps of steps
//   - cancelling from the progress callback throws Cancelled
//   - a progress callback, a cancelled run before, or the thread count do not change the result
//   - samples and fit_curves are consistent with the summary
//
// usage: fit_api_test PHOT_CSV   (a `transitscope reduce` CSV of the 2014-09-27 HD 189733 b set)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "transit/ensemble.h"
#include "transit/fit.h"

using namespace transit;

static int failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

static bool same(const FitResult& a, const FitResult& b) {
    if (a.q.size() != b.q.size() || a.best != b.best || a.samples != b.samples) return false;
    for (size_t i = 0; i < a.q.size(); ++i)
        if (a.q[i].med != b.q[i].med || a.q[i].lo != b.q[i].lo || a.q[i].hi != b.q[i].hi) return false;
    return a.acceptance == b.acceptance;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s PHOT_CSV\n", argv[0]);
        return 2;
    }
    // light curve, in memory and through the CSV
    const PhotTable phot = load_phot_csv(argv[1]);
    EnsembleConfig lc_cfg;
    lc_cfg.eph = {2454279.436714, 2.21857567, 0.07527};
    lc_cfg.score_model.per = lc_cfg.eph.P;
    lc_cfg.exclude_stars = {10, 32, 36};
    lc_cfg.gain = 0.26;
    lc_cfg.read_noise_e = 4;
    lc_cfg.aperture_cm = 20;
    lc_cfg.site_height_m = 50;
    const EnsembleResult lc = build_light_curve(phot, lc_cfg);
    const FitData mem = fit_data_from(phot, lc);
    const std::string csv = "fit_api_test_lc.csv";
    write_light_curve_csv(csv, phot, lc);
    const FitData file = load_fit_data(csv);
    std::remove(csv.c_str());
    bool eq = mem.t.size() == file.t.size() && mem.frame == file.frame && std::fabs(mem.tc_pred - file.tc_pred) < 1e-9;
    for (size_t i = 0; eq && i < mem.t.size(); ++i)
        eq = std::fabs(mem.t[i] - file.t[i]) < 1e-9 && std::fabs(mem.f[i] - file.f[i]) < 1e-14 &&
             std::fabs(mem.sigma[i] - file.sigma[i]) < 1e-15;
    check(eq, "fit_data_from equals load_fit_data of the CSV (" + std::to_string(mem.t.size()) + " points)");

    FitConfig cfg;
    cfg.steps = 1500;
    cfg.discard = 300;
    cfg.thin = 10;
    const FitResult plain = fit_transit(mem, cfg);

    // progress
    std::vector<size_t> seen;
    const FitResult with_progress = fit_transit(mem, cfg, [&](size_t done, size_t total) {
        seen.push_back(done);
        return total == cfg.steps;
    });
    bool in_order = seen.size() == cfg.steps;
    for (size_t i = 0; in_order && i < seen.size(); ++i) in_order = seen[i] == i + 1;
    check(in_order, "progress reported once per step, 1.." + std::to_string(cfg.steps));
    check(same(plain, with_progress), "a progress callback does not change the result");

    // cancellation
    size_t calls = 0;
    bool cancelled = false;
    try {
        fit_transit(mem, cfg, [&](size_t done, size_t) { ++calls; return done < 300; });
    } catch (const Cancelled&) {
        cancelled = true;
    }
    check(cancelled && calls == 300, "cancel at step 300 throws Cancelled");
    check(same(fit_transit(mem, cfg), plain), "a fit after a cancelled one equals a clean fit");

    // thread count
    FitConfig one = cfg;
    one.threads = 1;
    check(same(fit_transit(mem, one), plain), "1 thread gives the same chain as all cores");

    // samples and curves
    const size_t expect = cfg.walkers * ((cfg.steps - cfg.discard) / cfg.thin);
    check(plain.samples.size() == expect && plain.sample_lnprob.size() == expect,
          "samples: walkers x (steps - discard) / thin = " + std::to_string(expect));
    std::vector<double> tc;
    for (const auto& s : plain.samples) tc.push_back(s[0]);
    std::sort(tc.begin(), tc.end());
    const double pos = 0.5 * (tc.size() - 1);
    const double med = tc[size_t(pos)] + (pos - size_t(pos)) * (tc[size_t(pos) + 1] - tc[size_t(pos)]);
    check(med == plain.q[0].med, "median of the Tc samples is the reported Tc");
    const FitCurves c = fit_curves(mem.t, mem.f, plain, cfg);
    double ss = 0, mean = 0;
    for (double r : c.residual) mean += r;
    mean /= c.residual.size();
    for (double r : c.residual) ss += (r - mean) * (r - mean);
    check(c.t.size() == mem.t.size() && std::fabs(std::sqrt(ss / c.residual.size()) - plain.rms_resid) < 1e-15,
          "fit_curves residuals give the reported rms");

    std::printf("Tc = %.6f +- %.6f, Rp/Rs = %.4f (1500 steps)\n", plain.q[0].med, 0.5 * (plain.q[0].lo + plain.q[0].hi),
                plain.q[1].med);
    std::printf("%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
