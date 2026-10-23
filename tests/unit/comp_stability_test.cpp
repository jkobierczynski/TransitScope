// SPDX-License-Identifier: GPL-3.0-or-later
// Stability of the automatic comparison-star selection: excluding a star that was NOT selected
// (as happens when the star list or the exclusions are edited elsewhere) should not change which
// stars are selected. Reports, for every unselected candidate, whether the selection changed.
//
// usage: comp_stability_test PHOT_CSV [max_flips=0]   (2014-09-27 HD 189733 b settings)
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

#include "transit/ensemble.h"

using namespace transit;

static std::string list(const std::vector<size_t>& v) {
    std::string s;
    for (size_t i : v) s += (s.empty() ? "" : ",") + std::to_string(i);
    return s;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: %s PHOT_CSV [max_flips] [min_improvement]\n", argv[0]); return 2; }
    const int max_flips = argc > 2 ? std::atoi(argv[2]) : 0;
    const PhotTable t = load_phot_csv(argv[1]);
    EnsembleConfig cfg;
    cfg.eph = {2454279.436714, 2.21857567, 0.07527};
    cfg.score_model.per = cfg.eph.P;
    cfg.exclude_stars = {10, 32, 36};
    cfg.gain = 0.26; cfg.read_noise_e = 4; cfg.aperture_cm = 20; cfg.site_height_m = 50;
    if (argc > 3) cfg.min_improvement = std::atof(argv[3]);
    std::printf("min_improvement %g\n", cfg.min_improvement);
    const EnsembleResult base = build_light_curve(t, cfg);
    auto sorted = base.comps;
    std::sort(sorted.begin(), sorted.end());
    for (const auto& sc : base.scan) std::printf("  r = %4g: %2zu stars, model rms %.6f\n", sc.radius, sc.ncomps, sc.model_rms);
    std::printf("baseline: r = %g px, %zu comparison stars: %s\n", base.radius, sorted.size(), list(sorted).c_str());
    for (size_t i = 0; i < base.comps.size(); ++i) std::printf("  star %2zu weight %.3f\n", base.comps[i], base.weights[i]);
    int flips = 0, tried = 0;
    for (size_t s = 0; s < t.ns; ++s) {
        if (s == cfg.target || std::count(cfg.exclude_stars.begin(), cfg.exclude_stars.end(), s) ||
            std::count(sorted.begin(), sorted.end(), s))
            continue;
        EnsembleConfig c2 = cfg;
        c2.exclude_stars.push_back(s);
        const EnsembleResult r = build_light_curve(t, c2);
        auto got = r.comps;
        std::sort(got.begin(), got.end());
        ++tried;
        const bool same = got == sorted && r.radius == base.radius;
        if (!same) ++flips;
        std::printf("exclude unselected star %2zu: %s%s\n", s, same ? "same selection" : "CHANGED -> ",
                    same ? "" : (list(got) + " r=" + std::to_string(r.radius)).c_str());
    }
    std::printf("%d of %d exclusions of unselected stars changed the selection\n", flips, tried);

    // tiny changes in the photometry (another backend, reference frame or aperture placement give
    // differences of this order): multiply every flux by (1 + N(0, eps))
    int pflips = 0, ptried = 0, large = 0;
    for (double eps : {1e-4, 3e-4, 1e-3}) {
        for (int seed = 1; seed <= 6; ++seed) {
            PhotTable p = t;
            std::mt19937_64 rng(seed);
            std::normal_distribution<double> N(0, eps);
            for (auto& f : p.flux) f *= 1 + N(rng);
            const EnsembleResult r = build_light_curve(p, cfg);
            auto got = r.comps;
            std::sort(got.begin(), got.end());
            const bool same = got == sorted && r.radius == base.radius;
            // 1e-3 is far beyond real photometry differences: reported, not counted as a failure
            if (eps < 5e-4) { ++ptried; if (!same) ++pflips; }
            else if (!same) ++large;
            std::printf("flux noise %.0e seed %d: %s%s\n", eps, seed, same ? "same selection" : "CHANGED -> ",
                        same ? "" : (list(got) + " r=" + std::to_string(r.radius)).c_str());
        }
    }
    std::printf("%d of %d perturbed photometry runs (noise 1e-4, 3e-4) changed the selection; "
                "%d of 6 with noise 1e-3 (informational)\n", pflips, ptried, large);
    return flips + pflips > max_flips ? 1 : 0;
}
