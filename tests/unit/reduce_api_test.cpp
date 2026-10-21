// SPDX-License-Identifier: GPL-3.0-or-later
// Checks of the libtransit reduce API (the calls a GUI makes) on a real light sequence:
//   - progress is reported once per frame, in order, ending at n of n
//   - cancelling from the progress callback throws Cancelled, after the queued work is done
//   - a run after a cancelled one gives the same result as a clean run
//   - 1 thread and N threads give identical results
//   - to_phot_table equals load_phot_csv of write_phot_csv (to the 16 digits written)
//
// usage: reduce_api_test LIGHTS_DIR PREFIX DARK FLAT BIAS STARS_CSV REF_INDEX [N_FRAMES]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "transit/reduce.h"

using namespace transit;

static int failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

static bool same(const ReduceResult& a, const ReduceResult& b) {
    if (a.frames.size() != b.frames.size() || a.stars != b.stars) return false;
    for (size_t f = 0; f < a.frames.size(); ++f) {
        const auto &x = a.frames[f], &y = b.frames[f];
        if (x.file != y.file || x.bjd_tdb != y.bjd_tdb || x.dx != y.dx || x.dy != y.dy) return false;
        for (size_t s = 0; s < x.stars.size(); ++s) {
            const auto &p = x.stars[s], &q = y.stars[s];
            if (p.valid != q.valid || p.x != q.x || p.y != q.y || p.flux != q.flux || p.bg != q.bg) return false;
        }
    }
    return true;
}

static bool close_or_both_nan(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return std::fabs(a - b) <= 1e-15 * std::max(1.0, std::fabs(a));
}

int main(int argc, char** argv) {
    if (argc < 8) {
        std::fprintf(stderr, "usage: %s LIGHTS_DIR PREFIX DARK FLAT BIAS STARS_CSV REF_INDEX [N_FRAMES]\n", argv[0]);
        return 2;
    }
    ReduceConfig cfg;
    cfg.lights = find_light_frames(argv[1], argv[2]);
    const size_t nmax = argc > 8 ? std::strtoul(argv[8], nullptr, 10) : 120;
    if (cfg.lights.size() > nmax) cfg.lights.resize(nmax);
    cfg.stars = read_star_list(argv[6]);
    cfg.ref_index = std::strtoul(argv[7], nullptr, 10);
    if (cfg.ref_index >= cfg.lights.size()) cfg.ref_index = cfg.lights.size() / 2;
    cfg.target = {parse_ra("20h00m43.713s"), parse_dec("+22d42m39.07s")};
    cfg.site = {50.93, 5.40, 50};
    const Calibration cal = load_calibration(cfg.backend, argv[3], argv[4], argv[5], 50);
    const size_t n = cfg.lights.size();

    // progress
    std::vector<size_t> seen;
    cfg.threads = 0;
    const ReduceResult full = reduce_sequence(cfg, cal, [&](size_t done, size_t total) {
        seen.push_back(done);
        return total == n;
    });
    bool in_order = seen.size() == n;
    for (size_t i = 0; in_order && i < n; ++i) in_order = seen[i] == i + 1;
    check(in_order, "progress reported once per frame, 1.." + std::to_string(n));
    check(full.frames.size() == n, "one result per frame");

    // cancellation
    size_t calls = 0;
    bool cancelled = false;
    try {
        reduce_sequence(cfg, cal, [&](size_t done, size_t) { ++calls; return done < n / 4; });
    } catch (const Cancelled&) {
        cancelled = true;
    }
    check(cancelled && calls == n / 4, "cancel at frame " + std::to_string(n / 4) + " throws Cancelled");
    check(same(reduce_sequence(cfg, cal), full), "a run after a cancelled one equals a clean run");

    // thread count
    cfg.threads = 1;
    check(same(reduce_sequence(cfg, cal), full), "1 thread gives the same result as all cores");

    // in-memory table vs CSV round trip
    const std::string csv = "reduce_api_test.csv";
    write_phot_csv(csv, full);
    const PhotTable a = to_phot_table(full), b = load_phot_csv(csv);
    bool eq = a.nf == b.nf && a.ns == b.ns && a.nr == b.nr && a.radii == b.radii && a.file == b.file;
    for (size_t i = 0; eq && i < a.flux.size(); ++i) eq = close_or_both_nan(a.flux[i], b.flux[i]) && close_or_both_nan(a.npix[i], b.npix[i]);
    for (size_t i = 0; eq && i < a.x.size(); ++i)
        eq = close_or_both_nan(a.x[i], b.x[i]) && close_or_both_nan(a.y[i], b.y[i]) && close_or_both_nan(a.bg[i], b.bg[i]) &&
             close_or_both_nan(a.r90[i], b.r90[i]) && close_or_both_nan(a.peak[i], b.peak[i]);
    for (size_t f = 0; eq && f < a.nf; ++f) eq = close_or_both_nan(a.bjd[f], b.bjd[f]) && close_or_both_nan(a.airmass[f], b.airmass[f]);
    check(eq, "to_phot_table equals load_phot_csv(write_phot_csv)");
    std::remove(csv.c_str());

    // errors
    bool threw = false;
    try { ReduceConfig bad = cfg; bad.target_star = 999; reduce_sequence(bad, cal); } catch (const std::out_of_range&) { threw = true; }
    check(threw, "target star beyond the star list is rejected");

    std::printf("%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
