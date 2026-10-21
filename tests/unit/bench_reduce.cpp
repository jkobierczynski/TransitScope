// SPDX-License-Identifier: GPL-3.0-or-later
// Times the stages of `reduce` per frame: bench_reduce BACKEND DARK FLAT BIAS STARS FILE...
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include "transit/calib.h"
#include "transit/image_reader.h"
#include "transit/photometry.h"
using namespace transit;
using clk = std::chrono::steady_clock;
int main(int argc, char** argv) {
    if (argc < 7) return 2;
    auto r = make_reader(argv[1]);
    auto load = [&](const char* p) { r->open(p); return r->pixels(0); };
    Calibration cal; cal.dark = load(argv[2]); cal.flat = load(argv[3]); cal.build_hot_mask(load(argv[4]), 50);
    std::vector<std::pair<double, double>> st;
    std::ifstream in(argv[5]); std::string l;
    while (std::getline(in, l)) { double x, y; if (l[0] != '#' && std::sscanf(l.c_str(), "%lf,%lf", &x, &y) == 2) st.emplace_back(x, y); }
    ApertureConfig cfg;
    double t_read = 0, t_cal = 0, t_cen = 0, t_meas = 0, t_tgt = 0; int n = 0;
    auto sec = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); };
    double dx = 0, dy = 0;
    for (int i = 6; i < argc; ++i, ++n) {
        auto t0 = clk::now(); r->open(argv[i]); Image raw = r->pixels(0); auto hdr = r->header(0);
        auto t1 = clk::now(); Image img = cal.apply(raw);
        auto t2 = clk::now(); double x = st[0].first + dx, y = st[0].second + dy; centroid(img, x, y, cfg); dx = x - st[0].first; dy = y - st[0].second;
        auto t3 = clk::now();
        std::vector<std::pair<double, double>> pos;
        for (auto& s : st) { double xx = s.first + dx, yy = s.second + dy; if (xx > 60 && xx < img.nx - 60 && yy > 60 && yy < img.ny - 60) { centroid(img, xx, yy, cfg); pos.emplace_back(xx, yy); } }
        auto t4 = clk::now();
        double chk = 0; for (auto& p : pos) chk += measure(img, raw, p.first, p.second, cfg).flux[1];
        auto t5 = clk::now();
        t_read += sec(t0, t1); t_cal += sec(t1, t2); t_tgt += sec(t2, t3); t_cen += sec(t3, t4); t_meas += sec(t4, t5);
        if (chk == -1) std::puts("");
    }
    std::printf("%d frames, ms/frame: read %.2f  calibrate %.2f  target-track %.2f  centroids(%zu stars) %.2f  measure %.2f  total %.2f\n",
                n, 1e3 * t_read / n, 1e3 * t_cal / n, 1e3 * t_tgt / n, st.size(), 1e3 * t_cen / n, 1e3 * t_meas / n,
                1e3 * (t_read + t_cal + t_tgt + t_cen + t_meas) / n);
}
