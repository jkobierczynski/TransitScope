// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/reduce.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <future>
#include <memory>
#include <mutex>
#include <sstream>
#include <tuple>

#include "transit/image_reader.h"
#include "transit/thread_pool.h"

namespace fs = std::filesystem;

namespace transit {

std::vector<std::string> find_light_frames(const std::string& dir, const std::string& prefix) {
    std::vector<std::string> files;
    for (const auto& e : fs::directory_iterator(dir)) {
        const auto name = e.path().filename().string();
        const auto ext = e.path().extension().string();
        if ((ext == ".fit" || ext == ".fits" || ext == ".fts") && name.rfind(prefix, 0) == 0)
            files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

StarList read_star_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    StarList v;
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] == '#') continue;
        double x, y;
        if (std::sscanf(line.c_str(), "%lf,%lf", &x, &y) == 2) v.emplace_back(x, y);
    }
    return v;
}

namespace {
Image load_image(ImageReader& r, const std::string& path) {
    r.open(path);
    return r.pixels(0);
}
}  // namespace

Calibration load_calibration(const std::string& backend, const std::string& dark, const std::string& flat,
                             const std::string& bias, float hot_adu) {
    auto r = make_reader(backend);
    Calibration cal;
    cal.dark = load_image(*r, dark);
    cal.flat = load_image(*r, flat);
    if (!bias.empty()) cal.build_hot_mask(load_image(*r, bias), hot_adu);
    return cal;
}

ReduceResult reduce_sequence(const ReduceConfig& cfg, const Calibration& cal, const ReduceProgress& progress) {
    const auto& files = cfg.lights;
    const size_t n = files.size();
    if (n == 0) throw std::runtime_error("no light frames");
    if (cfg.ref_index >= n) throw std::out_of_range("reference frame index beyond the last frame");
    const ApertureConfig& ap = cfg.aperture;
    const std::string backend = cfg.backend;

    ReduceResult res;
    res.backend = backend;
    res.radii = ap.radii;
    res.ref_index = cfg.ref_index;
    res.target_star = cfg.target_star;
    res.stars = cfg.stars;
    if (res.stars.empty()) {
        auto r = make_reader(backend);
        const Image ref = cal.apply(load_image(*r, files[cfg.ref_index]));
        const double m = cfg.detect_margin;
        for (const auto& d : find_stars(ref))
            if (d.x > m && d.x < ref.nx - m && d.y > m && d.y < ref.ny - m) res.stars.emplace_back(d.x, d.y);
    }
    if (cfg.target_star >= res.stars.size()) throw std::out_of_range("target star index beyond the star list");
    const StarList& stars = res.stars;
    const size_t tgt = cfg.target_star, ref_idx = cfg.ref_index;

    ThreadPool pool(cfg.threads);
    const bool reentrant = backend_thread_safe(backend);
    std::mutex io_mutex;  // only used if the backend library is not reentrant
    struct Frame { Image raw, img; std::string date; double exptime = NAN, temp = NAN; };

    auto load_frame = [&](size_t i) {
        thread_local std::unique_ptr<ImageReader> rd;
        thread_local std::string rd_backend;
        if (!rd || rd_backend != backend) { rd = make_reader(backend); rd_backend = backend; }
        Frame fr;
        {
            std::unique_lock<std::mutex> lock(io_mutex, std::defer_lock);
            if (!reentrant) lock.lock();
            rd->open(files[i]);
            fr.raw = rd->pixels(0);
            const auto hdr = rd->header(0);
            fr.date = card_string(hdr, "DATE-OBS").value_or("");
            fr.exptime = card_double(hdr, "EXPTIME").value_or(NAN);
            fr.temp = card_double(hdr, "CCD-TEMP").value_or(NAN);
        }
        fr.img = cal.apply(fr.raw);
        return fr;
    };
    auto measure_frame = [&](size_t i, const Frame& fr, double dx, double dy) {
        const Image& img = fr.img;
        const FrameTime ft = frame_time(fr.date, fr.exptime, cfg.target, cfg.site);
        FrameResult out;
        out.file = fs::path(files[i]).filename().string();
        out.date_obs = fr.date;
        out.exptime = fr.exptime;
        out.ccd_temp = fr.temp;
        out.bjd_tdb = ft.bjd_tdb_mid;
        out.airmass = ft.airmass;
        out.dx = dx;
        out.dy = dy;
        out.stars.resize(stars.size());
        for (size_t s = 0; s < stars.size(); ++s) {
            double x = stars[s].first + dx, y = stars[s].second + dy;
            const bool inside = x > ap.edge_margin && x < img.nx - ap.edge_margin && y > ap.edge_margin &&
                                y < img.ny - ap.edge_margin;
            if (inside) {
                centroid(img, x, y, ap);
                out.stars[s] = measure(img, fr.raw, x, y, ap);
            }
        }
        return out;
    };

    std::vector<size_t> order;                       // reference, then forward, then backward
    order.push_back(ref_idx);
    for (size_t i = ref_idx + 1; i < n; ++i) order.push_back(i);
    for (size_t i = ref_idx; i-- > 0;) order.push_back(i);
    std::vector<std::pair<double, double>> offset(n, {NAN, NAN});
    std::vector<std::future<Frame>> loads(n);
    std::vector<std::future<FrameResult>> rows(n);
    const size_t lookahead = pool.size() + 2;        // bounds memory: ~2 x lookahead frames in flight
    size_t issued = 0;

    // every task refers to this function's locals: wait for all of them before leaving, also on error
    auto drain = [&] {
        for (auto& f : loads) if (f.valid()) f.wait();
        for (auto& f : rows) if (f.valid()) f.wait();
    };
    try {
        for (size_t pos = 0; pos < n; ++pos) {
            while (issued < n && issued < pos + lookahead) {
                const size_t j = order[issued++];
                loads[j] = pool.submit([&, j] { return load_frame(j); });
            }
            const size_t i = order[pos];
            auto fr = std::make_shared<Frame>(loads[i].get());
            double dx0 = 0, dy0 = 0;
            if (pos > 0) std::tie(dx0, dy0) = offset[i > ref_idx ? i - 1 : i + 1];
            double x = stars[tgt].first + dx0, y = stars[tgt].second + dy0;
            centroid(fr->img, x, y, ap);
            offset[i] = {x - stars[tgt].first, y - stars[tgt].second};
            const auto [dx, dy] = offset[i];
            // the packaged task (and so this lambda's captures) lives as long as its future: release the
            // frame as soon as it has been measured, or every frame stays in memory until the end
            rows[i] = pool.submit([&, i, fr = std::move(fr), dx, dy]() mutable {
                FrameResult r = measure_frame(i, *fr, dx, dy);
                fr.reset();
                return r;
            });
            if (progress && !progress(pos + 1, n)) throw Cancelled();
        }
        res.frames.reserve(n);
        for (auto& f : rows) res.frames.push_back(f.get());
    } catch (...) {
        drain();
        throw;
    }
    return res;
}

void write_phot_csv(const std::string& path, const ReduceResult& r) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out.precision(16);
    out << "frame,file,date_obs,exptime,ccd_temp,bjd_tdb,airmass,star,valid,x,y,bg,bg_sd,peak_raw,r90";
    for (double R : r.radii) out << ",flux_r" << R;
    for (double R : r.radii) out << ",npix_r" << R;
    out << "\n";
    for (size_t i = 0; i < r.frames.size(); ++i) {
        const FrameResult& f = r.frames[i];
        for (size_t s = 0; s < f.stars.size(); ++s) {
            const StarMeasure& m = f.stars[s];
            out << i << ',' << f.file << ',' << f.date_obs << ',' << f.exptime << ',' << f.ccd_temp << ',' << f.bjd_tdb
                << ',' << f.airmass << ',' << s << ',' << (m.valid ? 1 : 0) << ',' << m.x << ',' << m.y << ',' << m.bg
                << ',' << m.bg_sd << ',' << m.peak_raw << ',' << m.r90;
            for (size_t k = 0; k < r.radii.size(); ++k) out << ',' << (m.valid ? m.flux[k] : NAN);
            for (size_t k = 0; k < r.radii.size(); ++k) out << ',' << (m.valid ? m.npix[k] : 0);
            out << '\n';
        }
    }
    if (!out) throw std::runtime_error("error writing " + path);
}

PhotTable to_phot_table(const ReduceResult& r) {
    PhotTable t;
    t.nf = r.frames.size();
    t.ns = r.stars.size();
    t.nr = r.radii.size();
    t.radii = r.radii;
    t.file.resize(t.nf);
    t.date_obs.resize(t.nf);
    for (auto* v : {&t.exptime, &t.temp, &t.bjd, &t.airmass}) v->assign(t.nf, NAN);
    for (auto* v : {&t.x, &t.y, &t.bg, &t.bg_sd, &t.peak, &t.r90}) v->assign(t.nf * t.ns, NAN);
    t.flux.assign(t.nf * t.ns * t.nr, NAN);
    t.npix.assign(t.nf * t.ns * t.nr, NAN);
    for (size_t f = 0; f < t.nf; ++f) {
        const FrameResult& fr = r.frames[f];
        t.file[f] = fr.file;
        t.date_obs[f] = fr.date_obs;
        t.exptime[f] = fr.exptime;
        t.temp[f] = fr.ccd_temp;
        t.bjd[f] = fr.bjd_tdb;
        t.airmass[f] = fr.airmass;
        for (size_t s = 0; s < t.ns && s < fr.stars.size(); ++s) {
            const StarMeasure& m = fr.stars[s];
            if (!m.valid) continue;
            const size_t k = t.fs(f, s);
            t.x[k] = m.x; t.y[k] = m.y; t.bg[k] = m.bg; t.bg_sd[k] = m.bg_sd; t.peak[k] = m.peak_raw; t.r90[k] = m.r90;
            for (size_t i = 0; i < t.nr; ++i) {
                t.flux[t.fsr(f, s, i)] = m.flux[i];
                t.npix[t.fsr(f, s, i)] = m.npix[i];
            }
        }
    }
    return t;
}

}  // namespace transit
