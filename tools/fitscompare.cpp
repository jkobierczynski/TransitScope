// SPDX-License-Identifier: GPL-3.0-or-later
// fitscompare — open the same file with every compiled-in backend and diff the results.
//   fitscompare FILE...
// Per image: dims, header cards missing on either side (structural cards reported separately,
// libxisfconv leaves them out by design), type/value mismatches, pixel max/mean |difference| and
// NaN mismatches. Exit status 1 if anything differs.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "transit/image_reader.h"

using namespace transit;

// Cards that describe how a FITS HDU stores its data. CFITSIO lists them; libxisfconv keeps them
// out of the keyword list (they are regenerated on write), so their absence is expected.
static bool structural(const std::string& k) {
    static const std::set<std::string> s{"SIMPLE", "BITPIX", "EXTEND", "BZERO", "BSCALE", "PCOUNT", "GCOUNT",
                                         "XTENSION", "END", "BLANK", "TFIELDS", "THEAP", "LONGSTRN", "EXTNAME"};
    for (const char* p : {"NAXIS", "ZIMAGE", "ZNAXIS", "ZBITPIX", "ZCMPTYPE", "ZTILE", "ZNAME", "ZVAL", "ZQUANTIZ",
                          "ZDITHER", "ZSIMPLE", "ZEXTEND", "ZTENSION", "ZPCOUNT", "ZGCOUNT", "ZHECKSUM", "ZDATASUM",
                          "TTYPE", "TFORM"})
        if (k.rfind(p, 0) == 0) return true;
    return s.count(k) > 0;
}

static std::map<std::string, HeaderCard> by_key(const std::vector<HeaderCard>& h) {
    std::map<std::string, HeaderCard> m;
    for (auto& c : h)
        if (!c.key.empty() && c.key != "COMMENT" && c.key != "HISTORY") m.emplace(c.key, c);
    return m;
}

static int compare_image(ImageReader& A, ImageReader& B, int img) {
    int diffs = 0;
    const auto ha = A.header(img), hb = B.header(img);
    const auto ma = by_key(ha), mb = by_key(hb);
    std::vector<std::string> struct_only;
    for (auto& [key, c] : ma) {
        auto it = mb.find(key);
        if (it == mb.end()) {
            if (structural(key)) struct_only.push_back(key);
            else { std::cout << "    missing in " << B.backend() << ": " << key << "\n"; ++diffs; }
            continue;
        }
        if (it->second.type != c.type) { std::cout << "    type mismatch " << key << "\n"; ++diffs; }
        auto va = card_double(ha, key), vb = card_double(hb, key);
        const bool same = (va && vb) ? (*va == *vb) : (card_string(ha, key) == card_string(hb, key));
        if (!same) { std::cout << "    value mismatch " << key << ": [" << c.raw_value << "] vs [" << it->second.raw_value << "]\n"; ++diffs; }
    }
    for (auto& [key, c] : mb)
        if (!ma.count(key)) {
            if (structural(key)) continue;
            std::cout << "    missing in " << A.backend() << ": " << key << "\n"; ++diffs;
        }
    const size_t nca = std::count_if(ha.begin(), ha.end(), [](auto& c) { return c.key == "COMMENT" || c.key == "HISTORY"; });
    const size_t ncb = std::count_if(hb.begin(), hb.end(), [](auto& c) { return c.key == "COMMENT" || c.key == "HISTORY"; });
    if (nca != ncb) { std::cout << "    COMMENT/HISTORY count " << nca << " vs " << ncb << "\n"; ++diffs; }
    if (!struct_only.empty()) {
        std::cout << "    structural cards only in " << A.backend() << " (expected):";
        for (auto& k : struct_only) std::cout << ' ' << k;
        std::cout << "\n";
    }
    if (A.dims(img) != B.dims(img)) { std::cout << "    dims mismatch\n"; return diffs + 1; }
    const Image pa = A.pixels(img), pb = B.pixels(img);
    double mx = 0, sum = 0;
    size_t nan_a = 0, nan_b = 0;
    for (size_t i = 0; i < pa.px.size(); ++i) {
        const bool na = std::isnan(pa.px[i]), nb = std::isnan(pb.px[i]);
        nan_a += na && !nb; nan_b += nb && !na;
        if (na || nb) continue;
        const double d = std::fabs(double(pa.px[i]) - pb.px[i]);
        mx = std::max(mx, d); sum += d;
    }
    std::cout << "    pixels: max|d|=" << mx << " mean|d|=" << sum / pa.px.size() << "  NaN only in " << A.backend()
              << "=" << nan_a << ", only in " << B.backend() << "=" << nan_b << "\n";
    return diffs + (mx > 0 || nan_a || nan_b);
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: fitscompare FILE...\n"; return 2; }
    const auto backends = available_backends();
    int failures = 0;
    for (int f = 1; f < argc; ++f) {
        const std::string path = argv[f];
        std::cout << "== " << path << "\n";
        std::vector<std::unique_ptr<ImageReader>> rd;
        for (const auto& b : backends) {
            auto r = make_reader(b);
            try { r->open(path); rd.push_back(std::move(r)); }
            catch (const std::exception& e) { std::cout << "  [" << b << "] open failed: " << e.what() << "\n"; ++failures; }
        }
        if (rd.empty()) continue;
        for (auto& r : rd) std::cout << "  [" << r->backend() << "] images=" << r->hdu_count() << "\n";
        for (size_t k = 1; k < rd.size(); ++k) {
            auto &A = *rd[0], &B = *rd[k];
            if (A.hdu_count() != B.hdu_count()) { std::cout << "  image count differs\n"; ++failures; }
            for (int i = 0; i < std::min(A.hdu_count(), B.hdu_count()); ++i) {
                std::cout << "  image " << i << " (dims";
                for (long d : A.dims(i)) std::cout << ' ' << d;
                std::cout << ")\n";
                try { failures += compare_image(A, B, i) ? 1 : 0; }
                catch (const std::exception& e) { std::cout << "    error: " << e.what() << "\n"; ++failures; }
            }
        }
    }
    std::cout << (failures ? "DIFFERENCES FOUND" : "ALL IDENTICAL") << " (" << failures << ")\n";
    return failures ? 1 : 0;
}
