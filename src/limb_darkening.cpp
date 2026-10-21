// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/limb_darkening.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace transit {

namespace {

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream is(s);
    while (std::getline(is, cur, sep)) out.push_back(cur);
    return out;
}

bool numeric_start(const std::string& s) {
    const auto p = s.find_first_not_of(" \t");
    return p != std::string::npos && (std::isdigit(static_cast<unsigned char>(s[p])) || s[p] == '-' || s[p] == '+' || s[p] == '.');
}

// index i and fraction f with g[i] <= v <= g[i+1]; throws outside the grid
std::pair<size_t, double> bracket(const std::vector<double>& g, double v, const char* what) {
    if (g.size() < 2 || v < g.front() - 1e-9 || v > g.back() + 1e-9)
        throw std::out_of_range(std::string("limb darkening: ") + what + " outside the table grid");
    size_t i = std::upper_bound(g.begin(), g.end(), v) - g.begin();
    i = std::clamp<size_t>(i, 1, g.size() - 1) - 1;
    return {i, (v - g[i]) / (g[i + 1] - g[i])};
}

double interp_xy(const std::vector<double>& x, const std::vector<double>& y, double v) {  // flat beyond the ends
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    const size_t i = std::upper_bound(x.begin(), x.end(), v) - x.begin() - 1;
    const double f = (v - x[i]) / (x[i + 1] - x[i]);
    return y[i] + f * (y[i + 1] - y[i]);
}

// photons per unit wavelength from a blackbody: B_lambda * lambda  (arbitrary scale)
double planck_photons(double nm, double teff) {
    const double l = nm * 1e-9;
    return std::pow(l, -4.0) / std::expm1(1.438776877e-2 / (l * teff));
}

}  // namespace

double Curve::at(double lambda_nm) const {
    if (nm.empty() || lambda_nm < nm.front() || lambda_nm > nm.back()) return 0;
    return interp_xy(nm, y, lambda_nm);
}

std::string find_data_file(const std::string& name, const std::string& sub) {
    namespace fs = std::filesystem;
    if (fs::exists(name)) return name;
    std::vector<std::string> bases;
    if (const char* env = std::getenv("TRANSIT_DATA_DIR")) bases.emplace_back(env);
#ifdef TRANSIT_DATA_DIR
    bases.emplace_back(TRANSIT_DATA_DIR);
#endif
    bases.emplace_back("data");
    std::string looked;
    for (const auto& base : bases) {
        for (const std::string ext : {"", ".csv"}) {
            const auto p = fs::path(base) / sub / (name + ext);
            if (fs::exists(p)) return p.string();
        }
        looked += (looked.empty() ? "" : ", ") + (fs::path(base) / sub).string();
    }
    throw std::runtime_error("data file not found: " + name + " (looked in " + looked + ")");
}

Curve load_curve_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Curve c;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || !numeric_start(line)) continue;
        auto f = split(line, ',');
        if (f.size() < 2) continue;
        c.nm.push_back(std::stod(f[0]));
        c.y.push_back(std::stod(f[1]));
    }
    if (c.nm.size() < 2) throw std::runtime_error(path + ": fewer than two data rows");
    for (size_t i = 1; i < c.nm.size(); ++i)
        if (c.nm[i] <= c.nm[i - 1]) throw std::runtime_error(path + ": wavelengths must increase");
    return c;
}

void ClaretTable::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    struct Row { std::string key; double feh, teff, logg, a, b; };
    std::vector<Row> rows;
    std::set<double> fe, te, lg;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto f = split(line, ',');
        if (f.size() != 7 || f[0] == "method") continue;
        Row r{f[0] + f[1], std::stod(f[2]), std::stod(f[3]), std::stod(f[4]), std::stod(f[5]), std::stod(f[6])};
        fe.insert(r.feh); te.insert(r.teff); lg.insert(r.logg);
        rows.push_back(std::move(r));
    }
    if (rows.empty()) throw std::runtime_error(path + ": no limb-darkening rows");
    feh_.assign(fe.begin(), fe.end());
    teff_.assign(te.begin(), te.end());
    logg_.assign(lg.begin(), lg.end());
    const size_t n = feh_.size() * teff_.size() * logg_.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto pos = [](const std::vector<double>& g, double v) { return std::lower_bound(g.begin(), g.end(), v) - g.begin(); };
    for (const auto& r : rows) {
        auto& A = a_[r.key];
        auto& B = b_[r.key];
        if (A.empty()) { A.assign(n, nan); B.assign(n, nan); }
        const size_t k = idx(pos(feh_, r.feh), pos(teff_, r.teff), pos(logg_, r.logg));
        A[k] = r.a;
        B[k] = r.b;
    }
}

std::vector<std::string> ClaretTable::bands() const {
    std::set<std::string> s;
    for (const auto& [k, v] : a_) s.insert(k.substr(1));
    return {s.begin(), s.end()};
}

std::pair<double, double> ClaretTable::at(char method, const std::string& band, const StarParams& s) const {
    const auto key = std::string(1, method) + band;
    const auto ia = a_.find(key);
    if (ia == a_.end()) throw std::invalid_argument("limb darkening: no table for method " + key.substr(0, 1) + ", band " + band);
    const auto& A = ia->second;
    const auto& B = b_.at(key);
    const auto [i, fi] = bracket(feh_, s.feh, "[Fe/H]");
    const auto [j, fj] = bracket(teff_, s.teff, "Teff");
    const auto [k, fk] = bracket(logg_, s.logg, "log g");
    double a = 0, b = 0;
    for (int di = 0; di < 2; ++di)
        for (int dj = 0; dj < 2; ++dj)
            for (int dk = 0; dk < 2; ++dk) {
                const double w = (di ? fi : 1 - fi) * (dj ? fj : 1 - fj) * (dk ? fk : 1 - fk);
                if (w == 0) continue;
                const size_t n = idx(i + di, j + dj, k + dk);
                if (std::isnan(A[n])) throw std::out_of_range("limb darkening: table has no entry next to this star (band " + band + ")");
                a += w * A[n];
                b += w * B[n];
            }
    return {a, b};
}

double band_wavelength_nm(const std::string& band) {
    // effective wavelengths for a solar-type spectrum (Bessell 2005, ARA&A 43, 293; Fukugita et al. 1996)
    static const std::map<std::string, double> w = {
        {"U", 366}, {"B", 438}, {"V", 545}, {"R", 641}, {"I", 798}, {"J", 1220}, {"H", 1630}, {"K", 2190},
        {"Sloanu", 356}, {"Sloang", 483}, {"Sloanr", 626}, {"Sloani", 767}, {"Sloanz", 910}};
    const auto it = w.find(band);
    if (it == w.end()) throw std::invalid_argument("limb darkening: unknown band " + band);
    return it->second;
}

namespace {

struct BandSet { std::vector<double> lam, a, b; };

BandSet band_set(const ClaretTable& tab, char method, const StarParams& s, const std::vector<std::string>& bands,
                 const std::string& skip = {}) {
    std::vector<std::pair<double, std::pair<double, double>>> v;
    for (const auto& b : bands)
        if (b != skip) v.push_back({band_wavelength_nm(b), tab.at(method, b, s)});
    std::sort(v.begin(), v.end(), [](auto& x, auto& y) { return x.first < y.first; });
    BandSet r;
    for (auto& [l, ab] : v) { r.lam.push_back(l); r.a.push_back(ab.first); r.b.push_back(ab.second); }
    return r;
}

struct Avg { double a, b, lam, lo, hi; };

// photon-weighted passband average of the wavelength-interpolated coefficients
Avg passband_average(const BandSet& bs, const Curve& filter, const Curve& qe, double teff, double shift_nm,
                     double tilt, double step) {
    const double l0 = std::max(filter.nm.front() + shift_nm, qe.nm.front());
    const double l1 = std::min(filter.nm.back() + shift_nm, qe.nm.back());
    if (l1 <= l0) throw std::runtime_error("limb darkening: filter and QE curves do not overlap");
    std::vector<double> lam, w;
    for (double l = l0; l <= l1 + 1e-9; l += step) {
        const double ww = filter.at(l - shift_nm) * qe.at(l) * (1 + tilt * (l - 700) / 300) * planck_photons(l, teff);
        lam.push_back(l);
        w.push_back(std::max(ww, 0.0));
    }
    double sw = 0, sa = 0, sb = 0, sl = 0;
    for (size_t i = 0; i < lam.size(); ++i) {
        sw += w[i];
        sa += w[i] * interp_xy(bs.lam, bs.a, lam[i]);
        sb += w[i] * interp_xy(bs.lam, bs.b, lam[i]);
        sl += w[i] * lam[i];
    }
    if (!(sw > 0)) throw std::runtime_error("limb darkening: passband has zero throughput");
    double c = 0, lo = lam.front(), hi = lam.back();
    bool got_lo = false;
    for (size_t i = 0; i < lam.size(); ++i) {
        c += w[i] / sw;
        if (!got_lo && c >= 0.05) { lo = lam[i]; got_lo = true; }
        if (c >= 0.95) { hi = lam[i]; break; }
    }
    return {sa / sw, sb / sw, sl / sw, lo, hi};
}

}  // namespace

LdResult passband_limb_darkening(const ClaretTable& tab, const StarParams& star, const Curve& filter, const Curve& qe,
                                 const LdConfig& cfg) {
    LdResult r;
    auto run = [&](char m, const StarParams& s, double shift, double tilt) {
        return passband_average(band_set(tab, m, s, cfg.bands), filter, qe, s.teff, shift, tilt, cfg.grid_step_nm);
    };
    const Avg L = run('L', star, 0, 0), F = run('F', star, 0, 0);
    r.u1_L = L.a; r.u2_L = L.b; r.u1_F = F.a; r.u2_F = F.b;
    r.u1 = 0.5 * (L.a + F.a);
    r.u2 = 0.5 * (L.b + F.b);
    r.lambda_eff_nm = L.lam;
    r.lambda_lo_nm = L.lo;
    r.lambda_hi_nm = L.hi;

    auto mean_lf = [&](const StarParams& s, double shift, double tilt) {
        const Avg l = run('L', s, shift, tilt), f = run('F', s, shift, tilt);
        return std::pair{0.5 * (l.a + f.a), 0.5 * (l.b + f.b)};
    };

    // 1. fitting method
    r.s1.method = 0.5 * std::abs(L.a - F.a);
    r.s2.method = 0.5 * std::abs(L.b - F.b);

    // 2. wavelength interpolation: leave each interior band out, predict it from the others.
    //    RMS over the interior bands whose wavelength lies inside the passband's 5-95 % range.
    const BandSet all = band_set(tab, 'L', star, cfg.bands);
    double ea = 0, eb = 0;
    int n_in = 0;
    for (const auto& b : cfg.bands) {
        const double lam = band_wavelength_nm(b);
        const auto [aL, bL] = tab.at('L', b, star);
        const auto [aF, bF] = tab.at('F', b, star);
        LdResult::Band rb{b, lam, aL, bL, aF, bF, NAN, NAN, false};
        const bool interior = lam > all.lam.front() && lam < all.lam.back();
        if (interior) {
            const BandSet rest = band_set(tab, 'L', star, cfg.bands, b);
            rb.loo_a = interp_xy(rest.lam, rest.a, lam) - aL;
            rb.loo_b = interp_xy(rest.lam, rest.b, lam) - bL;
            rb.in_passband = lam >= r.lambda_lo_nm && lam <= r.lambda_hi_nm;
            if (rb.in_passband) { ea += rb.loo_a * rb.loo_a; eb += rb.loo_b * rb.loo_b; ++n_in; }
        }
        r.per_band.push_back(rb);
    }
    // leaving a band out doubles the gap, so the error with all bands present is about a quarter of it
    // (linear interpolation error scales with gap^2); keep half as a conservative estimate
    if (n_in) { r.s1.interp = 0.5 * std::sqrt(ea / n_in); r.s2.interp = 0.5 * std::sqrt(eb / n_in); }

    // 3. passband: filter edge +- shift, QE tilt +- qe_tilt (largest excursion of each)
    for (double sh : {-cfg.filter_shift_nm, cfg.filter_shift_nm}) {
        const auto [a, b] = mean_lf(star, sh, 0);
        r.s1.passband = std::max(r.s1.passband, std::abs(a - r.u1));
        r.s2.passband = std::max(r.s2.passband, std::abs(b - r.u2));
    }
    double ta = 0, tb = 0;
    for (double t : {-cfg.qe_tilt, cfg.qe_tilt}) {
        const auto [a, b] = mean_lf(star, 0, t);
        ta = std::max(ta, std::abs(a - r.u1));
        tb = std::max(tb, std::abs(b - r.u2));
    }
    r.s1.passband = std::hypot(r.s1.passband, ta);
    r.s2.passband = std::hypot(r.s2.passband, tb);

    // 4. stellar parameters: symmetric finite differences, summed in quadrature
    auto stellar = [&](StarParams d) {
        StarParams p = star, m = star;
        p.teff += d.teff; p.logg += d.logg; p.feh += d.feh;
        m.teff -= d.teff; m.logg -= d.logg; m.feh -= d.feh;
        const auto [ap, bp] = mean_lf(p, 0, 0);
        const auto [am, bm] = mean_lf(m, 0, 0);
        return std::pair{0.5 * std::abs(ap - am), 0.5 * std::abs(bp - bm)};
    };
    for (StarParams d : {StarParams{cfg.star_err.teff, 0, 0}, StarParams{0, cfg.star_err.logg, 0},
                         StarParams{0, 0, cfg.star_err.feh}}) {
        if (d.teff == 0 && d.logg == 0 && d.feh == 0) continue;
        const auto [da, db] = stellar(d);
        r.s1.stellar = std::hypot(r.s1.stellar, da);
        r.s2.stellar = std::hypot(r.s2.stellar, db);
    }

    // 5. model floor
    r.s1.floor = r.s2.floor = cfg.model_floor;
    for (LdBudget* s : {&r.s1, &r.s2})
        s->total = std::sqrt(s->method * s->method + s->interp * s->interp + s->passband * s->passband +
                             s->stellar * s->stellar + s->floor * s->floor);
    return r;
}

void write_ld_json(const std::string& path, const LdResult& r, const StarParams& star, const LdConfig& cfg) {
    std::ofstream o(path);
    if (!o) throw std::runtime_error("cannot write " + path);
    o.precision(6);
    auto budget = [&](const LdBudget& b) {
        std::ostringstream s;
        s.precision(4);
        s << "{\"method\": " << b.method << ", \"interp\": " << b.interp << ", \"passband\": " << b.passband
          << ", \"stellar\": " << b.stellar << ", \"floor\": " << b.floor << ", \"total\": " << b.total << "}";
        return s.str();
    };
    o << "{\n  \"source\": \"Claret & Bloemen 2011, A&A 529, A75; ATLAS, vt = 2 km/s\",\n"
      << "  \"star\": {\"teff\": " << star.teff << ", \"logg\": " << star.logg << ", \"feh\": " << star.feh
      << ", \"teff_err\": " << cfg.star_err.teff << ", \"logg_err\": " << cfg.star_err.logg
      << ", \"feh_err\": " << cfg.star_err.feh << "},\n"
      << "  \"u1\": " << r.u1 << ", \"u2\": " << r.u2 << ",\n"
      << "  \"u1_sigma\": " << r.s1.total << ", \"u2_sigma\": " << r.s2.total << ",\n"
      << "  \"u1_L\": " << r.u1_L << ", \"u2_L\": " << r.u2_L << ", \"u1_F\": " << r.u1_F << ", \"u2_F\": " << r.u2_F << ",\n"
      << "  \"lambda_eff_nm\": " << r.lambda_eff_nm << ", \"lambda_5pct_nm\": " << r.lambda_lo_nm
      << ", \"lambda_95pct_nm\": " << r.lambda_hi_nm << ",\n"
      << "  \"u1_budget\": " << budget(r.s1) << ",\n  \"u2_budget\": " << budget(r.s2) << ",\n"
      << "  \"bands\": [\n";
    for (size_t i = 0; i < r.per_band.size(); ++i) {
        const auto& b = r.per_band[i];
        o << "    {\"band\": \"" << b.name << "\", \"lambda_nm\": " << b.lambda_nm << ", \"aL\": " << b.aL
          << ", \"bL\": " << b.bL << ", \"aF\": " << b.aF << ", \"bF\": " << b.bF;
        if (!std::isnan(b.loo_a)) o << ", \"loo_da\": " << b.loo_a << ", \"loo_db\": " << b.loo_b;
        o << ", \"in_passband\": " << (b.in_passband ? "true" : "false") << "}" << (i + 1 < r.per_band.size() ? "," : "") << "\n";
    }
    o << "  ]\n}\n";
}

}  // namespace transit
