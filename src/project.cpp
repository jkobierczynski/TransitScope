// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/project.h"

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <type_traits>
#include <ostream>
#include <set>
#include <sstream>

#define TOML_EXCEPTIONS 1
#include <toml.hpp>

#include "transit/calib.h"
#include "transit/diagnostics.h"
#include "transit/provenance.h"

namespace fs = std::filesystem;

namespace transit {

// --------------------------------------------------------------------------------------------
// paths

std::string Project::resolve(const std::string& p) const {
    if (p.empty()) return p;
    // forward slashes on every system (Windows accepts them), so paths compare and print alike
    if (p == "~" || p.rfind("~/", 0) == 0) {
        const char* home = std::getenv("HOME");
        if (!home || !*home) home = std::getenv("USERPROFILE");
        if (home) return (fs::path(home) / p.substr(p.size() > 1 ? 2 : 1)).lexically_normal().generic_string();
    }
    const fs::path q(p);
    // "/data/x" counts as given on Windows too (rooted, without a drive letter)
    if (q.is_absolute() || q.has_root_directory() || base_dir.empty()) return q.generic_string();
    return (fs::path(base_dir) / q).lexically_normal().generic_string();
}

std::string Project::output(const std::string& file) const { return (fs::path(resolve(output_dir)) / file).generic_string(); }

// --------------------------------------------------------------------------------------------
// reading

namespace {

const char* mode_name(FitConfig::Mode m) { return m == FitConfig::Mode::Free ? "free" : "constrained"; }

class Reader {
public:
    explicit Reader(const toml::table& root) : root_(root) {
        for (auto&& [k, v] : root_) {
            if (!v.is_table()) error("top level: '" + std::string(k.str()) + "' is not a [section]", v);
        }
    }

    struct Section {
        Reader* r;
        std::string name;
        const toml::table* t;   // null if the section is absent
        std::set<std::string> used;

        const toml::node* get(const std::string& key) {
            used.insert(key);
            return t ? t->get(key) : nullptr;
        }
        bool has(const std::string& key) const { return t && t->contains(key); }
        std::string where(const std::string& key) const { return "[" + name + "] " + key; }

        template <class T>
        bool number(const std::string& key, T& out, bool required = false) {
            const toml::node* n = get(key);
            if (!n) { if (required) r->error(where(key) + " is required"); return false; }
            if (auto v = n->value<double>(); v && (n->is_number())) {
                if constexpr (std::is_integral_v<T>) {
                    if (!n->is_integer() || *v < 0) { r->error(where(key) + " must be a whole number >= 0", *n); return false; }
                }
                out = static_cast<T>(*v);
                return true;
            }
            r->error(where(key) + " must be a number", *n);
            return false;
        }
        bool string(const std::string& key, std::string& out, bool required = false) {
            const toml::node* n = get(key);
            if (!n) { if (required) r->error(where(key) + " is required"); return false; }
            if (auto v = n->value<std::string>(); v && n->is_string()) { out = *v; return true; }
            r->error(where(key) + " must be a string", *n);
            return false;
        }
        // ra/dec: a sexagesimal string or a number of degrees
        bool angle(const std::string& key, std::string& out, bool required) {
            const toml::node* n = get(key);
            if (!n) { if (required) r->error(where(key) + " is required"); return false; }
            if (n->is_string()) { out = *n->value<std::string>(); return true; }
            if (n->is_number()) {
                char buf[64];
                auto res = std::to_chars(buf, buf + sizeof buf, *n->value<double>());
                out.assign(buf, res.ptr);
                return true;
            }
            r->error(where(key) + " must be a string such as \"20h00m43.713s\" or a number of degrees", *n);
            return false;
        }
        template <class T>
        bool numbers(const std::string& key, std::vector<T>& out, size_t exact = 0) {
            const toml::node* n = get(key);
            if (!n) return false;
            const toml::array* a = n->as_array();
            if (!a) { r->error(where(key) + " must be an array of numbers", *n); return false; }
            std::vector<T> v;
            for (const auto& el : *a) {
                if (!el.is_number() || (std::is_integral_v<T> && (!el.is_integer() || *el.value<double>() < 0))) {
                    r->error(where(key) + (std::is_integral_v<T> ? " must hold whole numbers >= 0" : " must hold numbers"), el);
                    return false;
                }
                v.push_back(static_cast<T>(*el.value<double>()));
            }
            if (exact && v.size() != exact) { r->error(where(key) + " must have " + std::to_string(exact) + " values", *n); return false; }
            out = std::move(v);
            return true;
        }
        bool gaussian(const std::string& key, Gaussian& g) {
            std::vector<double> v;
            if (!numbers(key, v, 2)) return false;
            if (!(v[1] > 0)) { r->error(where(key) + ": the sigma (second value) must be > 0"); return false; }
            g = {v[0], v[1]};
            return true;
        }
        bool stars(const std::string& key, StarList& out) {
            const toml::node* n = get(key);
            if (!n) return false;
            const toml::array* a = n->as_array();
            bool ok = a != nullptr;
            StarList v;
            if (a)
                for (const auto& el : *a) {
                    const toml::array* xy = el.as_array();
                    if (!xy || xy->size() != 2 || !(*xy)[0].is_number() || !(*xy)[1].is_number()) { ok = false; break; }
                    v.emplace_back(*(*xy)[0].value<double>(), *(*xy)[1].value<double>());
                }
            if (!ok) { r->error(where(key) + " must be an array of [x, y] pairs", *n); return false; }
            out = std::move(v);
            return true;
        }
        void finish(const std::set<std::string>& also_allowed = {}) {
            if (!t) return;
            for (auto&& [k, v] : *t) {
                const std::string key(k.str());
                if (!used.count(key) && !also_allowed.count(key)) r->error("unknown key " + where(key), v);
            }
        }
    };

    Section section(const std::string& name) {
        known_.insert(name);
        const toml::node* n = root_.get(name);
        return Section{this, name, n ? n->as_table() : nullptr, {}};
    }

    void error(const std::string& msg) { errors_.push_back(msg); }
    void error(const std::string& msg, const toml::node& n) {
        const auto& src = n.source();
        errors_.push_back(src.begin.line ? msg + " (line " + std::to_string(src.begin.line) + ")" : msg);
    }

    void finish() {
        for (auto&& [k, v] : root_)
            if (!known_.count(std::string(k.str()))) error("unknown section [" + std::string(k.str()) + "]", v);
        if (!errors_.empty()) {
            std::string all = "project file:";
            for (const auto& e : errors_) all += "\n  " + e;
            throw std::runtime_error(all);
        }
    }

private:
    const toml::table& root_;
    std::set<std::string> known_;
    std::vector<std::string> errors_;
};

}  // namespace

Project parse_project(const std::string& text, const std::string& base_dir) {
    toml::table root;
    try {
        root = toml::parse(text);
    } catch (const toml::parse_error& e) {
        std::ostringstream m;
        m << "project file: " << e.description() << " (line " << e.source().begin.line << ", column "
          << e.source().begin.column << ")";
        throw std::runtime_error(m.str());
    }
    Project p;
    p.base_dir = base_dir;
    Reader rd(root);

    {
        auto s = rd.section("project");
        s.string("name", p.name);
        s.string("output_dir", p.output_dir);
        s.finish();
    }
    {
        auto s = rd.section("target");
        s.string("name", p.target_name);
        s.string("planet", p.planet);
        s.angle("ra", p.ra, true);
        s.angle("dec", p.dec, true);
        s.finish();
        try {
            if (!p.ra.empty()) p.reduce.target.ra_deg = parse_ra(p.ra);
            if (!p.dec.empty()) p.reduce.target.dec_deg = parse_dec(p.dec);
        } catch (const std::exception& e) {
            rd.error(std::string("[target] ra/dec: ") + e.what());
        }
    }
    {
        auto s = rd.section("site");
        s.number("lat", p.reduce.site.lat_deg, true);
        s.number("lon", p.reduce.site.lon_deg, true);
        s.number("height_m", p.reduce.site.height_m);
        p.lightcurve.site_height_m = p.reduce.site.height_m;
        s.finish();
    }
    {
        auto s = rd.section("camera");
        s.number("gain_e_per_adu", p.lightcurve.gain);
        s.number("read_noise_e", p.lightcurve.read_noise_e);
        s.number("aperture_cm", p.lightcurve.aperture_cm);
        s.string("filter_name", p.filter_name);
        s.finish();
    }
    {
        auto s = rd.section("frames");
        auto& f = p.frames;
        s.string("backend", f.backend);
        s.string("lights_dir", f.lights_dir, true);
        s.string("lights_prefix", f.lights_prefix);
        s.string("flats_dir", f.flats_dir);
        s.string("bias_dir", f.bias_dir);
        s.string("darks_dir", f.darks_dir);
        s.string("master_bias", f.master_bias);
        s.string("master_dark", f.master_dark);
        s.string("master_flat", f.master_flat);
        s.number("hot_adu", f.hot_adu);
        s.finish();
        if (f.backend != "cfitsio" && f.backend != "xisfconv") rd.error("[frames] backend must be \"cfitsio\" or \"xisfconv\"");
        if (f.master_dark.empty() && f.darks_dir.empty()) rd.error("[frames] needs darks_dir or master_dark");
        if (f.master_flat.empty() && f.flats_dir.empty()) rd.error("[frames] needs flats_dir or master_flat");
        if (f.master_flat.empty() && f.master_bias.empty() && f.bias_dir.empty())
            rd.error("[frames] a master flat is built with the master bias: give bias_dir or master_bias");
        p.reduce.backend = f.backend;
    }
    {
        auto s = rd.section("reduce");
        auto& r = p.reduce;
        auto& ap = r.aperture;
        s.number("ref_index", r.ref_index);
        s.number("target_star", r.target_star);
        s.stars("stars", r.stars);
        s.string("stars_file", p.stars_file);
        if (!r.stars.empty() && !p.stars_file.empty()) rd.error("[reduce] give stars or stars_file, not both");
        s.numbers("radii", ap.radii);
        std::vector<double> ann;
        if (s.numbers("annulus", ann, 2)) { ap.ann_in = ann[0]; ap.ann_out = ann[1]; }
        s.number("centroid_radius", ap.centroid_radius);
        s.number("centroid_iterations", ap.centroid_iter);
        s.number("half_box", ap.half_box);
        s.number("edge_margin", ap.edge_margin);
        s.number("detect_margin", r.detect_margin);
        s.number("threads", r.threads);
        s.finish();
        if (ap.radii.empty()) rd.error("[reduce] radii must not be empty");
        if (!(ap.ann_in < ap.ann_out)) rd.error("[reduce] annulus: inner radius must be below the outer one");
        p.lightcurve.target = r.target_star;
    }
    {
        auto s = rd.section("ephemeris");
        auto& e = p.lightcurve.eph;
        s.number("T0", e.T0, true);
        s.number("P", e.P, true);
        s.number("T14", e.T14, true);
        s.finish();
        p.lightcurve.score_model.per = e.P;
        p.fit.period = e.P;
    }
    {
        auto s = rd.section("lightcurve");
        auto& c = p.lightcurve;
        s.numbers("exclude_stars", c.exclude_stars);
        s.numbers("comparison_stars", c.comps);
        s.numbers("radii", c.radii);
        s.number("min_comparison_stars", c.min_comps);
        s.number("min_improvement", c.min_improvement);
        s.number("oot_margin", c.oot_margin);
        s.number("cooler_tolerance", c.cooler_tol);
        s.number("saturation_adu", c.sat_adu);
        s.number("min_relative_flux", c.min_rel_flux);
        s.number("score_rp", c.score_model.rp);
        s.number("score_a", c.score_model.a);
        s.number("score_inc", c.score_model.inc);
        s.number("score_u1", c.score_model.u1);
        s.number("score_u2", c.score_model.u2);
        s.finish();
    }
    {
        auto s = rd.section("limb_darkening");
        auto& l = p.ld;
        std::string mode = "passband";
        s.string("mode", mode);
        if (mode == "passband") l.mode = LdSettings::Mode::Passband;
        else if (mode == "fixed") l.mode = LdSettings::Mode::Fixed;
        else rd.error("[limb_darkening] mode must be \"passband\" or \"fixed\"");
        s.number("teff", l.star.teff);
        s.number("logg", l.star.logg);
        s.number("feh", l.star.feh);
        s.number("teff_err", l.cfg.star_err.teff);
        s.number("logg_err", l.cfg.star_err.logg);
        s.number("feh_err", l.cfg.star_err.feh);
        s.string("table", l.table);
        s.string("filter", l.filter);
        s.string("qe", l.qe);
        s.number("filter_shift_nm", l.cfg.filter_shift_nm);
        s.number("qe_tilt", l.cfg.qe_tilt);
        s.number("model_floor", l.cfg.model_floor);
        s.number("u1", l.u1.mu);
        s.number("u2", l.u2.mu);
        s.number("u1_sigma", l.u1.sigma);
        s.number("u2_sigma", l.u2.sigma);
        s.finish();
        if (l.mode == LdSettings::Mode::Passband && (s.has("u1") || s.has("u2") || s.has("u1_sigma") || s.has("u2_sigma")))
            rd.error("[limb_darkening] u1/u2/u1_sigma/u2_sigma are for mode = \"fixed\"; with \"passband\" they are computed");
        if (l.mode == LdSettings::Mode::Passband && !s.has("teff"))
            rd.error("[limb_darkening] mode = \"passband\" needs teff, logg and feh of the star");
        if (!(l.u1.sigma > 0 && l.u2.sigma > 0)) rd.error("[limb_darkening] u1_sigma and u2_sigma must be > 0");
        if (l.mode == LdSettings::Mode::Fixed) { p.fit.u1 = l.u1; p.fit.u2 = l.u2; }
    }
    {
        auto s = rd.section("fit");
        auto& f = p.fit;
        std::vector<std::string> modes;
        if (const toml::node* n = s.get("modes")) {
            const toml::array* a = n->as_array();
            if (!a) rd.error("[fit] modes must be an array such as [\"constrained\", \"free\"]", *n);
            else {
                p.fit_modes.clear();
                for (const auto& el : *a) {
                    const auto v = el.value<std::string>();
                    if (v && *v == "constrained") p.fit_modes.push_back(FitConfig::Mode::Constrained);
                    else if (v && *v == "free") p.fit_modes.push_back(FitConfig::Mode::Free);
                    else rd.error("[fit] modes: each must be \"constrained\" or \"free\"", el);
                }
            }
        }
        s.number("steps", f.steps);
        s.number("discard", f.discard);
        s.number("thin", f.thin);
        s.number("walkers", f.walkers);
        s.number("seed", f.seed);
        s.number("threads", f.threads);
        if (s.gaussian("a_prior", f.a)) f.a0 = f.a.mu;
        if (s.gaussian("inc_prior", f.inc)) f.inc0 = f.inc.mu;
        s.number("rp_start", f.rp0);
        s.finish();
        if (f.discard >= f.steps) rd.error("[fit] discard must be below steps");
        if (f.thin == 0) rd.error("[fit] thin must be >= 1");
        if (f.walkers < 20 || f.walkers % 2) rd.error("[fit] walkers must be even and >= 20 (2 x the 10 parameters)");
    }
    rd.finish();
    return p;
}

Project load_project(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    const fs::path abs = fs::absolute(path);
    Project p = parse_project(ss.str(), abs.parent_path().generic_string());
    p.path = abs.generic_string();
    return p;
}

// --------------------------------------------------------------------------------------------
// writing

namespace {

std::string num(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[64];
    auto r = std::to_chars(buf, buf + sizeof buf, v);
    std::string s(buf, r.ptr);
    if (s.find_first_of(".eEn") == std::string::npos) s += ".0";   // keep it a TOML float
    return s;
}
std::string str(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\t': o += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); o += b; }
                else o += c;
        }
    }
    return o + "\"";
}
template <class T>
std::string list(const std::vector<T>& v) {
    std::string o = "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) o += ", ";
        if constexpr (std::is_integral_v<T>) o += std::to_string(v[i]);
        else o += num(v[i]);
    }
    return o + "]";
}
// a number of degrees stays a number, a sexagesimal string a string
std::string angle(const std::string& s) {
    // strtod rather than from_chars: floating-point from_chars is missing from older libc++ (macOS)
    char* end = nullptr;
    const double v = s.empty() ? 0 : std::strtod(s.c_str(), &end);
    return !s.empty() && end == s.c_str() + s.size() ? num(v) : str(s);
}

}  // namespace

std::string project_to_toml(const Project& p) {
    std::ostringstream o;
    const auto& f = p.frames;
    const auto& r = p.reduce;
    const auto& ap = r.aperture;
    const auto& c = p.lightcurve;
    const auto& l = p.ld;
    const auto& ft = p.fit;
    o << "# TransitScope project\n\n";
    o << "[project]\nname = " << str(p.name) << "\noutput_dir = " << str(p.output_dir) << "\n\n";
    o << "[target]\nname = " << str(p.target_name) << "\nplanet = " << str(p.planet) << "\nra = " << angle(p.ra) << "\ndec = " << angle(p.dec) << "\n\n";
    o << "[site]\nlat = " << num(r.site.lat_deg) << "\nlon = " << num(r.site.lon_deg) << "\nheight_m = " << num(r.site.height_m) << "\n\n";
    o << "[camera]\ngain_e_per_adu = " << num(c.gain) << "\nread_noise_e = " << num(c.read_noise_e)
      << "\naperture_cm = " << num(c.aperture_cm) << "\nfilter_name = " << str(p.filter_name) << "\n\n";
    o << "[frames]\nbackend = " << str(f.backend) << "\nlights_dir = " << str(f.lights_dir) << "\nlights_prefix = " << str(f.lights_prefix) << "\n";
    for (auto [k, v] : {std::pair{"flats_dir", &f.flats_dir}, {"bias_dir", &f.bias_dir}, {"darks_dir", &f.darks_dir},
                        {"master_bias", &f.master_bias}, {"master_dark", &f.master_dark}, {"master_flat", &f.master_flat}})
        if (!v->empty()) o << k << " = " << str(*v) << "\n";
    o << "hot_adu = " << num(f.hot_adu) << "\n\n";
    o << "[reduce]\nref_index = " << r.ref_index << "\ntarget_star = " << r.target_star << "\n";
    if (!p.stars_file.empty()) o << "stars_file = " << str(p.stars_file) << "\n";
    else if (!r.stars.empty()) {
        o << "stars = [\n";
        for (const auto& [x, y] : r.stars) o << "  [" << num(x) << ", " << num(y) << "],\n";
        o << "]\n";
    }
    o << "radii = " << list(ap.radii) << "\nannulus = [" << num(ap.ann_in) << ", " << num(ap.ann_out) << "]\n"
      << "centroid_radius = " << num(ap.centroid_radius) << "\ncentroid_iterations = " << ap.centroid_iter
      << "\nhalf_box = " << ap.half_box << "\nedge_margin = " << num(ap.edge_margin) << "\ndetect_margin = " << num(r.detect_margin)
      << "\nthreads = " << r.threads << "\n\n";
    o << "[ephemeris]\nT0 = " << num(c.eph.T0) << "\nP = " << num(c.eph.P) << "\nT14 = " << num(c.eph.T14) << "\n\n";
    o << "[lightcurve]\nexclude_stars = " << list(c.exclude_stars) << "\ncomparison_stars = " << list(c.comps)
      << "\nradii = " << list(c.radii) << "\nmin_comparison_stars = " << c.min_comps << "\nmin_improvement = " << num(c.min_improvement) << "\noot_margin = " << num(c.oot_margin)
      << "\ncooler_tolerance = " << num(c.cooler_tol) << "\nsaturation_adu = " << num(c.sat_adu)
      << "\nmin_relative_flux = " << num(c.min_rel_flux) << "\nscore_rp = " << num(c.score_model.rp)
      << "\nscore_a = " << num(c.score_model.a) << "\nscore_inc = " << num(c.score_model.inc)
      << "\nscore_u1 = " << num(c.score_model.u1) << "\nscore_u2 = " << num(c.score_model.u2) << "\n\n";
    o << "[limb_darkening]\n";
    if (l.mode == LdSettings::Mode::Passband) {
        o << "mode = \"passband\"\nteff = " << num(l.star.teff) << "\nlogg = " << num(l.star.logg) << "\nfeh = " << num(l.star.feh)
          << "\nteff_err = " << num(l.cfg.star_err.teff) << "\nlogg_err = " << num(l.cfg.star_err.logg)
          << "\nfeh_err = " << num(l.cfg.star_err.feh) << "\ntable = " << str(l.table) << "\nfilter = " << str(l.filter)
          << "\nqe = " << str(l.qe) << "\nfilter_shift_nm = " << num(l.cfg.filter_shift_nm) << "\nqe_tilt = " << num(l.cfg.qe_tilt)
          << "\nmodel_floor = " << num(l.cfg.model_floor) << "\n\n";
    } else {
        o << "mode = \"fixed\"\nu1 = " << num(l.u1.mu) << "\nu2 = " << num(l.u2.mu) << "\nu1_sigma = " << num(l.u1.sigma)
          << "\nu2_sigma = " << num(l.u2.sigma) << "\n\n";
    }
    o << "[fit]\nmodes = [";
    for (size_t i = 0; i < p.fit_modes.size(); ++i) o << (i ? ", " : "") << '"' << mode_name(p.fit_modes[i]) << '"';
    o << "]\nsteps = " << ft.steps << "\ndiscard = " << ft.discard << "\nthin = " << ft.thin << "\nwalkers = " << ft.walkers
      << "\nseed = " << ft.seed << "\nthreads = " << ft.threads << "\na_prior = [" << num(ft.a.mu) << ", " << num(ft.a.sigma)
      << "]\ninc_prior = [" << num(ft.inc.mu) << ", " << num(ft.inc.sigma) << "]\nrp_start = " << num(ft.rp0) << "\n";
    return o.str();
}

void save_project(const std::string& path, const Project& p) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream o(tmp);
        if (!o) throw std::runtime_error("cannot write " + tmp);
        o << project_to_toml(p);
        if (!o) throw std::runtime_error("error writing " + tmp);
    }
    // a crash while saving never leaves a half-written project
    std::error_code ec;
    fs::rename(tmp, path, ec);
    if (ec) {   // some Windows runtimes do not replace an existing file
        fs::remove(path);
        fs::rename(tmp, path);
    }
}

std::string project_template() {
    return R"(# TransitScope project file (TOML).  Run with:  transitscope run project.toml
# Relative paths are relative to this file's directory; "~/" is your home directory.
# Every key except the ones marked "required" can be left out; the value shown is then the default
# unless it says "example".

[project]
name = "HD 189733 b 2014-09-27"             # example
output_dir = "out"                          # masters, phot.csv, lc.csv, ld.json, fit_*.json go here

[target]
name = "HD 189733"                          # example
planet = "HD189733b"                        # as ExoClock / ETD name it (for the exports; example)
ra = "20h00m43.713s"                        # required (ICRS), or degrees as a number
dec = "+22d42m39.07s"                       # required

[site]
lat = 50.93                                 # required, degrees north
lon = 5.40                                  # required, degrees east
height_m = 50.0

[camera]                                    # for the light-curve noise model
gain_e_per_adu = 0.26
read_noise_e = 4.0
aperture_cm = 20.0                          # telescope aperture, for scintillation
filter_name = "Astrodon Exo-Planet BB"      # as you choose it in the ExoClock / ETD upload forms (example)

[frames]
backend = "cfitsio"                         # or "xisfconv"
lights_dir = "2014-09-27 HD189733b"         # required (example)
lights_prefix = "HD189733b_"                # light frames: names starting with this (example)
flats_dir = "2014-09-27 HD189733b"          # files named *flat*; in another folder: all FITS files
bias_dir = "bias"                           # files named *bias*, else all FITS files
darks_dir = "darks"                         # files named *dark*, else all FITS files
# master_bias = "masters/master_bias.fits"  # ready-made masters instead of the folders above
# master_dark = "masters/master_dark_10s.fits"
# master_flat = "masters/master_flat.fits"
hot_adu = 50.0                              # hot pixel: master dark - master bias above this (ADU)

[reduce]
ref_index = 99                              # reference frame (0-based): star list and start of tracking (example)
target_star = 0                             # index of the target in the star list
# stars_file = "stars.csv"                  # "x,y" per line on the reference frame ...
# stars = [[627.8, 648.9], [511.1, 979.7]]  # ... or inline; neither: stars are detected on the reference frame
radii = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]   # aperture radii, px
annulus = [26.0, 40.0]                      # background annulus, px
centroid_radius = 16.0
centroid_iterations = 4
half_box = 45
edge_margin = 60.0                          # stars closer to an edge are not measured
detect_margin = 70.0                        # detected stars closer to an edge are dropped
threads = 0                                 # 0 = all cores

[ephemeris]                                 # required
T0 = 2454279.436714                         # BJD_TDB (example: Agol et al. 2010)
P = 2.21857567                              # days
T14 = 0.07527                               # transit duration, days

[lightcurve]
exclude_stars = [10, 32, 36]                # never comparison stars (example)
comparison_stars = []                       # empty: chosen automatically
radii = []                                  # radii to try; empty: all of [reduce] radii
min_comparison_stars = 3
min_improvement = 0.01                      # a star is dropped only if that lowers the residual RMS by more than 1 %
oot_margin = 0.003                          # days beyond T14/2 still treated as in transit
cooler_tolerance = 0.35                     # CCD-TEMP off its median by more: frame flagged
saturation_adu = 45000.0                    # target peak above this: frame flagged
min_relative_flux = 0.03                    # comparison candidates: at least this fraction of the target
score_rp = 0.155                            # transit shape used to score comparison sets
score_a = 8.84
score_inc = 85.7
score_u1 = 0.5
score_u2 = 0.19

[limb_darkening]
mode = "passband"                           # computed for filter x QE; or "fixed" with u1, u2, u1_sigma, u2_sigma
teff = 5050.0                               # required for "passband" (example)
logg = 4.59
feh = -0.03
teff_err = 50.0
logg_err = 0.02
feh_err = 0.08
table = "claret2011_atlas_quad"             # names in data/ld and data/passbands, or file paths
filter = "astrodon_xop_bb"
qe = "sony_icx285al_qe"
filter_shift_nm = 10.0
qe_tilt = 0.1
model_floor = 0.05
# u1 = 0.5                                  # mode = "fixed": the priors themselves (Gaussian mean ...
# u2 = 0.19
# u1_sigma = 0.1                            # ... and sigma)
# u2_sigma = 0.1

[fit]
modes = ["constrained", "free"]             # constrained: a/Rs and inclination priors below
steps = 16000
discard = 2000
thin = 10
walkers = 40
seed = 1
threads = 0
a_prior = [8.863, 0.05]                     # a/Rs: mean, sigma (example: Agol et al. 2010)
inc_prior = [85.71, 0.05]                   # inclination, degrees
rp_start = 0.155                            # walker start for Rp/Rs
)";
}

// --------------------------------------------------------------------------------------------
// running

const char* stage_name(Stage s) {
    switch (s) {
        case Stage::Masters: return "masters";
        case Stage::Reduce: return "reduce";
        case Stage::LightCurve: return "lightcurve";
        case Stage::LimbDarkening: return "ld";
        case Stage::Fit: return "fit";
    }
    return "?";
}

Stage parse_stage(const std::string& s) {
    for (Stage st : {Stage::Masters, Stage::Reduce, Stage::LightCurve, Stage::LimbDarkening, Stage::Fit})
        if (s == stage_name(st)) return st;
    throw std::invalid_argument("unknown stage '" + s + "' (masters, reduce, lightcurve, ld, fit)");
}

namespace {

struct Masters { std::string bias, dark, flat; };

Masters master_paths(const Project& p) {
    const auto& f = p.frames;
    Masters m;
    m.bias = !f.master_bias.empty() ? p.resolve(f.master_bias) : (f.bias_dir.empty() ? "" : p.output("master_bias.fits"));
    m.dark = !f.master_dark.empty() ? p.resolve(f.master_dark) : p.output("master_dark.fits");
    m.flat = !f.master_flat.empty() ? p.resolve(f.master_flat) : p.output("master_flat.fits");
    return m;
}

}  // namespace

std::vector<std::string> check_project_inputs(const Project& p, const RunOptions& opt) {
    std::vector<std::string> bad;
    auto on = [&](Stage s) { return s >= opt.from && s <= opt.to; };
    const auto& f = p.frames;
    const Masters m = master_paths(p);
    const std::string lights = p.resolve(f.lights_dir);
    auto need_file = [&](const std::string& path, const std::string& what) {
        if (!fs::is_regular_file(path)) bad.push_back(what + " not found: " + path);
    };
    auto need_frames = [&](const std::string& dir, const char* word) {
        const std::string d = p.resolve(dir);
        if (!fs::is_directory(d)) bad.push_back(std::string(word) + " folder not found: " + d);
        else if (find_calibration_frames(d, word, lights).empty()) bad.push_back(std::string("no ") + word + " frames in " + d);
    };
    // masters: raw frames to build from, or the given master files
    const bool build = on(Stage::Masters);
    if (on(Stage::Masters) || on(Stage::Reduce)) {
        if (!f.master_bias.empty()) need_file(m.bias, "master bias");
        else if (!f.bias_dir.empty()) { if (build) need_frames(f.bias_dir, "bias"); else need_file(m.bias, "master bias (from the masters stage)"); }
        if (!f.master_dark.empty()) need_file(m.dark, "master dark");
        else if (build) need_frames(f.darks_dir, "dark"); else need_file(m.dark, "master dark (from the masters stage)");
        if (!f.master_flat.empty()) need_file(m.flat, "master flat");
        else if (build) need_frames(f.flats_dir, "flat"); else need_file(m.flat, "master flat (from the masters stage)");
    }
    if (on(Stage::Reduce)) {
        if (!fs::is_directory(lights)) bad.push_back("lights folder not found: " + lights);
        else {
            const auto frames = find_light_frames(lights, f.lights_prefix);
            if (frames.empty()) bad.push_back("no light frames named " + f.lights_prefix + "* in " + lights);
            else if (p.reduce.ref_index >= frames.size())
                bad.push_back("[reduce] ref_index " + std::to_string(p.reduce.ref_index) + " is beyond the last of the " +
                              std::to_string(frames.size()) + " light frames");
        }
        if (!p.stars_file.empty()) need_file(p.resolve(p.stars_file), "[reduce] stars_file");
    }
    // later stages started on their own need the earlier outputs
    if (on(Stage::LightCurve) && !on(Stage::Reduce)) need_file(p.output("phot.csv"), "phot.csv (from the reduce stage)");
    const bool passband = p.ld.mode == LdSettings::Mode::Passband;
    if (on(Stage::LimbDarkening) && passband)
        for (auto [name, sub] : {std::pair{&p.ld.table, "ld"}, {&p.ld.filter, "passbands"}, {&p.ld.qe, "passbands"}})
            try { find_data_file(*name, sub); } catch (const std::exception& e) { bad.push_back(e.what()); }
    if (on(Stage::Fit)) {
        if (!on(Stage::LightCurve)) need_file(p.output("lc.csv"), "lc.csv (from the lightcurve stage)");
        if (passband && !on(Stage::LimbDarkening)) need_file(p.output("ld.json"), "ld.json (from the ld stage)");
    }
    return bad;
}

void run_project(const Project& p, const RunOptions& opt, const StageProgress& progress) {
    if (const auto bad = check_project_inputs(p, opt); !bad.empty()) {
        std::string all = "cannot run the project:";
        for (const auto& b : bad) all += "\n  " + b;
        throw std::runtime_error(all);
    }
    auto log = [&](const std::string& s) { if (opt.log) *opt.log << s << "\n"; };
    auto on = [&](Stage s) { return s >= opt.from && s <= opt.to; };
    auto cb = [&](Stage s) -> ProgressFn {
        if (!progress) return {};
        return [&progress, s](size_t d, size_t t) { return progress(s, d, t); };
    };
    fs::create_directories(p.resolve(p.output_dir));
    const Masters m = master_paths(p);
    const auto& f = p.frames;
    const std::string settings = project_to_toml(p);
    // provenance/<stage>.json: build, times, settings, size + SHA-256 of every input and output
    auto record = [&](Stage s, const std::string& started, const std::vector<std::pair<std::string, std::string>>& in,
                      const std::vector<std::pair<std::string, std::string>>& out) {
        ProvenanceRecord r;
        r.stage = stage_name(s);
        r.started_utc = started;
        r.finished_utc = utc_now_iso();
        r.project_file = p.path;
        r.build = build_info();
        r.settings_toml = settings;
        r.inputs = digest_files(in, p.reduce.threads);
        r.outputs = digest_files(out, p.reduce.threads);
        write_provenance(provenance_path(p, s), r);
        uint64_t bytes = 0;
        for (const auto& d : r.inputs) bytes += d.size;
        std::ostringstream msg;
        msg.precision(2);
        msg << std::fixed << "provenance: " << r.inputs.size() << " inputs (" << bytes / 1e9 << " GB) and "
            << r.outputs.size() << " outputs hashed -> " << provenance_path(p, s);
        log(msg.str());
    };
    using Files = std::vector<std::pair<std::string, std::string>>;

    if (on(Stage::Masters)) {
        const std::string started = utc_now_iso();
        Files in, out;
        const std::string lights = p.resolve(f.lights_dir);
        auto build = [&](const char* type, const std::string& dir, const std::string& outfile, const std::string& bias) {
            const auto files = find_calibration_frames(p.resolve(dir), type, lights);
            if (files.empty()) throw std::runtime_error(std::string("no ") + type + " frames in " + p.resolve(dir));
            const Image img = build_master(f.backend, type, files, bias, cb(Stage::Masters));
            write_master(outfile, img, type, files.size());
            log(std::string("master ") + type + ": " + std::to_string(files.size()) + " frames -> " + outfile);
            for (const auto& fr : files) in.emplace_back(type, fr);
            out.emplace_back(std::string("master_") + type, outfile);
        };
        if (f.master_bias.empty() && !f.bias_dir.empty()) build("bias", f.bias_dir, m.bias, "");
        if (f.master_dark.empty()) build("dark", f.darks_dir, m.dark, "");
        if (f.master_flat.empty()) build("flat", f.flats_dir, m.flat, m.bias);
        if (!out.empty()) record(Stage::Masters, started, in, out);
    }

    if (on(Stage::Reduce)) {
        const std::string started = utc_now_iso();
        ReduceConfig rc = p.reduce;
        rc.lights = find_light_frames(p.resolve(f.lights_dir), f.lights_prefix);
        if (rc.lights.empty())
            throw std::runtime_error("no light frames named " + f.lights_prefix + "* in " + p.resolve(f.lights_dir));
        if (!p.stars_file.empty()) rc.stars = read_star_list(p.resolve(p.stars_file));
        const Calibration cal = load_calibration(f.backend, m.dark, m.flat, m.bias, f.hot_adu);
        const ReduceResult r = reduce_sequence(rc, cal, cb(Stage::Reduce));
        write_phot_csv(p.output("phot.csv"), r);
        log("reduce: " + std::to_string(r.frames.size()) + " frames, " + std::to_string(r.stars.size()) + " stars -> " +
            p.output("phot.csv"));
        Files in;
        for (const auto& l : rc.lights) in.emplace_back("light", l);
        in.emplace_back("master_dark", m.dark);
        in.emplace_back("master_flat", m.flat);
        if (!m.bias.empty() && fs::exists(m.bias)) in.emplace_back("master_bias", m.bias);
        if (!p.stars_file.empty()) in.emplace_back("stars", p.resolve(p.stars_file));
        record(Stage::Reduce, started, in, {{"phot", p.output("phot.csv")}});
    }

    if (on(Stage::LightCurve)) {
        const std::string started = utc_now_iso();
        const PhotTable t = load_phot_csv(p.output("phot.csv"));
        const EnsembleResult lc = build_light_curve(t, p.lightcurve);
        write_light_curve_csv(p.output("lc.csv"), t, lc);
        size_t good = 0;
        for (bool g : lc.good) good += g;
        std::ostringstream s;
        s << "lightcurve: radius " << lc.radius << " px, " << lc.comps.size() << " comparison stars, " << good << "/" << t.nf
          << " good frames -> " << p.output("lc.csv");
        log(s.str());
        record(Stage::LightCurve, started, {{"phot", p.output("phot.csv")}}, {{"lightcurve", p.output("lc.csv")}});
        if (progress && !progress(Stage::LightCurve, 1, 1)) throw Cancelled();
    }

    const bool passband = p.ld.mode == LdSettings::Mode::Passband;
    if (on(Stage::LimbDarkening) && passband) {
        const std::string started = utc_now_iso();
        const std::string tabf = find_data_file(p.ld.table, "ld"), filf = find_data_file(p.ld.filter, "passbands"),
                          qef = find_data_file(p.ld.qe, "passbands");
        ClaretTable tab;
        tab.load(tabf);
        const Curve filter = load_curve_csv(filf);
        const Curve qe = load_curve_csv(qef);
        const LdResult r = passband_limb_darkening(tab, p.ld.star, filter, qe, p.ld.cfg);
        write_ld_json(p.output("ld.json"), r, p.ld.star, p.ld.cfg);
        std::ostringstream s;
        s.precision(3);
        s << std::fixed << "ld: u1 = " << r.u1 << " +- " << r.s1.total << ", u2 = " << r.u2 << " +- " << r.s2.total << " -> "
          << p.output("ld.json");
        log(s.str());
        record(Stage::LimbDarkening, started, {{"ld_table", tabf}, {"filter", filf}, {"qe", qef}}, {{"ld", p.output("ld.json")}});
        if (progress && !progress(Stage::LimbDarkening, 1, 1)) throw Cancelled();
    }

    if (on(Stage::Fit)) {
        const std::string started = utc_now_iso();
        Files in{{"lightcurve", p.output("lc.csv")}}, out;
        if (passband) in.emplace_back("ld", p.output("ld.json"));
        const FitData d = load_fit_data(p.output("lc.csv"));
        // phot.csv supplies the per-frame quantities for the systematics diagnostics
        std::optional<PhotTable> phot;
        if (std::filesystem::exists(p.output("phot.csv"))) {
            phot = load_phot_csv(p.output("phot.csv"));
            in.emplace_back("phot", p.output("phot.csv"));
        }
        for (FitConfig::Mode mode : p.fit_modes) {
            FitConfig fc = p.fit;
            fc.mode = mode;
            if (passband) read_ld_json(p.output("ld.json"), fc);
            const FitResult r = fit_transit(d, fc, cb(Stage::Fit));
            const std::string name = mode_name(mode);
            write_fit_json(p.output("fit_" + name + ".json"), r, fc);
            write_fit_model_csv(p.output("model_" + name + ".csv"), d.t, d.f, r, fc);
            std::ostringstream s;
            s.precision(6);
            s << std::fixed << "fit " << name << ": Tc = " << r.q[0].med << " +- " << 0.5 * (r.q[0].lo + r.q[0].hi)
              << " BJD_TDB, Rp/Rs = " << std::setprecision(4) << r.q[1].med << " -> " << p.output("fit_" + name + ".json");
            log(s.str());
            out.emplace_back("fit_" + name, p.output("fit_" + name + ".json"));
            out.emplace_back("model_" + name, p.output("model_" + name + ".csv"));
            try {
                const FitCurves fcur = fit_curves(d.t, d.f, r, fc);
                const Diagnostics dg = diagnose_fit(phot ? &*phot : nullptr, p.output("lc.csv"), d, fcur.residual,
                                                    p.lightcurve.target);
                write_diagnostics_json(p.output("diagnostics_" + name + ".json"), dg);
                std::istringstream sum(diagnostics_summary(dg));
                for (std::string line; std::getline(sum, line);) log("  " + line);
                out.emplace_back("diagnostics_" + name, p.output("diagnostics_" + name + ".json"));
            } catch (const std::exception& e) {
                log(std::string("  diagnostics skipped: ") + e.what());
            }
        }
        record(Stage::Fit, started, in, out);
    }
}

std::string provenance_path(const Project& p, Stage s) { return p.output(std::string("provenance/") + stage_name(s) + ".json"); }

VerifyReport verify_provenance(const Project& p, bool hash_frames) {
    VerifyReport rep;
    struct Rec { Stage stage; ProvenanceRecord r; };
    std::vector<Rec> recs;
    for (Stage s : {Stage::Masters, Stage::Reduce, Stage::LightCurve, Stage::LimbDarkening, Stage::Fit}) {
        const std::string path = provenance_path(p, s);
        if (!fs::exists(path)) {
            rep.lines.push_back({VerifyReport::Info, std::string(stage_name(s)) + ": no record"});
            continue;
        }
        try {
            recs.push_back({s, read_provenance(path)});
        } catch (const std::exception& e) {
            rep.lines.push_back({VerifyReport::Problem, std::string(stage_name(s)) + ": " + e.what()});
        }
    }
    auto is_frame = [](const std::string& role) { return role == "light" || role == "bias" || role == "dark" || role == "flat"; };
    for (const auto& [stage, r] : recs) {
        const std::string sn = stage_name(stage);
        size_t checked = 0, size_only = 0, bad = 0;
        auto check = [&](const FileDigest& d, bool output) {
            std::error_code ec;
            if (!fs::exists(d.path, ec)) {
                rep.lines.push_back({VerifyReport::Problem, sn + ": " + (output ? "output " : "input ") + d.path + " is missing"});
                ++bad;
                return;
            }
            const uint64_t size = fs::file_size(d.path, ec);
            const bool hash = hash_frames || !is_frame(d.role);
            if (size != d.size || (hash && sha256_file(d.path) != d.sha256)) {
                rep.lines.push_back({VerifyReport::Problem,
                                     sn + ": " + (output ? "output " : "input ") + d.path + " has changed since the stage ran (" +
                                         r.finished_utc + ")" + (output ? "" : ": rerun from " + sn)});
                ++bad;
            }
            ++checked;
            if (!hash) ++size_only;
        };
        for (const auto& d : r.inputs) check(d, false);
        for (const auto& d : r.outputs) check(d, true);
        if (!r.settings_toml.empty()) {
            // only the sections this stage reads count
            static const std::map<Stage, std::vector<std::string>> used = {
                {Stage::Masters, {"frames"}},
                {Stage::Reduce, {"target", "site", "frames", "reduce"}},
                {Stage::LightCurve, {"ephemeris", "lightcurve", "camera", "site", "reduce"}},
                {Stage::LimbDarkening, {"limb_darkening"}},
                {Stage::Fit, {"ephemeris", "limb_darkening", "fit"}}};
            auto section = [](const std::string& toml, const std::string& name) {
                const std::string head = "[" + name + "]\n";
                const size_t a = toml.find(head);
                if (a == std::string::npos) return std::string();
                const size_t b = toml.find("\n[", a + head.size());
                return toml.substr(a, b == std::string::npos ? std::string::npos : b - a);
            };
            const std::string now = project_to_toml(p);
            std::string changed;
            for (const auto& sec : used.at(stage))
                if (section(now, sec) != section(r.settings_toml, sec)) changed += (changed.empty() ? "[" : ", [") + sec + "]";
            if (!changed.empty())
                rep.lines.push_back({VerifyReport::Info, sn + ": settings in " + changed + " have changed since this stage ran"});
        }
        if (!bad)
            rep.lines.push_back({VerifyReport::Ok, sn + ": " + std::to_string(checked) + " files as recorded (" + r.finished_utc +
                                                       ", TransitScope " + r.build.version + ")" +
                                                       (size_only ? ", " + std::to_string(size_only) + " frames by size only" : "")});
    }
    // chain: every input that an earlier stage produced must be the version that stage recorded
    for (size_t i = 0; i < recs.size(); ++i)
        for (size_t j = i + 1; j < recs.size(); ++j)
            for (const auto& in : recs[j].r.inputs)
                for (const auto& out : recs[i].r.outputs)
                    if (in.path == out.path && in.sha256 != out.sha256)
                        rep.lines.push_back({VerifyReport::Problem,
                                             std::string(stage_name(recs[j].stage)) + " used a different " + in.path + " than " +
                                                 stage_name(recs[i].stage) + " last wrote: rerun from " + stage_name(recs[j].stage)});
    return rep;
}

bool VerifyReport::ok() const {
    for (const auto& l : lines) if (l.kind == Problem) return false;
    return true;
}

}  // namespace transit
