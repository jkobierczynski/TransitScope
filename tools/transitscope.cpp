// SPDX-License-Identifier: GPL-3.0-or-later
// transitscope — command-line front end of libtransit.
//
//   transitscope master  --type bias|dark|flat [--bias B.fits] --out M.fits FILE...
//   transitscope reduce  --lights DIR [--prefix NAME_] --dark D.fits --flat F.fits [--bias B.fits --hot-adu 50]
//                       --ra 20h00m43.713s --dec +22d42m39.07s --lat 50.93 --lon 5.40 [--height 50]
//                       [--stars stars.csv | --ref-index N] [--target I] [--radii 8,10,12] [--backend cfitsio]
//                       [--threads N]   (0 = all cores; output is identical for any N)
//                       --out phot.csv
//
// reduce writes one row per (frame, star): positions, background, peak, r90, flux per radius,
// plus BJD_TDB (mid-exposure) and airmass per frame. Later stages (ensemble, fit) read that CSV.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <tuple>
#include <sstream>
#include <string>
#include <vector>

#include "transit/calib.h"
#include "transit/diagnostics.h"
#include "transit/ensemble.h"
#include "transit/export.h"
#include "transit/fit.h"
#include "transit/image_reader.h"
#include "transit/limb_darkening.h"
#include "transit/phot_table.h"
#include "transit/photometry.h"
#include "transit/project.h"
#include "transit/reduce.h"
#include "transit/timing.h"

namespace fs = std::filesystem;
using namespace transit;

namespace {

struct Args {
    std::map<std::string, std::string> opt;
    std::vector<std::string> pos;
    std::string get(const std::string& k, const std::string& def = "") const {
        auto it = opt.find(k);
        return it == opt.end() ? def : it->second;
    }
    std::string need(const std::string& k) const {
        auto it = opt.find(k);
        if (it == opt.end()) throw std::invalid_argument("missing --" + k);
        return it->second;
    }
};

Args parse(int argc, char** argv, int from) {
    Args a;
    for (int i = from; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--", 0) == 0 && i + 1 < argc) a.opt[s.substr(2)] = argv[++i];
        else a.pos.push_back(s);
    }
    return a;
}

int cmd_master(const Args& a) {
    const std::string type = a.need("type"), out = a.need("out");
    const Image m = build_master(a.get("backend", "cfitsio"), type, a.pos, a.get("bias"));
    write_master(out, m, type, a.pos.size());
    std::cerr << "master " << type << ": " << a.pos.size() << " frames -> " << out << "\n";
    return 0;
}

std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v;
    std::stringstream ss(s);
    for (std::string t; std::getline(ss, t, ',');) v.push_back(std::stod(t));
    return v;
}

int cmd_reduce(const Args& a) {
    ReduceConfig cfg;
    cfg.backend = a.get("backend", "cfitsio");
    if (!a.get("radii").empty()) cfg.aperture.radii = parse_list(a.get("radii"));
    cfg.lights = find_light_frames(a.need("lights"), a.get("prefix", ""));
    if (cfg.lights.empty()) throw std::runtime_error("no light frames found");
    cfg.target = {parse_ra(a.need("ra")), parse_dec(a.need("dec"))};
    cfg.site = {std::stod(a.need("lat")), std::stod(a.need("lon")), std::stod(a.get("height", "0"))};
    if (!a.get("stars").empty()) cfg.stars = read_star_list(a.get("stars"));
    cfg.ref_index = std::stoul(a.get("ref-index", "0"));
    cfg.target_star = std::stoul(a.get("target", "0"));
    cfg.threads = static_cast<unsigned>(std::stoul(a.get("threads", "0")));
    const Calibration cal = load_calibration(cfg.backend, a.need("dark"), a.need("flat"), a.get("bias"),
                                             std::stof(a.get("hot-adu", "50")));
    const std::string out = a.need("out");
    const ReduceResult r = reduce_sequence(cfg, cal);
    std::cerr << r.frames.size() << " frames, " << r.stars.size() << " stars, target #" << r.target_star
              << ", backend " << r.backend << "\n";
    write_phot_csv(out, r);
    std::cerr << "wrote " << out << "\n";
    return 0;
}


std::vector<size_t> parse_index_list(const std::string& s) {
    std::vector<size_t> v;
    std::stringstream ss(s);
    for (std::string t; std::getline(ss, t, ',');) if (!t.empty()) v.push_back(std::stoul(t));
    return v;
}

int cmd_lightcurve(const Args& a) {
    const PhotTable t = load_phot_csv(a.need("phot"));
    EnsembleConfig cfg;
    cfg.target = std::stoul(a.get("target", "0"));
    cfg.eph = {std::stod(a.need("T0")), std::stod(a.need("P")), std::stod(a.need("T14"))};
    cfg.oot_margin = std::stod(a.get("oot-margin", "0.003"));
    cfg.cooler_tol = std::stod(a.get("cooler-tol", "0.35"));
    cfg.sat_adu = std::stod(a.get("sat-adu", "45000"));
    cfg.min_rel_flux = std::stod(a.get("min-rel-flux", "0.03"));
    cfg.min_improvement = std::stod(a.get("min-improvement", "0.01"));
    cfg.exclude_stars = parse_index_list(a.get("exclude-stars", ""));
    if (!a.get("radii").empty()) cfg.radii = parse_list(a.get("radii"));
    cfg.comps = parse_index_list(a.get("comps", ""));
    cfg.score_model.per = cfg.eph.P;
    cfg.gain = std::stod(a.get("gain", "1"));
    cfg.read_noise_e = std::stod(a.get("read-noise", "0"));
    cfg.aperture_cm = std::stod(a.get("aperture-cm", "20"));
    cfg.site_height_m = std::stod(a.get("height", "0"));
    const auto r = build_light_curve(t, cfg);
    write_light_curve_csv(a.need("out"), t, r);
    size_t ngood = 0;
    for (bool g : r.good) ngood += g;
    std::cout << "aperture scan (radius, n_comps, model_rms, oot_rms):\n";
    for (const auto& s : r.scan) std::printf("  r=%-4g n=%-3zu %.5f %.5f\n", s.radius, s.ncomps, s.model_rms, s.oot_rms);
    std::cout << "chosen radius " << r.radius << " px, comparison stars";
    for (size_t i = 0; i < r.comps.size(); ++i) std::printf(" %zu(w=%.3f)", r.comps[i], r.weights[i]);
    std::printf("\ngood frames %zu/%zu (cooler flags %d, saturation flags %d), Tc predicted %.6f (epoch %ld)\n",
                ngood, t.nf, r.n_flag_cooler, r.n_flag_saturated, r.tc_pred, r.epoch);
    std::printf("noise model: photon median %.4f, scintillation %.4f..%.4f\n", r.sig_photon_median, r.sig_scint_min, r.sig_scint_max);
    return 0;
}

int cmd_fit(const Args& a) {
    FitData d = load_fit_data(a.need("lc"));
    if (!a.get("tc-pred").empty()) d.tc_pred = std::stod(a.get("tc-pred"));
    if (!std::isfinite(d.tc_pred)) throw std::runtime_error("no predicted Tc (lc header or --tc-pred)");
    FitConfig cfg;
    const std::string mode = a.get("mode", "constrained");
    if (mode == "free") cfg.mode = FitConfig::Mode::Free;
    else if (mode != "constrained") throw std::invalid_argument("--mode free|constrained");
    cfg.period = std::stod(a.get("P", "2.21857567"));
    cfg.steps = std::stoul(a.get("steps", "16000"));
    cfg.discard = std::stoul(a.get("discard", "2000"));
    cfg.thin = std::stoul(a.get("thin", "10"));
    cfg.walkers = std::stoul(a.get("walkers", "40"));
    cfg.seed = std::stoull(a.get("seed", "1"));
    cfg.threads = std::stoul(a.get("threads", "0"));
    if (!a.get("ld-json").empty()) read_ld_json(a.get("ld-json"), cfg);   // output of `transitscope ld`
    if (!a.get("ld").empty()) {                            // u1,u2,sigma  or  u1,u2,sigma1,sigma2
        auto v = parse_list(a.get("ld"));
        if (v.size() != 3 && v.size() != 4) throw std::invalid_argument("--ld u1,u2,sigma[,sigma2]");
        cfg.u1 = {v[0], v[2]};
        cfg.u2 = {v[1], v.size() == 4 ? v[3] : v[2]};
    }
    if (!a.get("a-prior").empty()) { auto v = parse_list(a.get("a-prior")); cfg.a = {v.at(0), v.at(1)}; cfg.a0 = v[0]; }
    if (!a.get("inc-prior").empty()) { auto v = parse_list(a.get("inc-prior")); cfg.inc = {v.at(0), v.at(1)}; cfg.inc0 = v[0]; }
    std::printf("limb-darkening priors: u1 = %.3f +- %.3f, u2 = %.3f +- %.3f\n", cfg.u1.mu, cfg.u1.sigma, cfg.u2.mu, cfg.u2.sigma);
    const auto r = fit_transit(d, cfg);
    write_fit_json(a.need("out"), r, cfg);
    if (!a.get("model-csv").empty()) write_fit_model_csv(a.get("model-csv"), d.t, d.f, r, cfg);
    std::printf("%s fit, %zu points: Tc = %.6f +%.6f/-%.6f BJD_TDB  (O-C %+.2f min, +/- %.2f min incl. beta=%.2f)\n",
                mode.c_str(), d.t.size(), r.q[0].med, r.q[0].hi, r.q[0].lo, r.o_c_min, r.tc_err_inflated_min, r.beta);
    std::printf("Rp/Rs = %.4f +%.4f/-%.4f  a/Rs = %.3f  inc = %.2f  b = %.3f  depth = %.3f %%  rms = %.2f mmag\n",
                r.q[1].med, r.q[1].hi, r.q[1].lo, r.q[2].med, r.q[3].med, r.b, r.depth_pct, r.rms_resid * 1000);
    std::printf("acceptance %.3f, max autocorrelation time %.0f steps (%.0f tau kept)\n", r.acceptance, r.autocorr_max,
                (cfg.steps - cfg.discard) / r.autocorr_max);
    return 0;
}


int cmd_ld(const Args& a) {
    ClaretTable tab;
    tab.load(find_data_file(a.get("table", "claret2011_atlas_quad"), "ld"));
    const Curve filter = load_curve_csv(find_data_file(a.get("filter", "astrodon_xop_bb"), "passbands"));
    const Curve qe = load_curve_csv(find_data_file(a.get("qe", "sony_icx285al_qe"), "passbands"));
    const StarParams star{std::stod(a.need("teff")), std::stod(a.need("logg")), std::stod(a.need("feh"))};
    LdConfig cfg;
    cfg.star_err = {std::stod(a.get("teff-err", "0")), std::stod(a.get("logg-err", "0")), std::stod(a.get("feh-err", "0"))};
    cfg.filter_shift_nm = std::stod(a.get("filter-shift", "10"));
    cfg.qe_tilt = std::stod(a.get("qe-tilt", "0.10"));
    cfg.model_floor = std::stod(a.get("model-floor", "0.05"));
    const auto r = passband_limb_darkening(tab, star, filter, qe, cfg);
    if (!a.get("out").empty()) write_ld_json(a.get("out"), r, star, cfg);

    std::printf("passband: photon-weighted mean %.0f nm, 90%% of photons in %.0f-%.0f nm\n", r.lambda_eff_nm,
                r.lambda_lo_nm, r.lambda_hi_nm);
    std::printf("%-8s %7s %7s %7s %7s %7s  %s\n", "band", "nm", "a(L)", "b(L)", "a(F)", "b(F)", "leave-one-out da, db");
    for (const auto& b : r.per_band) {
        std::printf("%-8s %7.0f %7.4f %7.4f %7.4f %7.4f", b.name.c_str(), b.lambda_nm, b.aL, b.bL, b.aF, b.bF);
        if (!std::isnan(b.loo_a)) std::printf("  %+.4f %+.4f%s", b.loo_a, b.loo_b, b.in_passband ? "  (in passband)" : "");
        std::printf("\n");
    }
    std::printf("\n        L       F       adopted   sigma: method  interp  passband stellar floor  total\n");
    std::printf("u1  %7.4f %7.4f %9.4f        %7.4f %7.4f %8.4f %7.4f %5.3f %6.4f\n", r.u1_L, r.u1_F, r.u1, r.s1.method,
                r.s1.interp, r.s1.passband, r.s1.stellar, r.s1.floor, r.s1.total);
    std::printf("u2  %7.4f %7.4f %9.4f        %7.4f %7.4f %8.4f %7.4f %5.3f %6.4f\n", r.u2_L, r.u2_F, r.u2, r.s2.method,
                r.s2.interp, r.s2.passband, r.s2.stellar, r.s2.floor, r.s2.total);
    std::printf("\nfit priors: --ld %.3f,%.3f,%.3f,%.3f\n", r.u1, r.u2, r.s1.total, r.s2.total);
    return 0;
}

int cmd_project(const Args& a) {
    if (a.pos.size() != 2) throw std::invalid_argument("usage: transitscope project new|check FILE");
    const std::string& what = a.pos[0];
    const std::string& file = a.pos[1];
    if (what == "new") {
        if (fs::exists(file)) throw std::runtime_error(file + " exists; not overwriting it");
        std::ofstream(file) << project_template();
        std::cerr << "wrote " << file << " (an example: edit the paths, target and ephemeris)\n";
        return 0;
    }
    if (what != "check") throw std::invalid_argument("usage: transitscope project new|check FILE");
    const Project p = load_project(file);
    const auto& f = p.frames;
    auto count = [](const std::vector<std::string>& v) { return std::to_string(v.size()); };
    const std::string lights = p.resolve(f.lights_dir);
    std::cout << "project  " << (p.name.empty() ? "(unnamed)" : p.name) << "\n"
              << "target   " << p.target_name << "  ra " << p.ra << "  dec " << p.dec << "\n"
              << "output   " << p.resolve(p.output_dir) << "\n";
    if (fs::is_directory(lights))
        std::cout << "lights   " << count(find_light_frames(lights, f.lights_prefix)) << " frames in " << lights << "\n";
    else
        std::cout << "lights   MISSING folder " << lights << "\n";
    auto cal = [&](const char* label, const char* word, const std::string& dir, const std::string& master) {
        if (!master.empty()) {
            const std::string m = p.resolve(master);
            std::cout << label << "master " << m << (fs::exists(m) ? "" : "  (MISSING)") << "\n";
        } else if (!dir.empty()) {
            const std::string d = p.resolve(dir);
            if (fs::is_directory(d)) std::cout << label << count(find_calibration_frames(d, word, lights)) << " frames in " << d << "\n";
            else std::cout << label << "MISSING folder " << d << "\n";
        } else {
            std::cout << label << "none\n";
        }
    };
    cal("bias     ", "bias", f.bias_dir, f.master_bias);
    cal("darks    ", "dark", f.darks_dir, f.master_dark);
    cal("flats    ", "flat", f.flats_dir, f.master_flat);
    if (!p.stars_file.empty()) std::cout << "stars    " << p.resolve(p.stars_file) << "\n";
    else std::cout << "stars    " << (p.reduce.stars.empty() ? std::string("detected on the reference frame")
                                                            : std::to_string(p.reduce.stars.size()) + " in the file") << "\n";
    std::cout << "ld       " << (p.ld.mode == LdSettings::Mode::Passband ? "passband (" + p.ld.filter + " x " + p.ld.qe + ")" : std::string("fixed")) << "\n";
    std::cout << "fit      " << p.fit_modes.size() << " mode(s), " << p.fit.steps << " steps\n";
    const auto bad = check_project_inputs(p);
    if (bad.empty()) { std::cout << "\nready to run\n"; return 0; }
    std::cout << "\ncannot run yet:\n";
    for (const auto& b : bad) std::cout << "  " << b << "\n";
    return 1;
}

int cmd_run(const Args& a) {
    if (a.pos.size() != 1) throw std::invalid_argument("usage: transitscope run project.toml [--from STAGE] [--to STAGE] [--only STAGE]");
    const Project p = load_project(a.pos[0]);
    RunOptions opt;
    opt.log = &std::cerr;
    if (!a.get("from").empty()) opt.from = parse_stage(a.get("from"));
    if (!a.get("to").empty()) opt.to = parse_stage(a.get("to"));
    if (!a.get("only").empty()) opt.from = opt.to = parse_stage(a.get("only"));
    if (opt.from > opt.to) throw std::invalid_argument("--from is after --to");
    run_project(p, opt);
    return 0;
}

int cmd_export(const Args& a) {
    if (a.pos.size() != 1) throw std::invalid_argument("usage: transitscope export project.toml --to exoclock|etd [--out FILE]");
    const Project p = load_project(a.pos[0]);
    const ExportTarget t = parse_export_target(a.need("to"));
    const std::string out = a.get("out").empty() ? default_export_path(p, t) : a.get("out");
    const ExportResult r = export_light_curve(p, t, out);
    std::cout << export_target_name(t) << ": " << r.points << " points -> " << r.path << "\n\n" << r.instructions << "\n";
    return 0;
}

int cmd_verify(const Args& a) {
    if (a.pos.size() != 1) throw std::invalid_argument("usage: transitscope verify project.toml [--frames size]");
    const Project p = load_project(a.pos[0]);
    const bool hash_frames = a.get("frames", "hash") != "size";
    const VerifyReport r = verify_provenance(p, hash_frames);
    for (const auto& l : r.lines)
        std::cout << (l.kind == VerifyReport::Ok ? "ok       " : l.kind == VerifyReport::Info ? "--       " : "PROBLEM  ") << l.text << "\n";
    std::cout << (r.ok() ? "\nall recorded results match their inputs\n" : "\nsome results do not match the recorded inputs\n");
    return r.ok() ? 0 : 1;
}

int cmd_diagnose(const Args& a) {
    // project form: the files of a run; file form: --lc lc.csv --model model.csv [--phot phot.csv --target N]
    std::string lc, model, phot, out;
    size_t target = std::stoul(a.get("target", "0"));
    const std::string mode = a.get("mode", "constrained");
    if (a.pos.size() == 1) {
        const Project p = load_project(a.pos[0]);
        lc = p.output("lc.csv");
        model = p.output("model_" + mode + ".csv");
        phot = p.output("phot.csv");
        out = p.output("diagnostics_" + mode + ".json");
        target = p.lightcurve.target;
        if (!std::filesystem::exists(phot)) phot.clear();
    } else if (a.pos.empty()) {
        lc = a.need("lc");
        model = a.need("model");
        phot = a.get("phot");
    } else {
        throw std::invalid_argument("usage: transitscope diagnose project.toml [--mode constrained|free]\n"
                                    "       transitscope diagnose --lc lc.csv --model model.csv [--phot phot.csv --target N] [--out F]");
    }
    if (!a.get("out").empty()) out = a.get("out");
    const FitData d = load_fit_data(lc);
    const std::vector<double> res = load_residuals(model);
    if (res.size() != d.t.size())
        throw std::runtime_error(model + " has " + std::to_string(res.size()) + " points, " + lc + " " +
                                 std::to_string(d.t.size()) + ": the fit is older than the light curve; rerun the fit");
    std::optional<PhotTable> t;
    if (!phot.empty()) t = load_phot_csv(phot);
    const Diagnostics dg = diagnose_fit(t ? &*t : nullptr, lc, d, res, target);
    std::cout << diagnostics_summary(dg);
    std::cout << "\nbinned rms (ppt):   n  bins   measured  white-noise  beta\n";
    for (const auto& b : dg.bins) {
        char buf[128];
        std::snprintf(buf, sizeof buf, "               %5zu %5zu %10.3f %12.3f %5.2f\n", b.n, b.bins, b.rms * 1000, b.white * 1000, b.beta);
        std::cout << buf;
    }
    if (!out.empty()) {
        write_diagnostics_json(out, dg);
        std::cout << "-> " << out << "\n";
    }
    return 0;
}

void usage() {
    std::cerr << "usage: transitscope run project.toml [--from STAGE] [--to STAGE] [--only STAGE]\n"
                 "                  stages: masters reduce lightcurve ld fit\n"
                 "       transitscope verify project.toml [--frames size]   (results vs their recorded inputs)\n"
                 "       transitscope diagnose project.toml [--mode constrained|free]   (red noise, systematics)\n"
                 "       transitscope export project.toml --to exoclock|etd [--out FILE]\n"
                 "       transitscope project new FILE      (write a commented example project)\n"
                 "       transitscope project check FILE    (validate it and show what it finds)\n\n"
                 "       transitscope master --type bias|dark|flat [--bias B] --out M.fits FILES...\n"
                 "       transitscope reduce --lights DIR --dark D --flat F --ra RA --dec DEC --lat LAT --lon LON "
                 "--out phot.csv [options]\n"
                 "       transitscope lightcurve --phot phot.csv --T0 BJD --P DAYS --T14 DAYS [--target 0] [--exclude-stars I,J]\n"
                 "                  [--radii R,...] [--comps I,J,...] --gain E/ADU --read-noise E --aperture-cm D --height M --out lc.csv\n"
                 "       transitscope ld --teff K --logg CGS --feh DEX [--teff-err --logg-err --feh-err] [--filter astrodon_xop_bb]\n"
                 "                  [--qe sony_icx285al_qe] [--filter-shift 10] [--qe-tilt 0.10] [--model-floor 0.05] [--out ld.json]\n"
                 "       transitscope fit --lc lc.csv [--mode constrained|free] [--ld u1,u2,s1[,s2] | --ld-json ld.json]\n"
                 "                  [--steps 16000] [--seed 1] --out fit.json [--model-csv m.csv]\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    try {
        const std::string cmd = argv[1];
        const Args a = parse(argc, argv, 2);
        if (cmd == "master") return cmd_master(a);
        if (cmd == "reduce") return cmd_reduce(a);
        if (cmd == "lightcurve") return cmd_lightcurve(a);
        if (cmd == "fit") return cmd_fit(a);
        if (cmd == "ld") return cmd_ld(a);
        if (cmd == "run") return cmd_run(a);
        if (cmd == "export") return cmd_export(a);
        if (cmd == "verify") return cmd_verify(a);
        if (cmd == "diagnose") return cmd_diagnose(a);
        if (cmd == "project") return cmd_project(a);
        usage();
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "transitscope: " << e.what() << "\n";
        return 1;
    }
}
