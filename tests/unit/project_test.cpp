// SPDX-License-Identifier: GPL-3.0-or-later
// Project file: the template parses, save -> load gives the same project, paths resolve,
// problems are all reported at once, fixed limb darkening reaches the fit priors.
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "transit/project.h"

using namespace transit;

static int failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}
static std::string error_of(const std::string& text) {
    try { parse_project(text, "/base"); } catch (const std::exception& e) { return e.what(); }
    return "";
}
static bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

int main() {
    const Project t = parse_project(project_template(), "/data/obs");
    check(t.name == "HD 189733 b 2014-09-27" && t.reduce.ref_index == 99 && t.lightcurve.exclude_stars.size() == 3 &&
              t.fit_modes.size() == 2 && t.ld.mode == LdSettings::Mode::Passband,
          "the template parses");
    check(std::abs(t.reduce.target.ra_deg - 300.18213750) < 1e-6 && std::abs(t.reduce.target.dec_deg - 22.71085278) < 1e-6,
          "ra/dec strings are parsed");
    check(t.fit.period == t.lightcurve.eph.P && t.lightcurve.score_model.per == t.lightcurve.eph.P &&
              t.lightcurve.site_height_m == t.reduce.site.height_m && t.lightcurve.target == t.reduce.target_star,
          "shared values are copied to every stage");
    check(t.fit.a0 == t.fit.a.mu && t.fit.inc0 == t.fit.inc.mu, "fit walkers start at the a/Rs and inclination prior means");

    // round trip
    const std::string once = project_to_toml(t);
    const Project back = parse_project(once, "/data/obs");
    check(project_to_toml(back) == once, "save -> load -> save gives the same text");
    Project inl = t;
    inl.stars_file.clear();
    inl.reduce.stars = {{627.7731739714469, 648.9507160567141}, {1.5, -2.25e-7}};
    const Project inl2 = parse_project(project_to_toml(inl), "/x");
    check(inl2.reduce.stars == inl.reduce.stars, "inline star list survives a round trip bit for bit");
    Project deg = t;
    deg.ra = "300.182137";
    check(parse_project(project_to_toml(deg), "/x").reduce.target.ra_deg == 300.182137, "ra in degrees stays a number");

    // paths
    check(t.resolve("lights") == "/data/obs/lights" && t.resolve("/abs/x") == "/abs/x" &&
              t.resolve("../other/./f.fits") == "/data/other/f.fits" && t.output("lc.csv") == "/data/obs/out/lc.csv",
          "relative paths resolve against the project's folder");
    if (const char* home = std::getenv("HOME")) check(t.resolve("~/astro") == std::string(home) + "/astro", "~/ is the home folder");

    // errors
    const std::string e = error_of(R"(
[target]
ra = "20h00m43.713s"
dekl = 1
[site]
lat = "north"
lon = 5.4
[frames]
lights_dir = "l"
darks_dir = "d"
flats_dir = "f"
[ephemeris]
T0 = 1
P = 2
[fitt]
)");
    check(contains(e, "[target] dec is required") && contains(e, "unknown key [target] dekl (line 4)") &&
              contains(e, "[site] lat must be a number (line 6)") && contains(e, "[ephemeris] T14 is required") &&
              contains(e, "unknown section [fitt]") && contains(e, "give bias_dir or master_bias"),
          "all problems are reported at once, with line numbers");
    check(contains(error_of("[target\n"), "line 1"), "a TOML syntax error names the line");
    const std::string base = "[target]\nra=1\ndec=2\n[site]\nlat=1\nlon=2\n[frames]\nlights_dir='l'\nmaster_dark='d'\nmaster_flat='f'\n"
                             "[ephemeris]\nT0=1\nP=2\nT14=0.1\n";
    check(contains(error_of(base + "[fit]\nsteps = 10.5\n"), "[fit] steps must be a whole number"), "integers are checked");
    check(contains(error_of(base + "[limb_darkening]\nmode='passband'\nteff=5000\nu1=0.4\n"), "are for mode = \"fixed\""),
          "u1 with passband mode is rejected");
    const Project fixed = parse_project(base + "[limb_darkening]\nmode='fixed'\nu1=0.4\nu2=0.2\nu1_sigma=0.05\nu2_sigma=0.06\n", "/b");
    check(fixed.fit.u1.mu == 0.4 && fixed.fit.u2.sigma == 0.06, "fixed limb darkening sets the fit priors");
    check(error_of(base + "[limb_darkening]\nmode='fixed'\n").empty(), "a minimal project with fixed limb darkening parses");

    // save / load
    const std::string path = "project_test.toml";
    save_project(path, t);
    const Project loaded = load_project(path);
    check(project_to_toml(loaded) == once && !std::filesystem::exists(path + ".tmp"), "save_project + load_project");
    check(loaded.base_dir == std::filesystem::absolute(path).parent_path().string(), "load_project: paths relative to the file");
    std::filesystem::remove(path);

    std::printf("%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
