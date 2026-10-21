// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/timing.h"

#include <erfa.h>
#include <erfam.h>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace transit {
namespace {

void check(int status, const char* what) {
    if (status < 0) throw std::runtime_error(std::string("ERFA error in ") + what);
}

double parse_sexa(const std::string& s, bool hours) {
    // accepts "20h00m43.713s", "20:00:43.713", "20 00 43.713", "+22d42m39.07s" or plain decimal degrees
    std::string t;
    for (char c : s) t += (c == 'h' || c == 'd' || c == 'm' || c == 's' || c == ':' || c == '\'' || c == '"') ? ' ' : c;
    double a = 0, b = 0, c = 0;
    const int n = std::sscanf(t.c_str(), "%lf %lf %lf", &a, &b, &c);
    if (n <= 0) throw std::invalid_argument("cannot parse angle: " + s);
    if (n == 1) return a;  // decimal degrees
    const bool neg = s.find('-') != std::string::npos;
    const double v = std::fabs(a) + b / 60.0 + c / 3600.0;
    return (neg ? -v : v) * (hours ? 15.0 : 1.0);
}

}  // namespace

double parse_ra(const std::string& s) { return parse_sexa(s, true); }
double parse_dec(const std::string& s) { return parse_sexa(s, false); }

FrameTime frame_time(const std::string& date_obs, double exptime_s, const SkyPos& target, const Site& site) {
    int Y, M, D, h = 0, mi = 0;
    double sec = 0;
    if (std::sscanf(date_obs.c_str(), "%d-%d-%dT%d:%d:%lf", &Y, &M, &D, &h, &mi, &sec) < 3)
        throw std::invalid_argument("bad DATE-OBS: " + date_obs);

    FrameTime ft;
    double u1, u2;
    check(eraDtf2d("UTC", Y, M, D, h, mi, sec, &u1, &u2), "eraDtf2d");
    ft.jd_utc_start = u1 + u2;
    u2 += 0.5 * exptime_s / ERFA_DAYSEC;                 // mid-exposure (UTC)

    double a1, a2, t1, t2;
    check(eraUtctai(u1, u2, &a1, &a2), "eraUtctai");
    check(eraTaitt(a1, a2, &t1, &t2), "eraTaitt");
    // TDB-TT (periodic, < 2 ms); UT1 approximated by UTC, geocentric observer term ignored
    const double jd = u1 + u2;
    const double ut_frac = (jd + 0.5) - std::floor(jd + 0.5);   // fraction of the UT day
    const double dtdb = eraDtdb(t1, t2, ut_frac, site.lon_deg * ERFA_DD2R, 0, 0);
    const double b1 = t1, b2 = t2 + dtdb / ERFA_DAYSEC;

    double pvh[2][3], pvb[2][3];
    check(eraEpv00(b1, b2, pvh, pvb), "eraEpv00");      // Earth, barycentric, AU (BCRS ~ ICRS axes)
    double n[3];
    eraS2c(target.ra_deg * ERFA_DD2R, target.dec_deg * ERFA_DD2R, n);
    const double lt_s = eraPdp(pvb[0], n) * ERFA_AULT;   // Roemer delay
    ft.bjd_tdb_mid = b1 + b2 + lt_s / ERFA_DAYSEC;

    // observed alt/az without refraction (phpa = 0), for airmass
    double aob, zob, hob, dob, rob, eo;
    check(eraAtco13(target.ra_deg * ERFA_DD2R, target.dec_deg * ERFA_DD2R, 0, 0, 0, 0, u1, u2, 0.0,
                    site.lon_deg * ERFA_DD2R, site.lat_deg * ERFA_DD2R, site.height_m, 0, 0, 0, 0, 0, 0.55,
                    &aob, &zob, &hob, &dob, &rob, &eo),
          "eraAtco13");
    ft.alt_deg = 90.0 - zob * ERFA_DR2D;
    ft.airmass = 1.0 / std::sin(ft.alt_deg * ERFA_DD2R);
    return ft;
}

}  // namespace transit
