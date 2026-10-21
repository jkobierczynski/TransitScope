// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <string>

namespace transit {

struct Site { double lat_deg = 0, lon_deg = 0, height_m = 0; };
struct SkyPos { double ra_deg = 0, dec_deg = 0; };   // ICRS

struct FrameTime {
    double jd_utc_start = 0;   // as recorded (DATE-OBS)
    double bjd_tdb_mid = 0;    // barycentric dynamical time at mid-exposure
    double airmass = 0;        // plane-parallel 1/sin(alt), no refraction
    double alt_deg = 0;
};

/// DATE-OBS (ISO, UTC, exposure start) + exposure -> BJD_TDB at mid-exposure.
/// Earth barycentric position from ERFA eraEpv00 (light-time accuracy ~15 us);
/// the observer's geocentric offset (< 21 ms) is neglected.
FrameTime frame_time(const std::string& date_obs_utc, double exptime_s, const SkyPos& target, const Site& site);

/// Parse "20h00m43.713s" / "20:00:43.713" / decimal degrees.
double parse_ra(const std::string& s);
double parse_dec(const std::string& s);

}  // namespace transit
