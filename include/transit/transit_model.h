// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <vector>

namespace transit {

/// Circular-orbit transit with quadratic limb darkening, I(mu) = 1 - u1(1-mu) - u2(1-mu)^2.
struct TransitParams {
    double t0 = 0;     // mid-transit time (same time scale as t)
    double per = 1;    // period, days
    double rp = 0.1;   // planet/star radius ratio
    double a = 10;     // semi-major axis / stellar radius
    double inc = 90;   // inclination, degrees
    double u1 = 0.4, u2 = 0.2;
};

/// Relative flux for a planet of radius p (stellar radii) at projected separation z.
/// Mandel & Agol (2002), case structure as in Eastman et al. (2013) occultquad; complete
/// elliptic integrals by AGM (K, E) and Bulirsch's cel (Pi) to ~1e-14.
double occult_quad(double z, double p, double u1, double u2);

/// Projected separation in stellar radii; returns a large value when the planet is behind the star.
double sky_separation(double t, const TransitParams& tp);

std::vector<double> light_curve(const TransitParams& tp, const std::vector<double>& t);
void light_curve(const TransitParams& tp, const std::vector<double>& t, std::vector<double>& out);

}  // namespace transit
