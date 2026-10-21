// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/transit_model.h"

#include <algorithm>
#include <cmath>

namespace transit {
namespace {

constexpr double kPi = 3.14159265358979323846;

// Complete elliptic integrals of the first and second kind, modulus k, via the AGM.
void ellke(double k, double& E, double& K) {
    double a = 1.0, b = std::sqrt(std::max(0.0, 1.0 - k * k)), c = k;
    double sum = 0.5 * c * c, pw = 0.5;
    for (int i = 0; i < 60 && std::fabs(c) > 1e-16; ++i) {
        const double an = 0.5 * (a + b);
        c = 0.5 * (a - b);
        b = std::sqrt(a * b);
        a = an;
        pw *= 2.0;
        sum += pw * c * c;
    }
    K = kPi / (2.0 * a);
    E = K * (1.0 - sum);
}

// Complete elliptic integral of the third kind as used by Mandel & Agol: Bulirsch (1965) cel
// with kc = sqrt(1 - k^2), p = n + 1, a = b = 1.
double ellpic(double n, double k) {
    double kc = std::sqrt(std::max(1e-300, 1.0 - k * k));
    double p = std::sqrt(n + 1.0);
    double m0 = 1.0, c = 1.0, d = 1.0 / p, e = kc;
    for (int i = 0; i < 100; ++i) {
        const double f = c;
        c = d / p + c;
        double g = e / p;
        d = 2.0 * (f * g + d);
        p = g + p;
        g = m0;
        m0 = kc + m0;
        if (std::fabs(1.0 - kc / g) <= 1e-14) break;
        kc = 2.0 * std::sqrt(e);
        e = kc * m0;
    }
    return 0.5 * kPi * (c * m0 + d) / (m0 * (m0 + p));
}

}  // namespace

double occult_quad(double z, double p, double u1, double u2) {
    p = std::fabs(p);
    if (p <= 0.0 || z >= 1.0 + p) return 1.0;
    const double tol = 1e-14;
    if (std::fabs(p - z) < tol) z = p;
    if (std::fabs((p - 1.0) - z) < tol) z = p - 1.0;
    if (std::fabs((1.0 - p) - z) < tol) z = 1.0 - p;
    if (z < tol) z = 0.0;

    const double omega = 1.0 - u1 / 3.0 - u2 / 6.0;
    const double x1 = (p - z) * (p - z), x2 = (p + z) * (p + z), x3 = p * p - z * z;
    double lambdad = 0.0, etad = 0.0, lambdae = 0.0;

    // Case 11: source completely occulted
    if (p >= 1.0 && z <= p - 1.0) return 0.0;

    // Cases 2, 7, 8: ingress/egress, uniform-disk part and eta_1
    if (z >= std::fabs(1.0 - p) && z < 1.0 + p) {
        const double kap1 = std::acos(std::clamp((1.0 - p * p + z * z) / (2.0 * z), -1.0, 1.0));
        const double kap0 = std::acos(std::clamp((p * p + z * z - 1.0) / (2.0 * p * z), -1.0, 1.0));
        const double s = std::max(0.0, 4.0 * z * z - std::pow(1.0 + z * z - p * p, 2));
        lambdae = (p * p * kap0 + kap1 - 0.5 * std::sqrt(s)) / kPi;
        etad = 1.0 / (2.0 * kPi) *
               (kap1 + p * p * (p * p + 2.0 * z * z) * kap0 -
                (1.0 + 5.0 * p * p + z * z) / 4.0 * std::sqrt(std::max(0.0, (1.0 - x1) * (x2 - 1.0))));
    }

    if (z == p) {
        // Cases 5, 6, 7: edge of the planet at the centre of the star
        if (p < 0.5) {
            double E, K;
            ellke(2.0 * p, E, K);
            lambdad = 1.0 / 3.0 + 2.0 / (9.0 * kPi) * (4.0 * (2.0 * p * p - 1.0) * E + (1.0 - 4.0 * p * p) * K);
            etad = p * p / 2.0 * (p * p + 2.0 * z * z);
            lambdae = p * p;
        } else if (p > 0.5) {
            double E, K;
            ellke(0.5 / p, E, K);
            lambdad = 1.0 / 3.0 + 16.0 * p / (9.0 * kPi) * (2.0 * p * p - 1.0) * E -
                      (32.0 * p * p * p * p - 20.0 * p * p + 3.0) / (9.0 * kPi * p) * K;
        } else {
            lambdad = 1.0 / 3.0 - 4.0 / (kPi * 9.0);
            etad = 3.0 / 32.0;
        }
    } else if ((z > 0.5 + std::fabs(p - 0.5) && z < 1.0 + p) || (p > 0.5 && z > std::fabs(1.0 - p) && z < p)) {
        // Cases 2, 8: ingress/egress with limb darkening, lambda_1
        const double q = std::sqrt((1.0 - x1) / (x2 - x1));
        double E, K;
        ellke(q, E, K);
        const double n = 1.0 / x1 - 1.0;
        lambdad = 2.0 / (9.0 * kPi * std::sqrt(x2 - x1)) *
                  (((1.0 - x2) * (2.0 * x2 + x1 - 3.0) - 3.0 * x3 * (x2 - 2.0)) * K +
                   (x2 - x1) * (z * z + 7.0 * p * p - 4.0) * E - 3.0 * x3 / x1 * ellpic(n, q));
    }

    if (p < 1.0 && z <= 1.0 - p) {
        // Cases 3, 4, 9, 10: planet completely inside the stellar disk
        etad = p * p / 2.0 * (p * p + 2.0 * z * z);
        lambdae = p * p;
        if (z == 1.0 - p) {
            lambdad = 2.0 / (3.0 * kPi) * std::acos(1.0 - 2.0 * p) -
                      4.0 / (9.0 * kPi) * std::sqrt(p * (1.0 - p)) * (3.0 + 2.0 * p - 8.0 * p * p) -
                      (p > 0.5 ? 2.0 / 3.0 : 0.0);
        } else if (z == 0.0) {
            lambdad = -2.0 / 3.0 * std::pow(1.0 - p * p, 1.5);
        } else if (z != p) {
            const double q = std::sqrt((x2 - x1) / (1.0 - x1));
            const double n = x2 / x1 - 1.0;
            double E, K;
            ellke(q, E, K);
            lambdad = 2.0 / (9.0 * kPi * std::sqrt(1.0 - x1)) *
                      ((1.0 - 5.0 * z * z + p * p + x3 * x3) * K + (1.0 - x1) * (z * z + 7.0 * p * p - 4.0) * E -
                       3.0 * x3 / x1 * ellpic(n, q));
        }
    }

    const double heav = (p > z) ? 1.0 : 0.0;
    return 1.0 - ((1.0 - u1 - 2.0 * u2) * lambdae + (u1 + 2.0 * u2) * (lambdad + 2.0 / 3.0 * heav) + u2 * etad) / omega;
}

double sky_separation(double t, const TransitParams& tp) {
    const double phase = 2.0 * kPi * (t - tp.t0) / tp.per;
    const double ci = std::cos(tp.inc * kPi / 180.0);
    if (std::cos(phase) <= 0) return 1e6;  // planet behind the star (secondary eclipse not modelled)
    const double s = std::sin(phase), c = std::cos(phase);
    return tp.a * std::sqrt(s * s + ci * ci * c * c);
}

void light_curve(const TransitParams& tp, const std::vector<double>& t, std::vector<double>& out) {
    out.resize(t.size());
    for (size_t i = 0; i < t.size(); ++i) out[i] = occult_quad(sky_separation(t[i], tp), tp.rp, tp.u1, tp.u2);
}

std::vector<double> light_curve(const TransitParams& tp, const std::vector<double>& t) {
    std::vector<double> out;
    light_curve(tp, t, out);
    return out;
}

}  // namespace transit
