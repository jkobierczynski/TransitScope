// SPDX-License-Identifier: GPL-3.0-or-later
// Reads lines "t0 per rp a inc u1 u2 t" and prints the model flux (for comparison with batman).
#include <cstdio>
#include "transit/transit_model.h"
int main() {
    transit::TransitParams p; double t;
    while (std::scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &p.t0, &p.per, &p.rp, &p.a, &p.inc, &p.u1, &p.u2, &t) == 8)
        std::printf("%.17g\n", transit::occult_quad(transit::sky_separation(t, p), p.rp, p.u1, p.u2));
}
