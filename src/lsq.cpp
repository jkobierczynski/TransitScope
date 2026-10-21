// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/lsq.h"

#include <cmath>
#include <utility>

namespace transit {

bool solve_dense(std::vector<double> A, std::vector<double> b, std::vector<double>& x) {
    const size_t n = b.size();
    for (size_t c = 0; c < n; ++c) {
        size_t piv = c;
        for (size_t r = c + 1; r < n; ++r)
            if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) piv = r;
        if (std::fabs(A[piv * n + c]) < 1e-300) return false;
        if (piv != c) {
            for (size_t k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]);
            std::swap(b[c], b[piv]);
        }
        for (size_t r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] / A[c * n + c];
            for (size_t k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
            b[r] -= f * b[c];
        }
    }
    x.assign(n, 0.0);
    for (size_t i = n; i-- > 0;) {
        double s = b[i];
        for (size_t k = i + 1; k < n; ++k) s -= A[i * n + k] * x[k];
        x[i] = s / A[i * n + i];
    }
    return true;
}

std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y, const std::vector<bool>& mask,
                            int deg) {
    const size_t m = static_cast<size_t>(deg) + 1;
    std::vector<double> A(m * m, 0.0), b(m, 0.0), c;
    for (size_t i = 0; i < x.size(); ++i) {
        if (!mask[i]) continue;
        std::vector<double> pw(2 * m, 1.0);
        for (size_t k = 1; k < 2 * m; ++k) pw[k] = pw[k - 1] * x[i];
        for (size_t r = 0; r < m; ++r) {
            b[r] += pw[r] * y[i];
            for (size_t k = 0; k < m; ++k) A[r * m + k] += pw[r + k];
        }
    }
    if (!solve_dense(A, b, c)) c.assign(m, 0.0);
    return c;
}

double polyval(const std::vector<double>& c, double x) {
    double v = 0;
    for (size_t k = c.size(); k-- > 0;) v = v * x + c[k];
    return v;
}

std::vector<double> levenberg_marquardt(const ResidualFn& f, std::vector<double> p, const std::vector<double>& scale,
                                        int max_iter, double* final_cost) {
    const size_t np = p.size();
    std::vector<double> r, rt;
    f(p, r);
    auto cost = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return s; };
    double c0 = cost(r), lambda = 1e-3;
    std::vector<double> J;
    for (int it = 0; it < max_iter; ++it) {
        const size_t nr = r.size();
        J.assign(nr * np, 0.0);
        for (size_t j = 0; j < np; ++j) {
            auto q = p;
            const double h = scale[j] * 1e-4;
            q[j] += h;
            f(q, rt);
            for (size_t i = 0; i < nr; ++i) J[i * np + j] = (rt[i] - r[i]) / h;
        }
        std::vector<double> JTJ(np * np, 0.0), JTr(np, 0.0);
        for (size_t i = 0; i < nr; ++i)
            for (size_t a = 0; a < np; ++a) {
                JTr[a] += J[i * np + a] * r[i];
                for (size_t b = 0; b < np; ++b) JTJ[a * np + b] += J[i * np + a] * J[i * np + b];
            }
        bool improved = false;
        for (int tries = 0; tries < 20; ++tries) {
            auto A = JTJ;
            for (size_t a = 0; a < np; ++a) A[a * np + a] *= (1.0 + lambda);
            std::vector<double> negJTr(np), dp;
            for (size_t a = 0; a < np; ++a) negJTr[a] = -JTr[a];
            if (!solve_dense(A, negJTr, dp)) { lambda *= 10; continue; }
            auto q = p;
            for (size_t a = 0; a < np; ++a) q[a] += dp[a];
            f(q, rt);
            const double c1 = cost(rt);
            if (std::isfinite(c1) && c1 < c0) {
                const double rel = (c0 - c1) / c0;
                p = q; r = rt; c0 = c1; lambda = std::max(lambda / 10, 1e-12);
                improved = true;
                if (rel < 1e-12) it = max_iter;  // converged
                break;
            }
            lambda *= 10;
        }
        if (!improved) break;
    }
    if (final_cost) *final_cost = c0;
    return p;
}

}  // namespace transit
