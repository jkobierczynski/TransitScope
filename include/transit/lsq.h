// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <functional>
#include <vector>

namespace transit {

/// Solve A x = b (n x n, row-major) by Gaussian elimination with partial pivoting. Returns false if singular.
bool solve_dense(std::vector<double> A, std::vector<double> b, std::vector<double>& x);

/// Least-squares polynomial fit y ~ sum c_k x^k (k = 0..deg) over the points where mask is true.
std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y,
                            const std::vector<bool>& mask, int deg);
double polyval(const std::vector<double>& c, double x);

/// Levenberg-Marquardt on a residual vector function, numerical Jacobian (forward differences
/// with steps scale[i] * 1e-4). Returns the best parameters; `final_cost` = sum of squares.
using ResidualFn = std::function<void(const std::vector<double>& p, std::vector<double>& r)>;
std::vector<double> levenberg_marquardt(const ResidualFn& f, std::vector<double> p0, const std::vector<double>& scale,
                                        int max_iter = 200, double* final_cost = nullptr);

}  // namespace transit
