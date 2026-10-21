// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/mcmc.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>

namespace transit {
namespace {

void fft(std::vector<std::complex<double>>& a, bool inverse) {
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2 * M_PI / len * (inverse ? 1 : -1);
        const std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (size_t j = 0; j < len / 2; ++j) {
                auto u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse) for (auto& x : a) x /= static_cast<double>(n);
}

// normalised autocorrelation function of a 1-D series (emcee autocorr.function_1d)
std::vector<double> acf(const std::vector<double>& x) {
    const size_t n = x.size();
    size_t m = 1;
    while (m < 2 * n) m <<= 1;
    const double mean = std::accumulate(x.begin(), x.end(), 0.0) / n;
    std::vector<std::complex<double>> a(m, 0.0);
    for (size_t i = 0; i < n; ++i) a[i] = x[i] - mean;
    fft(a, false);
    for (auto& v : a) v = v * std::conj(v);
    fft(a, true);
    std::vector<double> r(n);
    const double r0 = a[0].real();
    for (size_t i = 0; i < n; ++i) r[i] = r0 > 0 ? a[i].real() / r0 : 0.0;
    return r;
}

}  // namespace

EnsembleSampler::EnsembleSampler(size_t nwalkers, size_t ndim, LogProb lnprob, uint64_t seed, double a, unsigned threads)
    : nw_(nwalkers), nd_(ndim), lnprob_(std::move(lnprob)), seed_(seed), a_(a),
      threads_(threads ? threads : std::max(1u, std::thread::hardware_concurrency())) {
    if (nw_ < 2 * nd_ || nw_ % 2) throw std::invalid_argument("need an even number of walkers >= 2 x ndim");
}

void EnsembleSampler::run(const std::vector<std::vector<double>>& p0, size_t nsteps, const ProgressFn& progress) {
    if (p0.size() != nw_) throw std::invalid_argument("p0 must have one row per walker");
    std::mt19937_64 rng(seed_);
    std::vector<std::vector<double>> pos = p0;
    std::vector<double> lp(nw_);
    for (size_t w = 0; w < nw_; ++w) {
        lp[w] = lnprob_(pos[w]);
        if (!std::isfinite(lp[w])) throw std::runtime_error("initial walker outside the prior");
    }
    chain_.assign(nsteps * nw_ * nd_, 0.0);
    lnp_.assign(nsteps * nw_, 0.0);
    steps_ = nsteps;
    std::vector<size_t> idx(nw_);
    std::iota(idx.begin(), idx.end(), 0);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    for (size_t s = 0; s < nsteps; ++s) {
        std::shuffle(idx.begin(), idx.end(), rng);  // random red/blue split each step
        for (int half = 0; half < 2; ++half) {
            const size_t h = nw_ / 2;
            std::vector<size_t> act(idx.begin() + half * h, idx.begin() + (half + 1) * h);
            std::vector<size_t> comp(idx.begin() + (1 - half) * h, idx.begin() + (2 - half) * h);
            // draw all random numbers serially so results do not depend on the thread count
            std::vector<double> zz(h), uu(h);
            std::vector<size_t> pick(h);
            for (size_t k = 0; k < h; ++k) {
                zz[k] = std::pow((a_ - 1.0) * U(rng) + 1.0, 2) / a_;
                pick[k] = comp[static_cast<size_t>(U(rng) * h) % h];
                uu[k] = U(rng);
            }
            std::vector<std::vector<double>> prop(h, std::vector<double>(nd_));
            std::vector<double> plp(h);
            for (size_t k = 0; k < h; ++k)
                for (size_t d = 0; d < nd_; ++d)
                    prop[k][d] = pos[pick[k]][d] + zz[k] * (pos[act[k]][d] - pos[pick[k]][d]);
            auto work = [&](size_t from, size_t to) { for (size_t k = from; k < to; ++k) plp[k] = lnprob_(prop[k]); };
            std::vector<std::thread> th;
            const size_t per = (h + threads_ - 1) / threads_;
            for (size_t b = 0; b < h; b += per) th.emplace_back(work, b, std::min(h, b + per));
            for (auto& t : th) t.join();
            for (size_t k = 0; k < h; ++k) {
                ++proposed_;
                const double lnq = (nd_ - 1.0) * std::log(zz[k]) + plp[k] - lp[act[k]];
                if (std::isfinite(plp[k]) && std::log(uu[k]) < lnq) {
                    pos[act[k]] = prop[k];
                    lp[act[k]] = plp[k];
                    ++accepted_;
                }
            }
        }
        for (size_t w = 0; w < nw_; ++w) {
            for (size_t d = 0; d < nd_; ++d) chain_[(s * nw_ + w) * nd_ + d] = pos[w][d];
            lnp_[s * nw_ + w] = lp[w];
        }
        if (progress && !progress(s + 1, nsteps)) {
            steps_ = s + 1;  // the steps that were run stay readable
            throw Cancelled();
        }
    }
}

std::vector<std::vector<double>> EnsembleSampler::flat_chain(size_t discard, size_t thin) const {
    std::vector<std::vector<double>> out;
    for (size_t s = discard + thin - 1; s < steps_; s += thin)  // emcee: chain[discard + thin - 1 :: thin]
        for (size_t w = 0; w < nw_; ++w) {
            std::vector<double> v(nd_);
            for (size_t d = 0; d < nd_; ++d) v[d] = at(s, w, d);
            out.push_back(std::move(v));
        }
    return out;
}

std::vector<double> EnsembleSampler::flat_lnprob(size_t discard, size_t thin) const {
    std::vector<double> out;
    for (size_t s = discard + thin - 1; s < steps_; s += thin)
        for (size_t w = 0; w < nw_; ++w) out.push_back(lnp_[s * nw_ + w]);
    return out;
}

double EnsembleSampler::acceptance_fraction() const {
    return proposed_ ? static_cast<double>(accepted_) / proposed_ : 0.0;
}

std::vector<double> EnsembleSampler::autocorr_time(size_t discard) const {
    std::vector<double> tau(nd_, NAN);
    if (steps_ <= discard + 10) return tau;
    const size_t n = steps_ - discard;
    for (size_t d = 0; d < nd_; ++d) {
        std::vector<double> f(n, 0.0);
        for (size_t w = 0; w < nw_; ++w) {
            std::vector<double> x(n);
            for (size_t s = 0; s < n; ++s) x[s] = at(discard + s, w, d);
            const auto r = acf(x);
            for (size_t s = 0; s < n; ++s) f[s] += r[s] / nw_;
        }
        double cum = 0;
        for (size_t m = 0; m < n; ++m) {  // auto window: smallest M with M >= c * tau(M), c = 5
            cum += f[m];
            const double t = 2.0 * cum - 1.0;
            if (m >= 5.0 * t) { tau[d] = t; break; }
            tau[d] = t;
        }
    }
    return tau;
}

}  // namespace transit
