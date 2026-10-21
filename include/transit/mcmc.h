// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <functional>
#include <vector>

#include "transit/progress.h"

namespace transit {

/// Affine-invariant ensemble sampler (Goodman & Weare 2010), the "stretch move" with the
/// red/blue split used by emcee, so results are directly comparable with emcee runs.
class EnsembleSampler {
public:
    using LogProb = std::function<double(const std::vector<double>&)>;

    EnsembleSampler(size_t nwalkers, size_t ndim, LogProb lnprob, uint64_t seed = 1, double a = 2.0,
                    unsigned threads = 0);

    /// Runs `nsteps` from the given initial positions (nwalkers x ndim). `progress(step, nsteps)`
    /// is called after every step on the calling thread; returning false stops the run and
    /// throws Cancelled. The chain is the same with or without a progress callback.
    void run(const std::vector<std::vector<double>>& p0, size_t nsteps, const ProgressFn& progress = {});

    size_t nsteps() const { return steps_; }
    /// Flattened samples after `discard` steps, every `thin`-th step: [sample][dim].
    std::vector<std::vector<double>> flat_chain(size_t discard, size_t thin) const;
    std::vector<double> flat_lnprob(size_t discard, size_t thin) const;
    double acceptance_fraction() const;
    /// Integrated autocorrelation time per parameter (emcee's estimator, window c = 5),
    /// computed on the steps after `discard`.
    std::vector<double> autocorr_time(size_t discard) const;

private:
    double at(size_t step, size_t w, size_t d) const { return chain_[(step * nw_ + w) * nd_ + d]; }
    size_t nw_, nd_;
    LogProb lnprob_;
    uint64_t seed_;
    double a_;
    unsigned threads_;
    size_t steps_ = 0;
    std::vector<double> chain_, lnp_;  // [step][walker][dim], [step][walker]
    size_t accepted_ = 0, proposed_ = 0;
};

}  // namespace transit
