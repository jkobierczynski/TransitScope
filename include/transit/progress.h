// SPDX-License-Identifier: GPL-3.0-or-later
// Progress reporting and cancelling for long library calls (reduce_sequence, fit_transit).
#pragma once
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace transit {

/// Called on the thread that made the library call with (units done, units in total): frames for
/// reduce_sequence, MCMC steps for fit_transit. Return false to cancel; the call then stops,
/// waits for work already queued, and throws Cancelled.
using ProgressFn = std::function<bool(size_t done, size_t total)>;

struct Cancelled : std::runtime_error {
    Cancelled() : std::runtime_error("cancelled") {}
};

}  // namespace transit
