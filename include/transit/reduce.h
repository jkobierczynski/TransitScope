// SPDX-License-Identifier: GPL-3.0-or-later
// Reduction of a light-frame sequence: read, calibrate, track the target, measure every star.
//
// Pipeline (each frame is read and calibrated once):
//   worker threads  read + calibrate frames ahead of time (any order)
//   calling thread  tracks the target frame by frame, outward from the reference frame
//                   (the only sequential step: one centroid per frame)
//   worker threads  measure every star once the frame's offset is known
// Results are stored per frame, so the output is identical for any thread count and backend.
#pragma once
#include <string>
#include <utility>
#include <vector>

#include "transit/calib.h"
#include "transit/phot_table.h"
#include "transit/progress.h"
#include "transit/photometry.h"
#include "transit/timing.h"

namespace transit {

using StarList = std::vector<std::pair<double, double>>;   // (x, y) on the reference frame, px

struct ReduceConfig {
    std::vector<std::string> lights;      // full paths, in frame order (see find_light_frames)
    std::string backend = "cfitsio";
    ApertureConfig aperture;
    SkyPos target{};                      // target coordinates, for BJD_TDB and airmass
    Site site{};
    StarList stars;                       // empty: detected on the reference frame
    size_t ref_index = 0;                 // reference frame (star list, start of tracking)
    size_t target_star = 0;               // index into stars of the tracked star
    double detect_margin = 70;            // detected stars closer than this to an edge are dropped
    unsigned threads = 0;                 // 0 = all cores
};

struct FrameResult {
    std::string file;                     // file name without directory
    std::string date_obs;
    double exptime = 0, ccd_temp = 0;     // NaN when the header has none
    double bjd_tdb = 0, airmass = 0;      // mid-exposure
    double dx = 0, dy = 0;                // target offset from the reference frame, px
    std::vector<StarMeasure> stars;       // one per star; valid = false outside the edge margin
};

struct ReduceResult {
    std::string backend;
    std::vector<double> radii;
    StarList stars;                       // the star list used (given or detected)
    size_t ref_index = 0, target_star = 0;
    std::vector<FrameResult> frames;      // in frame order
};

/// Progress: frames tracked and queued for measuring, of all frames (see progress.h).
using ReduceProgress = ProgressFn;

/// Light frames in `dir` whose name starts with `prefix` and ends in .fit/.fits/.fts, sorted.
std::vector<std::string> find_light_frames(const std::string& dir, const std::string& prefix = "");

/// "x,y" per line; '#' lines skipped.
StarList read_star_list(const std::string& path);

/// Master dark and flat from files; with a master bias, also the hot-pixel mask
/// (pixels with dark - bias > hot_adu).
Calibration load_calibration(const std::string& backend, const std::string& dark, const std::string& flat,
                             const std::string& bias = {}, float hot_adu = 50);

ReduceResult reduce_sequence(const ReduceConfig& cfg, const Calibration& cal, const ReduceProgress& progress = {});

/// One row per frame x star (the format load_phot_csv reads).
void write_phot_csv(const std::string& path, const ReduceResult& r);

/// The same table load_phot_csv gives for the written CSV, without the file (values are not
/// rounded to 16 digits, so they can differ from a CSV round trip in the last bit).
PhotTable to_phot_table(const ReduceResult& r);

}  // namespace transit
