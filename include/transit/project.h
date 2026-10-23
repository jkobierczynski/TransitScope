// SPDX-License-Identifier: GPL-3.0-or-later
// Project file (TOML): every setting of one observation's reduction, shared by the CLI
// (`transitscope run project.toml`) and the GUI. See project_template() for the format,
// or `transitscope project new FILE`.
//
// Relative paths in the file are relative to the file's own directory; "~/" is the home directory.
// Unknown keys and wrong types are errors (a typo must not silently fall back to a default).
#pragma once
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

#include "transit/ensemble.h"
#include "transit/fit.h"
#include "transit/limb_darkening.h"
#include "transit/reduce.h"

namespace transit {

struct FrameSources {
    std::string backend = "cfitsio";
    std::string lights_dir, lights_prefix;
    std::string flats_dir, bias_dir, darks_dir;          // raw calibration frames (masters built from these)
    std::string master_bias, master_dark, master_flat;   // or ready-made masters (used as given)
    float hot_adu = 50;                                  // hot pixel: dark - bias above this
};

struct LdSettings {
    enum class Mode { Passband, Fixed } mode = Mode::Passband;
    // Passband: computed by `ld` from the tables for this star and passband
    StarParams star{5050, 4.59, -0.03};
    LdConfig cfg;                                        // error budget settings and stellar parameter errors
    std::string table = "claret2011_atlas_quad", filter = "astrodon_xop_bb", qe = "sony_icx285al_qe";
    // Fixed: priors given directly
    Gaussian u1{0.50, 0.10}, u2{0.19, 0.10};
};

struct Project {
    std::string path;                    // the file this was loaded from (empty for a new project)
    std::string base_dir;                // relative paths are resolved against this
    std::string name;
    std::string output_dir = "out";      // all stage outputs go here
    std::string target_name;
    std::string planet;                  // as the databases name it, e.g. "HD189733b" (exports)
    std::string filter_name;             // filter as named in upload forms (exports)
    std::string ra, dec;                 // as written ("20h00m43.713s" or degrees)
    std::string stars_file;              // alternative to reduce.stars
    FrameSources frames;
    ReduceConfig reduce;                 // lights, target and site are filled from the sections above
    EnsembleConfig lightcurve;
    LdSettings ld;
    FitConfig fit;
    std::vector<FitConfig::Mode> fit_modes{FitConfig::Mode::Constrained, FitConfig::Mode::Free};

    /// A path from the file, resolved against base_dir.
    std::string resolve(const std::string& p) const;
    /// Output file <output_dir>/<name>.
    std::string output(const std::string& file) const;
};

/// Parse a project file; throws std::runtime_error listing every problem found.
Project load_project(const std::string& path);
/// The same from text (base_dir: where relative paths point).
Project parse_project(const std::string& toml_text, const std::string& base_dir);
/// TOML text of a project; load_project(save) gives the same project back.
std::string project_to_toml(const Project& p);
void save_project(const std::string& path, const Project& p);
/// A commented example project (the 2014-09-27 HD 189733 b set), the reference for the format.
std::string project_template();

enum class Stage { Masters, Reduce, LightCurve, LimbDarkening, Fit };
const char* stage_name(Stage s);           // "masters", "reduce", "lightcurve", "ld", "fit"
Stage parse_stage(const std::string& s);

/// progress(stage, done, total): frames for masters and reduce, MCMC steps for fit (per mode).
/// Return false to cancel (throws Cancelled).
using StageProgress = std::function<bool(Stage stage, size_t done, size_t total)>;

struct RunOptions {
    Stage from = Stage::Masters, to = Stage::Fit;
    std::ostream* log = nullptr;           // short summary per stage
};

/// Every missing input of a run (folders, frames, star list, master files, data files, and the
/// outputs of earlier stages when a run starts later); empty when the run can start.
/// run_project checks this first, so nothing is computed when something is missing.
std::vector<std::string> check_project_inputs(const Project& p, const RunOptions& opt = {});

/// Runs stages from..to. Each stage writes its output to output_dir and the next stage reads it
/// from there, so a run can start at any stage once the earlier outputs exist:
///   masters     master_bias.fits, master_dark.fits, master_flat.fits (skipped when masters are given)
///   reduce      phot.csv
///   lightcurve  lc.csv
///   ld          ld.json (only with limb_darkening.mode = "passband")
///   fit         fit_<mode>.json, model_<mode>.csv for each mode
void run_project(const Project& p, const RunOptions& opt = {}, const StageProgress& progress = {});

/// <output_dir>/provenance/<stage>.json
std::string provenance_path(const Project& p, Stage s);

struct VerifyReport {
    enum Kind { Ok, Info, Problem };
    struct Line { Kind kind; std::string text; };
    std::vector<Line> lines;
    bool ok() const;   // no problems
};
/// Checks the provenance records against the files: every recorded input and output must still
/// have its recorded size and SHA-256, and each stage must have used the outputs the previous
/// stage recorded. hash_frames = false compares raw frames by size only (much faster).
VerifyReport verify_provenance(const Project& p, bool hash_frames = true);

}  // namespace transit
