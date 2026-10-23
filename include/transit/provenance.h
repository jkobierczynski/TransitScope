// SPDX-License-Identifier: GPL-3.0-or-later
// Provenance: what produced each result. Every stage of run_project writes
// <output_dir>/provenance/<stage>.json with the program version and build, the start and end
// time, the complete settings, and the size and SHA-256 of every input and output file.
// verify_provenance() checks those records against the files on disk, so a result can be traced
// back to the exact frames and settings, and a stale or edited output is detected.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace transit {

std::string sha256_hex(const void* data, size_t size);
std::string sha256_file(const std::string& path);   // throws if the file cannot be read

struct BuildInfo { std::string version, git, compiler, system; };
BuildInfo build_info();

struct FileDigest {
    std::string role;      // "light", "bias", "master_dark", "phot", ...
    std::string path;      // as resolved when the stage ran
    uint64_t size = 0;
    std::string sha256;
};

struct ProvenanceRecord {
    std::string stage;
    std::string started_utc, finished_utc;   // ISO 8601
    std::string project_file;                // the project file, if it was loaded from one
    BuildInfo build;
    std::string settings_toml;               // the complete settings (project_to_toml)
    std::vector<FileDigest> inputs, outputs;
};

/// Digests of files, computed on up to `threads` threads (0 = all cores); same order as the input.
std::vector<FileDigest> digest_files(const std::vector<std::pair<std::string, std::string>>& role_path,
                                     unsigned threads = 0);

void write_provenance(const std::string& path, const ProvenanceRecord& r);
ProvenanceRecord read_provenance(const std::string& path);

std::string utc_now_iso();

}  // namespace transit
