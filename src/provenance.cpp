// SPDX-License-Identifier: GPL-3.0-or-later
#include "transit/provenance.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

#include "transit/thread_pool.h"

namespace transit {

// ---------------------------------------------------------------------------------------------
// SHA-256 (FIPS 180-4)

namespace {

class Sha256 {
public:
    Sha256() { reset(); }
    void reset() {
        h_ = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
        len_ = 0;
        fill_ = 0;
    }
    void update(const uint8_t* p, size_t n) {
        len_ += n;
        if (fill_) {
            const size_t k = std::min(n, size_t(64) - fill_);
            std::memcpy(buf_.data() + fill_, p, k);
            fill_ += k;
            p += k;
            n -= k;
            if (fill_ == 64) { block(buf_.data()); fill_ = 0; }
        }
        for (; n >= 64; p += 64, n -= 64) block(p);
        if (n) { std::memcpy(buf_.data(), p, n); fill_ = n; }
    }
    std::string hex() {
        const uint64_t bits = len_ * 8;
        const uint8_t pad = 0x80, zero = 0;
        update(&pad, 1);
        while (fill_ != 56) update(&zero, 1);
        uint8_t l[8];
        for (int i = 0; i < 8; ++i) l[i] = uint8_t(bits >> (56 - 8 * i));
        update(l, 8);
        char out[65];
        for (int i = 0; i < 8; ++i) std::snprintf(out + 8 * i, 9, "%08x", h_[i]);
        return std::string(out, 64);
    }

private:
    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
    void block(const uint8_t* p) {
        static constexpr uint32_t K[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) w[i] = uint32_t(p[4 * i]) << 24 | uint32_t(p[4 * i + 1]) << 16 | uint32_t(p[4 * i + 2]) << 8 | p[4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3], e = h_[4], f = h_[5], g = h_[6], h = h_[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i];
            const uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h_[0] += a; h_[1] += b; h_[2] += c; h_[3] += d; h_[4] += e; h_[5] += f; h_[6] += g; h_[7] += h;
    }
    std::array<uint32_t, 8> h_;
    std::array<uint8_t, 64> buf_;
    uint64_t len_ = 0;
    size_t fill_ = 0;
};

std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); o += b; }
                else o += c;
        }
    }
    return o;
}

std::string json_unescape(const std::string& s) {
    std::string o;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] != '\\' || i + 1 >= s.size()) { o += s[i]; continue; }
        const char c = s[++i];
        switch (c) {
            case 'n': o += '\n'; break;
            case 'r': o += '\r'; break;
            case 't': o += '\t'; break;
            case 'u':
                if (i + 4 < s.size()) { o += char(std::stoi(s.substr(i + 1, 4), nullptr, 16)); i += 4; }
                break;
            default: o += c;
        }
    }
    return o;
}

}  // namespace

std::string sha256_hex(const void* data, size_t size) {
    Sha256 s;
    s.update(static_cast<const uint8_t*>(data), size);
    return s.hex();
}

std::string sha256_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path);
    Sha256 s;
    std::vector<char> buf(1 << 20);
    while (in) {
        in.read(buf.data(), std::streamsize(buf.size()));
        const auto n = in.gcount();
        if (n > 0) s.update(reinterpret_cast<const uint8_t*>(buf.data()), size_t(n));
    }
    return s.hex();
}

BuildInfo build_info() {
    BuildInfo b;
#ifdef TRANSITSCOPE_VERSION
    b.version = TRANSITSCOPE_VERSION;
#endif
#ifdef TRANSITSCOPE_GIT
    b.git = TRANSITSCOPE_GIT;
#endif
#if defined(__clang__)
    b.compiler = std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    b.compiler = std::string("GCC ") + __VERSION__;
#elif defined(_MSC_VER)
    b.compiler = "MSVC " + std::to_string(_MSC_VER);
#endif
#if defined(_WIN32)
    b.system = "Windows";
#elif defined(__APPLE__)
    b.system = "macOS";
#elif defined(__linux__)
    b.system = "Linux";
#else
    b.system = "unknown";
#endif
    return b;
}

std::string utc_now_iso() {
    const std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::vector<FileDigest> digest_files(const std::vector<std::pair<std::string, std::string>>& role_path, unsigned threads) {
    std::vector<FileDigest> out(role_path.size());
    ThreadPool pool(threads);
    std::vector<std::future<void>> jobs;
    for (size_t i = 0; i < role_path.size(); ++i)
        jobs.push_back(pool.submit([&, i] {
            FileDigest& d = out[i];
            d.role = role_path[i].first;
            d.path = role_path[i].second;
            d.size = std::filesystem::file_size(d.path);
            d.sha256 = sha256_file(d.path);
        }));
    for (auto& j : jobs) j.get();
    return out;
}

void write_provenance(const std::string& path, const ProvenanceRecord& r) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("cannot write " + path);
    auto files = [&](const char* name, const std::vector<FileDigest>& v) {
        o << "  \"" << name << "\": [\n";
        for (size_t i = 0; i < v.size(); ++i)
            o << "    {\"role\": \"" << json_escape(v[i].role) << "\", \"path\": \"" << json_escape(v[i].path)
              << "\", \"size\": " << v[i].size << ", \"sha256\": \"" << v[i].sha256 << "\"}"
              << (i + 1 < v.size() ? "," : "") << "\n";
        o << "  ]";
    };
    o << "{\n"
      << "  \"stage\": \"" << json_escape(r.stage) << "\",\n"
      << "  \"started_utc\": \"" << r.started_utc << "\",\n"
      << "  \"finished_utc\": \"" << r.finished_utc << "\",\n"
      << "  \"software\": \"TransitScope\",\n"
      << "  \"version\": \"" << json_escape(r.build.version) << "\",\n"
      << "  \"git\": \"" << json_escape(r.build.git) << "\",\n"
      << "  \"compiler\": \"" << json_escape(r.build.compiler) << "\",\n"
      << "  \"system\": \"" << json_escape(r.build.system) << "\",\n"
      << "  \"project_file\": \"" << json_escape(r.project_file) << "\",\n";
    files("inputs", r.inputs);
    o << ",\n";
    files("outputs", r.outputs);
    o << ",\n  \"settings_toml\": \"" << json_escape(r.settings_toml) << "\"\n}\n";
    if (!o) throw std::runtime_error("error writing " + path);
}

ProvenanceRecord read_provenance(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    ProvenanceRecord r;
    static const std::regex file(
        R"re(^\s*\{"role": "((?:[^"\\]|\\.)*)", "path": "((?:[^"\\]|\\.)*)", "size": (\d+), "sha256": "([0-9a-f]{64})"\},?\s*$)re");
    std::vector<FileDigest>* list = nullptr;
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::smatch m;
        if (line.find("\"inputs\": [") != std::string::npos) { list = &r.inputs; continue; }
        if (line.find("\"outputs\": [") != std::string::npos) { list = &r.outputs; continue; }
        if (list && std::regex_match(line, m, file)) {
            list->push_back({json_unescape(m[1]), json_unescape(m[2]), std::stoull(m[3]), m[4]});
            continue;
        }
        // "key": "value" lines, scanned by hand (the settings value is long; std::regex recursion
        // on it could exhaust the stack)
        const size_t q0 = line.find('"');
        const size_t q1 = q0 == std::string::npos ? q0 : line.find('"', q0 + 1);
        if (q1 != std::string::npos && line.compare(q1 + 1, 3, ": \"") == 0) {
            const std::string k = line.substr(q0 + 1, q1 - q0 - 1);
            size_t e = q1 + 4;
            std::string raw;
            for (; e < line.size() && line[e] != '"'; ++e) {
                if (line[e] == '\\' && e + 1 < line.size()) raw += line[e++];
                raw += line[e];
            }
            const std::string v = json_unescape(raw);
            if (k == "stage") r.stage = v;
            else if (k == "started_utc") r.started_utc = v;
            else if (k == "finished_utc") r.finished_utc = v;
            else if (k == "version") r.build.version = v;
            else if (k == "git") r.build.git = v;
            else if (k == "compiler") r.build.compiler = v;
            else if (k == "system") r.build.system = v;
            else if (k == "project_file") r.project_file = v;
            else if (k == "settings_toml") r.settings_toml = v;
        }
    }
    if (r.stage.empty()) throw std::runtime_error(path + ": not a provenance record");
    return r;
}

}  // namespace transit
