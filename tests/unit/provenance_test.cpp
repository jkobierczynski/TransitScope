// SPDX-License-Identifier: GPL-3.0-or-later
// SHA-256 against the FIPS 180-4 test vectors; provenance records survive a write/read round trip
// (paths with spaces, quotes, backslashes; multi-line settings).
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "transit/provenance.h"

using namespace transit;

static int failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok) ++failures;
}

int main() {
    check(sha256_hex("", 0) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "SHA-256 of the empty string");
    check(sha256_hex("abc", 3) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "SHA-256 of \"abc\"");
    const std::string two = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    check(sha256_hex(two.data(), two.size()) == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
          "SHA-256 of the 448-bit test message");
    const std::string mil(1000000, 'a');
    check(sha256_hex(mil.data(), mil.size()) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
          "SHA-256 of one million 'a'");
    {
        std::ofstream("prov_test.bin", std::ios::binary) << mil;   // file path, read in 1 MiB blocks
    }
    check(sha256_file("prov_test.bin") == sha256_hex(mil.data(), mil.size()), "sha256_file equals sha256_hex");

    ProvenanceRecord r;
    r.stage = "reduce";
    r.started_utc = "2026-10-24T00:00:00Z";
    r.finished_utc = utc_now_iso();
    r.project_file = "C:\\Users\\PC\\My \"projects\"\\p.toml";
    r.build = build_info();
    r.settings_toml = "[project]\nname = \"HD 189733 b\"\n\ttab and \\ backslash\n";
    r.inputs = digest_files({{"light", "prov_test.bin"}, {"stars", "prov_test.bin"}});
    r.outputs = {{"phot", "out dir/phot \"x\".csv", 12345, std::string(64, 'a')}};
    write_provenance("prov_test/reduce.json", r);
    const ProvenanceRecord b = read_provenance("prov_test/reduce.json");
    bool same = b.stage == r.stage && b.started_utc == r.started_utc && b.finished_utc == r.finished_utc &&
                b.project_file == r.project_file && b.settings_toml == r.settings_toml && b.build.version == r.build.version &&
                b.inputs.size() == 2 && b.outputs.size() == 1;
    for (size_t i = 0; same && i < 2; ++i)
        same = b.inputs[i].path == r.inputs[i].path && b.inputs[i].sha256 == r.inputs[i].sha256 && b.inputs[i].size == 1000000;
    same = same && b.outputs[0].path == r.outputs[0].path && b.outputs[0].size == 12345;
    check(same, "a provenance record survives a write/read round trip");
    std::filesystem::remove("prov_test.bin");
    std::filesystem::remove_all("prov_test");
    std::printf("%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
