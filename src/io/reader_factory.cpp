// SPDX-License-Identifier: GPL-3.0-or-later
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "transit/image_reader.h"

namespace transit {

#ifdef TRANSIT_WITH_CFITSIO
std::unique_ptr<ImageReader> make_cfitsio_reader();
#endif
#ifdef TRANSIT_WITH_XISFCONV
std::unique_ptr<ImageReader> make_xisfconv_reader();
#endif

std::unique_ptr<ImageReader> make_reader(std::string_view backend) {
#ifdef TRANSIT_WITH_CFITSIO
    if (backend == "cfitsio") return make_cfitsio_reader();
#endif
#ifdef TRANSIT_WITH_XISFCONV
    if (backend == "xisfconv") return make_xisfconv_reader();
#endif
    throw std::invalid_argument("image backend not available in this build: " + std::string(backend));
}

#ifndef TRANSIT_WITH_CFITSIO
void write_fits_float(const std::string&, const Image&, const std::vector<HeaderCard>&) {
    throw std::runtime_error("FITS writing needs the CFITSIO backend");
}
#endif

#ifdef TRANSIT_WITH_CFITSIO
bool cfitsio_reentrant();
#endif

bool backend_thread_safe(std::string_view backend) {
#ifdef TRANSIT_WITH_CFITSIO
    if (backend == "cfitsio") return cfitsio_reentrant();   // build-time option of the CFITSIO library
#endif
    if (backend == "xisfconv") return true;                 // no global state; one context per reader
    return false;
}

std::vector<std::string> available_backends() {
    std::vector<std::string> v;
#ifdef TRANSIT_WITH_CFITSIO
    v.push_back("cfitsio");
#endif
#ifdef TRANSIT_WITH_XISFCONV
    v.push_back("xisfconv");
#endif
    return v;
}

// Type of a FITS-formatted value string, using the same categories as CFITSIO's fits_get_keytype.
HeaderCard::Type classify_fits_value(const std::string& v) {
    size_t a = v.find_first_not_of(' ');
    if (a == std::string::npos) return HeaderCard::Type::Undefined;
    if (v[a] == '\'') return HeaderCard::Type::String;
    std::string t = v.substr(a);
    while (!t.empty() && t.back() == ' ') t.pop_back();
    if (t == "T" || t == "F") return HeaderCard::Type::Bool;
    if (t.find_first_of(".EeDd") != std::string::npos) return HeaderCard::Type::Float;
    return HeaderCard::Type::Int;
}

static std::string unquote(const std::string& raw) {
    auto a = raw.find('\''), b = raw.rfind('\'');
    if (a == std::string::npos || b <= a) return raw;
    std::string s = raw.substr(a + 1, b - a - 1), out;
    for (size_t i = 0; i < s.size(); ++i) {           // '' -> '
        out += s[i];
        if (s[i] == '\'' && i + 1 < s.size() && s[i + 1] == '\'') ++i;
    }
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out;
}

std::optional<std::string> card_string(const std::vector<HeaderCard>& h, std::string_view key) {
    for (const auto& c : h)
        if (c.key == key) return c.type == HeaderCard::Type::String ? unquote(c.raw_value) : c.raw_value;
    return std::nullopt;
}

std::optional<double> card_double(const std::vector<HeaderCard>& h, std::string_view key) {
    auto s = card_string(h, key);
    if (!s) return std::nullopt;
    char* end = nullptr;
    const double v = std::strtod(s->c_str(), &end);
    if (end == s->c_str()) return std::nullopt;
    return v;
}

}  // namespace transit
